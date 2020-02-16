use nalgebra::*;
use std::fs::OpenOptions;
use std::io;

/// A type for custom mesh vertices. Initialize with [vertex][self::vertex].
pub type Vertex = Vector4<f32>;
pub type Mat4 = Matrix4<f32>;

/// Initializes a vertex:
pub fn vertex(x: f32, y: f32, z: f32) -> Vertex {
    Vertex::new(x, y, z, 1.0)
}

#[derive(Clone, Debug)]
enum Tag {
    Body(usize),
    Exit(usize, usize),  // (group, vertex)
}

#[derive(Clone, Debug)]
struct OpenMesh {
    // Vertices (in homogeneous coordinates).
    verts: Vec<Vertex>,
    // Triangles, taken as every 3 values, treated each as indices
    // into 'verts':
    faces: Vec<Tag>,
    exit_groups: Vec<usize>,
}

impl OpenMesh {
    
    fn transform(&self, xfm: Mat4) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            // TODO: Is the above faster if I pack vectors into a
            // bigger matrix, and transform that?
            faces: self.faces.clone(), // TODO: Use Rc?
            exit_groups: self.exit_groups.clone(),
        }
    }

    fn write_stl_file(&self, fname: &str) -> io::Result<()> {
        let mut file = OpenOptions::new().write(true).create(true).truncate(true).open(fname)?;
        self.write_stl(&mut file)
    }
    
    fn write_stl<W: std::io::Write>(&self, writer: &mut W) -> io::Result<()> {

        // Every group of 3 indices in self.faces is one triangle, so
        // pre-allocate in the format stl_io wants:
        let num_faces = self.faces.len() / 3;
        let mut triangles = vec![stl_io::Triangle {
            normal: [0.0; 3],
            vertices: [[0.0; 3]; 3],
        }; num_faces];

        let get_vert = |j| {
            match self.faces[j] {
                Tag::Body(n) => self.verts[n].xyz(),
                Tag::Exit(_, _) => panic!("Cannot write_stl() if mesh has exit groups!"),
            }
        };
        // TODO: Handle this behavior
        
        // Turn every face into an stl_io::Triangle:
        for i in 0..num_faces {
            let v0 = get_vert(3*i + 0);
            let v1 = get_vert(3*i + 1);
            let v2 = get_vert(3*i + 2);
            
            let normal = (v1-v0).cross(&(v2-v0));

            triangles[i].normal.copy_from_slice(&normal.as_slice());
            triangles[i].vertices[0].copy_from_slice(v0.as_slice());
            triangles[i].vertices[1].copy_from_slice(v1.as_slice());
            triangles[i].vertices[2].copy_from_slice(v2.as_slice());
            // TODO: Is there a cleaner way to do the above?
        }

        // I could also solve this with something like
        // https://doc.rust-lang.org/std/primitive.slice.html#method.chunks_exact
        // however I don't know what performance difference may be.

        stl_io::write_stl(writer, triangles.iter())
    }

    fn connect(&self, others: &Vec<OpenMesh>) -> OpenMesh {

        //println!("DEBUG: connect(), self has {} exit groups, others have {:?}",
        //         self.exit_groups.len(), others.iter().map(|o| o.exit_groups.len()).collect::<Vec<usize>>());
        //println!("DEBUG: connect(), self: verts.len()={} faces.len()={} max face={}", self.verts.len(), self.faces.len(), self.faces.iter().map(|f| match f { Tag::Body(n) => n, Tag::Exit(_,n) => n }).max().unwrap());
        
        // Copy body vertices & faces:
        let mut verts: Vec<Vertex> = self.verts.clone();
        let mut faces = self.faces.clone();

        let mut exit_groups: Vec<usize> = vec![];
        
        let mut body_offset = self.verts.len();
        let mut exit_offset = 0;
        let mut offsets: Vec<usize> = vec![0; others.len()];
        for (i,other) in others.iter().enumerate() {

            //let max_ = other.faces.iter().map(|f| match f { Tag::Body(n) => n, Tag::Exit(_,n) => n }).max().unwrap_or(&0);
            //println!("DEBUG: connect(), other[{}]: verts.len()={} faces.len()={} max face={}", i, other.verts.len(), other.faces.len(), max_);
            //println!("DEBUG: start body_offset={}", body_offset);
            //println!("DEBUG: start exit_offset={}", exit_offset);
            
            // Append body vertices & exit vertices directly:
            verts.append(&mut other.verts.clone());
            
            // Append faces, shifting each kind by respective offset:
            faces.extend(other.faces.iter().map(|t| {
                match t {
                    Tag::Body(n) => Tag::Body(n + body_offset),
                    Tag::Exit(g, n) => Tag::Exit(g + exit_groups.len(), n + exit_offset),
                }
            }));
            if i < self.exit_groups.len() {
                exit_offset += self.exit_groups[i];
            }
            exit_groups.append(&mut other.exit_groups.clone());

            offsets[i] = body_offset;
            // Increase offsets by the number of elements we appended:
            body_offset += other.verts.len();

            //println!("DEBUG: end body_offset={}", body_offset);
            //println!("DEBUG: end exit_offset={}", exit_offset);
        }

        //println!("DEBUG: offsets={:?}", offsets);

        // All of the Exit face indices from 'self' need to be
        // modified to refer to Body vertices of something in
        // 'others':
        //println!("DEBUG: initial faces={:?}", faces);
        for i in 0..faces.len() {
            match faces[i] {
                Tag::Exit(g, n) => {
                    faces[i] = Tag::Body(n + offsets[g]);
                },
                _ => { },
            };
        }
        //println!("DEBUG: final faces={:?}", faces);

        let m = OpenMesh {
            verts: verts,
            faces: faces,
            exit_groups: exit_groups,
        };

        // TODO: Why is this still ending up with Exit faces despite my loop above?
        //println!("DEBUG: Returning mesh with verts.len()={} faces.len()={} max face={}", m.verts.len(), m.faces.len(), m.faces.iter().map(|f| match f { Tag::Body(n) => n, Tag::Exit(_,n) => n }).max().unwrap());
        //println!("Returning: {:?}", m);
        return m;
    }
}

// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?
enum Rule {
    // Produce geometry, and possibly recurse further:
    Recurse(fn () -> RuleStep),
    // Stop recursing here:
    EmptyRule,
}
// TODO: Rename rules?
// TODO: It may be possible to have just a 'static' rule that requires
// no function call.

struct RuleStep {
    // The geometry generated by this rule on its own (not by any of
    // the child rules).
    geom: OpenMesh,

    // The "final" geometry, used only if recursion must be stopped.
    // This should be in the same coordinate space as 'geom', and
    // properly close any exit groups that it may have (and have no
    // exit groups of its own).
    final_geom: OpenMesh,

    // Child rules, paired with the transform that will be applied to
    // all of their geometry
    children: Vec<(Rule, Mat4)>,
}

impl Rule {

    // TODO: Do I want to make 'geom' shared somehow, maybe with Rc? I
    // could end up having a lot of identical geometry that need not be
    // duplicated until it is transformed into the global space.
    //
    // This might produce bigger gains if I rewrite rule_to_mesh so that
    // rather than repeatedly transforming meshes, it stacks
    // transformations and then applies them all at once.

    fn to_mesh(&self, iters_left: u32) -> (OpenMesh, u32) {

        let mut nodes: u32 = 1;

        if iters_left <= 0 {
            match self {
                Rule::Recurse(f) => {
                    let rs: RuleStep = f();
                    return (rs.final_geom, 1);
                }
                Rule::EmptyRule => {
                    return (empty_mesh(), nodes);
                }
            }
        }

        match self {
            Rule::Recurse(f) => {
                let rs: RuleStep = f();

                // Get sub-geometry (from child rules) and transform it:
                let subgeom: Vec<(OpenMesh, Mat4, u32)> = rs.children.iter().map(|(subrule, subxform)| {
                    let (m,n) = subrule.to_mesh(iters_left - 1);
                    (m, *subxform, n)
                }).collect();

                // Tally up node count:
                subgeom.iter().for_each(|(_,_,n)| nodes += n);

                let g: Vec<OpenMesh> = subgeom.iter().map(|(m,x,_)| m.transform(*x)).collect();

                // Connect geometry from this rule (not child rules):
                return (rs.geom.connect(&g), nodes);
            }
            Rule::EmptyRule => {
                return (empty_mesh(), nodes);
            }
        }
    }
}

// is there a better way to do this?
fn empty_mesh() -> OpenMesh {
    OpenMesh {
        verts: vec![],
        faces: vec![],
        exit_groups: vec![],
    }
}

fn cube() -> OpenMesh {
    OpenMesh {
        verts: vec![
            vertex(0.0, 0.0, 0.0),
            vertex(1.0, 0.0, 0.0),
            vertex(0.0, 1.0, 0.0),
            vertex(1.0, 1.0, 0.0),
            vertex(0.0, 0.0, 1.0),
            vertex(1.0, 0.0, 1.0),
            vertex(0.0, 1.0, 1.0),
            vertex(1.0, 1.0, 1.0),
        ],
        faces: vec![
            Tag::Body(0), Tag::Body(3), Tag::Body(1),
            Tag::Body(0), Tag::Body(2), Tag::Body(3),
            Tag::Body(1), Tag::Body(7), Tag::Body(5),
            Tag::Body(1), Tag::Body(3), Tag::Body(7),
            Tag::Body(5), Tag::Body(6), Tag::Body(4),
            Tag::Body(5), Tag::Body(7), Tag::Body(6),
            Tag::Body(4), Tag::Body(2), Tag::Body(0),
            Tag::Body(4), Tag::Body(6), Tag::Body(2),
            Tag::Body(2), Tag::Body(7), Tag::Body(3),
            Tag::Body(2), Tag::Body(6), Tag::Body(7),
            Tag::Body(0), Tag::Body(1), Tag::Body(5),
            Tag::Body(0), Tag::Body(5), Tag::Body(4),
        ],
        exit_groups: vec![],
    }.transform(geometry::Translation3::new(-0.5, -0.5, -0.5).to_homogeneous())
}

/*
fn curve_horn_start() -> RuleStep {
    // Seed is a square in XY, sidelength 1, centered at (0,0,0):
    let seed = {
        let m = OpenMesh {
            verts: vec![
                vertex(0.0, 0.0, 0.0),
                vertex(1.0, 0.0, 0.0),
                vertex(1.0, 1.0, 0.0),
                vertex(0.0, 1.0, 0.0),
            ],
            faces: vec![
                0, 1, 2,
                0, 2, 3,
            ],
            idxs_entrance: vec![0],
            idxs_exit: vec![0],
            idxs_body: (0, 0),
        };
        let xform = nalgebra::geometry::Translation3::new(-0.5, -0.5, 0.0).to_homogeneous();
        m.transform(xform)
    };
    vec![
        // Since neither of the other two rules *start* with geometry:
        RuleStep { geom: seed.clone(),
                   rule: Box::new(Rule::EmptyRule),
                   xform: nalgebra::geometry::Transform3::identity().to_homogeneous(),
        },
        // Recurse in both directions:
        RuleStep { geom: seed.clone(),
                   rule: Box::new(Rule::Recurse(curve_horn_thing_rule)),
                   xform: nalgebra::geometry::Transform3::identity().to_homogeneous(),
        },
        RuleStep { geom: seed.clone(),
                   rule: Box::new(Rule::Recurse(curve_horn_thing_rule)),
                   xform: nalgebra::geometry::Rotation3::from_axis_angle(
                       &nalgebra::Vector3::y_axis(),
                       std::f32::consts::FRAC_PI_2).to_homogeneous(),
        },
    ]
}

//use std::convert::TryFrom;
*/

fn curve_horn_thing_rule() -> RuleStep {

    let y = &Vector3::y_axis();
    
    let m: Mat4 = geometry::Rotation3::from_axis_angle(y, 0.1).to_homogeneous() *
        Matrix4::new_scaling(0.95) *
        geometry::Translation3::new(0.0, 0.0, 0.2).to_homogeneous();

    let verts = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(0.5, -0.5, 0.0),
        vertex(-0.5, 0.5, 0.0),
        vertex(0.5, 0.5, 0.0),
    ];
    let final_verts: Vec<Vertex> = verts.iter().map(|v| m * v).collect();
    
    let geom = OpenMesh {
        verts: verts,
        faces: vec![
            // Endcaps purposely left off for now.
            // TODO: I should really generate these, not hard-code them.
            Tag::Body(1), Tag::Exit(0, 3), Tag::Exit(0, 1),
            Tag::Body(1), Tag::Body(3), Tag::Exit(0, 3),
            Tag::Exit(0, 0), Tag::Body(2), Tag::Body(0),
            Tag::Exit(0, 0), Tag::Exit(0, 2), Tag::Body(2),
            Tag::Body(2), Tag::Exit(0, 3), Tag::Body(3),
            Tag::Body(2), Tag::Exit(0, 2), Tag::Exit(0, 3),
            Tag::Body(0), Tag::Body(1), Tag::Exit(0, 1),
            Tag::Body(0), Tag::Exit(0, 1), Tag::Exit(0, 0),
        ],
        exit_groups: vec![4],
    };

    let final_geom = OpenMesh {
        verts: final_verts,
        faces: vec![
            Tag::Body(0), Tag::Body(3), Tag::Body(1),
            Tag::Body(0), Tag::Body(2), Tag::Body(3),
        ],
        exit_groups: vec![],
    };
    
    RuleStep{
        geom: geom,
        final_geom: final_geom,
        children: vec![
            (Rule::Recurse(curve_horn_thing_rule), m),
        ],
    }

    /*
        // We need 3 indices per face, 2 faces per (boundary) vertex:
        let num_verts = seed.no_vertices();
        let mut idxs: Vec<u32> = vec![0; 2 * num_verts * 3];
        for i in 0..num_verts {
            let a1: u32 = i                                   as _;
            let a2: u32 = ((i + 1) % num_verts)               as _;
            let b1: u32 = (i + num_verts)                     as _;
            let b2: u32 = (((i + 1) % num_verts) + num_verts) as _;
            // Connect vertices into faces with a zig-zag pattern
            // (mind the winding order).  First face:
            
            idxs[6*i + 0] = a1;
            idxs[6*i + 1] = a2;
            idxs[6*i + 2] = b1;
            //println!("connect vert {}, face 1: ({}, {}, {}) = {}, {}, {}", i, a1, a2, b1, vert2str(a1), vert2str(a2), vert2str(b1));
            // Second face:
            idxs[6*i + 3] = b1;
            idxs[6*i + 4] = a2;
            idxs[6*i + 5] = b2;
            //println!("connect vert {}, face 2: ({}, {}, {}) = {}, {}, {}", i, b1, a2, b2, vert2str(b1), vert2str(a2), vert2str(b2));
        }
    */
}

fn cube_thing_rule() -> RuleStep {

    let mesh = cube();

    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;

    let y = &Vector3::y_axis();
    let z = &Vector3::z_axis();
    
    // Each element of this turns to a branch for the recursion:
    let turns: Vec<Mat4> = vec![
        geometry::Transform3::identity().to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr * 2.0).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr * 3.0).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(z, qtr).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(z, -qtr).to_homogeneous(),
    ];

    let gen_rulestep = |rot: &Mat4| -> (Rule, Mat4) {
        let m: Mat4 = rot *
            Matrix4::new_scaling(0.5) *
            geometry::Translation3::new(6.0, 0.0, 0.0).to_homogeneous();
        (Rule::Recurse(cube_thing_rule), m)
    };

    RuleStep {
        geom: mesh,
        final_geom: empty_mesh(), // no exit groups
        children: turns.iter().map(gen_rulestep).collect(),
    }
}

fn main() {

    let run_test = |r: Rule, iters, name| {
        println!("Running {}...", name);
        let (mesh, nodes) = r.to_mesh(iters);
        println!("Merged {} nodes", nodes);
        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        mesh.write_stl_file(&fname).unwrap();
    };

    run_test(Rule::Recurse(cube_thing_rule), 3, "cube_thing");
    run_test(Rule::Recurse(curve_horn_thing_rule), 100, "curve_horn_thing");
}
