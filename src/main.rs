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
struct OpenMesh {
    // Vertices (in homogeneous coordinates).  These must be in a
    // specific order: 'Entrance' loops, then 'body' vertices, then
    // 'exit' loops.
    verts: Vec<Vertex>,
    // Triangles, taken as every 3 values, treated each as indices
    // into 'verts':
    faces: Vec<usize>,
    // A list of indices into verts, telling the index at which each
    // 'entrance' vertex group begins.  The group implicitly ends
    // where the next one begins, or if it is the last group, at
    // idxs_body._1.  Thus, this has one element per vertex group, and
    // must go in ascending order.
    idxs_entrance: Vec<usize>,
    // The same as idxs_entrance, but for 'exit' vertex groups.  The
    // final loop is taken as ending at the end of the list.
    idxs_exit:  Vec<usize>,
    // The start and end (non-inclusive) of the 'body' vertices -
    // those that are neither an entrance nor an exit group.
    idxs_body: (usize, usize),
}
// TODO: Do I even use idxs_entrance?  Is it still valuable as a
// cross-check?

impl OpenMesh {
    
    fn transform(&self, xfm: Mat4) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            // TODO: Is the above faster if I pack vectors into a
            // bigger matrix?
            faces: self.faces.clone(), // TODO: Use Rc?
            idxs_entrance: self.idxs_entrance.clone(), // TODO: Use Rc?
            idxs_exit: self.idxs_exit.clone(), // TODO: Use Rc?
            idxs_body: self.idxs_body.clone(), // TODO: Use Rc?
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

        // Turn every face into an stl_io::Triangle:
        for i in 0..num_faces {
            let v0 = self.verts[self.faces[3*i + 0]].xyz();
            let v1 = self.verts[self.faces[3*i + 1]].xyz();
            let v2 = self.verts[self.faces[3*i + 2]].xyz();
            
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

    fn connect_single(&self, other: &OpenMesh) -> OpenMesh {

        // Imagine connecting two pieces of pipe together.  We are
        // fitting the exit of 'self' to the entrance of 'other' - and
        // producing a new piece of pipe which has the entrance of
        // 'self', but the exit of 'other'.
        
        let mut v: Vec<Vertex> = vec![vertex(0.0,0.0,0.0); self.idxs_body.1];
        // Start out by cloning just entrance & body vertices:
        v.copy_from_slice(&self.verts[0..self.idxs_body.1]);
        let mut f = self.faces.clone();

        // I already know what size v will be so I can pre-allocate
        // and then just clone_from_slice to the proper locations

        // We are offsetting all vertices in 'other' by everything
        // else in 'v', so we need to account for this when we copy
        // 'faces' (which has vector indices):
        let offset = self.idxs_body.1;
        f.extend(other.faces.iter().map(|f| *f + offset));
        v.extend(other.verts.iter());
        // The new exit groups are those in 'other', but likewise we
        // need to shift these indices:
        let idxs_exit = other.idxs_exit.iter().map(|f| *f + offset).collect();
        // Body vertices start in the same place, but end where the
        // body vertices in 'other' end (thus needing offset):
        let idxs_body = (self.idxs_body.0, other.idxs_body.1 + offset);

        OpenMesh {
            verts: v,
            faces: f,
            idxs_entrance: self.idxs_entrance.clone(),
            idxs_exit: idxs_exit,
            idxs_body: idxs_body,
        }
    }

    // Just assume this is broken
    fn connect(&self, others: &Vec<OpenMesh>) -> OpenMesh {

        if others.len() > 1 && self.idxs_exit.len() > 0 {
            panic!("connect() is implemented for only one mesh if exit groups are present")
        }

        if false {
            let mut v: Vec<Vertex> = vec![vertex(0.0,0.0,0.0); self.verts.len()];
            // Start out by cloning just entrance & body vertices:
            v.copy_from_slice(&self.verts[0..self.idxs_body.1]);
            let mut f = self.faces.clone();
            // TODO: Don't I need to patch up 'f'?  self.faces refers to
            // exit vertices which - if others.len() > 1 - need to be
            // manually patched up.  This patching up should consist
            // solely of an offset to all indices in a certain range.
            //
            // e.g. let idxs_exit be [e0, e1, e2, ... e_(n-1)]
            // indices in range [e0, e1-1] are for exit group 0.
            // indices in range [e1, e2-1] are for exit group 1.
            // indices in range [e2, e3-1] are for exit group 2, etc.
            //
            // exit group 0 requires no offset (we'll be putting entrance
            // group vertices of self.others[0] right over top of them).
            // 
            // exit group 1 requires an offset of the number of entrace &
            // body vertices of self.others[0] (because we have appended
            // this all)... with some additional adjustment maybe?  not
            // sure.
            //
            // exit group 2 requires an offset of the same for
            // self.others[0] and self.others[1].

            for other in others {
                // We are offsetting all vertices in 'other' by everything
                // else in 'v', so we need to account for this when we
                // copy 'faces' (which has vector indices):
                let offset = v.len();
                v.extend(other.verts[0..other.idxs_body.1].iter());
                f.extend(other.faces.iter().map(|f| *f + offset));
            }
            
            // - Connect up so that each of self's exit groups is an
            // entrance group from one of 'other'

            return OpenMesh {
                verts: v,
                faces: f,
                idxs_entrance: self.idxs_entrance.clone(),
                idxs_exit: self.idxs_exit.clone(), // TODO
                idxs_body: self.idxs_body.clone(), // TODO
            };
        }

        // This is wrong, but close enough for now;
        let mut mesh = self.clone();
        for other in others {
            mesh = mesh.connect_single(&other);
        }
        return mesh;
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
            idxs_entrance: vec![],
            idxs_exit: vec![],
            idxs_body: (0, 0),
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
            0, 3, 1,
            0, 2, 3,
            1, 7, 5,
            1, 3, 7,
            5, 6, 4,
            5, 7, 6,
            4, 2, 0,
            4, 6, 2,
            2, 7, 3,
            2, 6, 7,
            0, 1, 5,
            0, 5, 4,
        ],
        idxs_entrance: vec![],
        idxs_exit: vec![],
        idxs_body: (0, 8),
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

    let mut verts = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(0.5, -0.5, 0.0),
        vertex(-0.5, 0.5, 0.0),
        vertex(0.5, 0.5, 0.0),
    ];
    let mut v2: Vec<Vertex> = verts.iter().map(|v| m * v).collect();
    let final_verts: Vec<Vertex> = v2.clone();
    verts.append(&mut v2);
    
    let geom = OpenMesh {
        verts: verts,
        faces: vec![
            // Endcaps purposely left off for now.
            // TODO: I should really generate these, not hard-code them.
            1, 7, 5,
            1, 3, 7,
            4, 2, 0,
            4, 6, 2,
            2, 7, 3,
            2, 6, 7,
            0, 1, 5,
            0, 5, 4,
        ],
        idxs_entrance: vec![0],
        idxs_exit: vec![4],
        idxs_body: (4, 4),
    };

    let final_geom = OpenMesh {
        verts: final_verts,
        faces: vec![
            0, 3, 1,
            0, 2, 3,
        ],
        idxs_entrance: vec![0],
        idxs_exit: vec![],
        idxs_body: (4, 4),
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

    // Below is so far my only example that uses entrance/exit groups:
    println!("DEBUG-------------------------------");
    let m = OpenMesh {
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
            // End caps disabled for now to test connect_single
            // 0, 3, 1,
            // 0, 2, 3,
            1, 7, 5,
            1, 3, 7,
            // 5, 6, 4,
            // 5, 7, 6,
            4, 2, 0,
            4, 6, 2,
            2, 7, 3,
            2, 6, 7,
            0, 1, 5,
            0, 5, 4,
        ],
        idxs_entrance: vec![0],
        idxs_exit: vec![4],
        idxs_body: (4, 4),
    };

    let xform = geometry::Translation3::new(0.0, 0.0, 1.0).to_homogeneous();
    let m2 = m.transform(xform);
    let m3 = m.connect_single(&m2);
    let m4 = m3.connect_single(&m2.transform(xform));
    println!("m4 = {:?}", m4);

    m.write_stl_file("openmesh_cube.obj").unwrap();
    m2.write_stl_file("openmesh_cube2.obj").unwrap();
    m3.write_stl_file("openmesh_cube3.obj").unwrap();

    {
        let count = 10;
        let mut mesh = m.clone();
        let mut inc = m.clone();
        for _ in 0..count {
            inc = inc.transform(xform);
            mesh = mesh.connect_single(&inc);
        }
    }

    {
        let r = Rule::Recurse(cube_thing_rule);
        let max_iters = 4;
        println!("Running rules...");
        let (cubemesh, nodes) = r.to_mesh(max_iters);
        println!("Merged {} nodes", nodes);
        println!("Writing STL...");
        cubemesh.write_stl_file("cubemesh.stl").unwrap();
    }

    {
        let r = Rule::Recurse(curve_horn_thing_rule);
        let max_iters = 50;
        println!("Running rules...");
        let (cubemesh, nodes) = r.to_mesh(max_iters);
        //println!("cubemesh={:?}", cubemesh);
        println!("Merged {} nodes", nodes);
        println!("Writing STL...");
        cubemesh.write_stl_file("curve_horn_thing.stl").unwrap();
    }
}
