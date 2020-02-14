//use std::io;
use tri_mesh::prelude as tm;
use nalgebra::*;

/// A type for custom mesh vertices. Initialize with [vertex][self::vertex].
pub type Vertex = Vector4<f32>;

/// Initializes a vertex for a custom mesh.
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
// TODO: What is proper for an index, u32 or usize?

impl OpenMesh {
    
    fn transform(&self, xfm: Matrix4<f32>) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            faces: self.faces.clone(), // TODO: Use Rc?
            idxs_entrance: self.idxs_entrance.clone(), // TODO: Use Rc?
            idxs_exit: self.idxs_exit.clone(), // TODO: Use Rc?
            idxs_body: self.idxs_body.clone(), // TODO: Use Rc?
        }
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

    fn to_trimesh(&self) -> Result<tm::Mesh, tri_mesh::mesh_builder::Error> {
        let mut v: Vec<f64> = vec![0.0; self.verts.len() * 3];
        for (i, vert) in self.verts.iter().enumerate() {
            v[3*i] = vert[0].into();
            v[3*i+1] = vert[1].into();
            v[3*i+2] = vert[2].into();
        }
        let faces: Vec<u32> = self.faces.iter().map(|f| *f as _).collect();
        tm::MeshBuilder::new().with_indices(faces).with_positions(v).build()
    }

    // Just assume this is broken
    fn connect(&self, others: &Vec<OpenMesh>) -> OpenMesh {

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
}

// TODO: Does OpenMesh subsume both 'geom' and 'seeds' in RuleStep?
// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?

enum Rule {
    // Produce geometry, and possibly recurse further:
    Recurse(fn () -> Vec<RuleStep>),
    // Stop recursing here:
    EmptyRule,
}
// TODO: Rename rules?

struct RuleStep {
    // The geometry generated by this rule on its own - and none of
    // the child rules.
    geom: OpenMesh,
    
    // The next rule to run.  If EmptyRule, then stop here (and
    // 'xform' is irrelevant).
    rule: Box<Rule>,
    // The transformation to apply to geometry generated by 'rule' and
    // any child rules.
    xform: Matrix4<f32>,
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

/*
fn curve_horn_start() -> Vec<RuleStep> {
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

fn curve_horn_thing_rule() -> Vec<RuleStep> {
    
    let gen_geom = |seed: &Mesh| -> RuleStep {
        let mut mesh = seed.clone();

        let m: Mat4 = tm::Matrix4::from_angle_y(Rad(0.1)) *
            tm::Matrix4::from_scale(0.95) *
            tm::Matrix4::from_translation(vec3(0.0, 0.0, 0.2));

        let r = Rule::Recurse(curve_horn_thing_rule);
        mesh.apply_transformation(m);

        // TODO: Fix this horrible code below that is seemingly
        // correct, but shouldn't be run on every rule iteration!

        // Collect together all the vertices from the boundaries of
        // 'seed' and 'mesh':
        let edge2vert = |m: &Mesh, e: HalfEdgeID| {
            let v = m.vertex_position(m.edge_vertices(e).0);
            vec![v.x, v.y, v.z]
        };
        let i1 = MeshBound::new(&seed).unwrap().flat_map(|id| edge2vert(&seed, id));
        let i2 = MeshBound::new(&mesh).unwrap().flat_map(|id| edge2vert(&mesh, id));
        let verts: Vec<f64> = i1.chain(i2).collect();
        
        /*
        let vert2str = |idx: u32| {
            let i2: usize = idx as _;
            format!("({:.4},{:.4},{:.4})", verts[3*i2], verts[3*i2+1], verts[3*i2+2])
        };
        for i in 0..(seed.no_vertices() + mesh.no_vertices()) {
            println!("vert {}: {}", i, vert2str(i as _))
        }
        */
        
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
        // TODO: Something is *still* not quite right there.  I think
        // that I cannot use MeshBuilder this way and then append
        // meshes - it just leads to disconnected geometry
        
        let joined = match tm::MeshBuilder::new().
            with_positions(verts).
            with_indices(idxs).
            build()
        {
            Ok(m) => m,
            Err(error) => {
                panic!("Error building mesh: {:?}", error)
            },
        };
        
        RuleStep { geom: joined, rule: Box::new(r), xform: m, seeds: vec![seed.clone()] }
    };
    // Since 'mesh' is computed directly by applying 'm' to 'seed',
    // trivially, we follow the requirement in a RuleStep that
    // applying 'xform' to 'seeds' puts it into the same space as
    // 'geom'.

    v.iter().map(gen_geom).collect()
}

// Assume v0, v1, and v2 are non-collinear points.  This tries to
// produce a transform which treats v0 as the origin of a new
// coordinate system, the line from v0 to v1 as the new X axis, the Y
// axis perpendicular to this along the plane that (v0,v1,v2) forms,
// and the Z axis the normal of this same plane.
//
// Scale is taken into account (to the extent that the length of
// (v1-v0) is taken as distance 1 in the new coordinate system).
fn points_to_xform(v0: Point3<f64>, v1: Point3<f64>, v2: Point3<f64>) -> Mat4 {
    let x:  Vec3 = v1 - v0;
    let xn: Vec3 = x.normalize();
    let zn:  Vec3 = x.cross(v2 - v0).normalize();
    let yn:  Vec3 = zn.cross(xn);
    let s = x.magnitude();

    let _m: Mat4 = tm::Matrix4::from_cols(
        (xn*s).extend(0.0),   // new X
        (yn*s).extend(0.0),   // new Y
        (zn*s).extend(0.0),   // new Z
        v0.to_homogeneous(), // translation
    );
    return _m;
}
*/

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

fn cube_thing_rule() -> Vec<RuleStep> {

    let mesh = cube();

    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;

    let y = &Vector3::y_axis();
    let z = &Vector3::z_axis();
    
    // Each element of this turns to a branch for the recursion:
    let turns: Vec<Matrix4<f32>> = vec![
        geometry::Transform3::identity().to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr * 2.0).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(y, qtr * 3.0).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(z, qtr).to_homogeneous(),
        geometry::Rotation3::from_axis_angle(z, -qtr).to_homogeneous(),
    ];

    let gen_rulestep = |rot: &Matrix4<f32>| -> RuleStep {
        let m: Matrix4<f32> = rot *
            Matrix4::new_scaling(0.5) *
            geometry::Translation3::new(6.0, 0.0, 0.0).to_homogeneous();
        let r = Rule::Recurse(cube_thing_rule);

        let m2 = mesh.transform(m);
        RuleStep { geom: m2, rule: Box::new(r), xform: m }
    };

    turns.iter().map(gen_rulestep).collect()
}

// Have I any need of this after making OpenMesh?
/*
struct MeshBound<'a> {
    m: &'a Mesh,
    start: HalfEdgeID,
    cur: HalfEdgeID,
    done: bool,
}

impl<'a> MeshBound<'a> {
    fn new(m: &'a Mesh) -> Option<MeshBound> {
        for halfedge_id in m.edge_iter() {
            if m.is_edge_on_boundary(halfedge_id) {
                return Some(MeshBound {
                    m: m,
                    start: halfedge_id,
                    cur: halfedge_id,
                    done: false,
                });
            }
        }
        // TODO: Maybe just return an iterator that returns None
        // immediately if this mesh has no boundary?
        return None;
    }
}

impl<'a> Iterator for MeshBound<'a> {
    type Item = HalfEdgeID;

    fn next(&mut self) -> Option<Self::Item> {

        if self.done {
            return None;
        }

        // Start from our current half-edge:
        let (v1, _) = self.m.edge_vertices(self.cur);
        // Pick a vertex and walk around incident half-edges:
        for halfedge_id in self.m.vertex_halfedge_iter(v1) {

            // Avoid twin half-edge, which returns where we started:
            let w = self.m.walker_from_halfedge(halfedge_id);
            if w.twin_id().map_or(false, |twin| twin == self.cur) {
                continue;
            }
            // TODO: is there a quicker way to get the twin?

            // If this incident half-edge is a boundary, follow it:
            if self.m.is_edge_on_boundary(halfedge_id) {
                
                self.cur = halfedge_id;
                if self.start == self.cur {
                    // We have returned back to start:
                    self.done = true;
                }
                //println!("DEBUG:   MeshBound: edge {} is {:?}", halfedge_id, self.m.edge_positions(halfedge_id));
                return Some(halfedge_id);
            }
        }
        return None;
    }
    
}
*/

//fn mesh_boundary(m: &Mesh) -> Vec<tri_mesh::HalfEdgeID> {
//}

// TODO: Do I want to make 'geom' shared somehow, maybe with Rc? I
// could end up having a lot of identical geometry that need not be
// duplicated until it is transformed into the global space.
//
// This might produce bigger gains if I rewrite rule_to_mesh so that
// rather than repeatedly transforming meshes, it stacks
// transformations and then applies them all at once.

fn rule_to_mesh(rule: &Rule, iters_left: u32) -> (OpenMesh, u32) {

    let mut mesh = empty_mesh();

    let mut nodes: u32 = 1;
    
    if iters_left <= 0 {
        return (mesh, nodes);
    }

    match rule {
        Rule::Recurse(func) => {
            for step in func() {
                let subrule: Rule = *step.rule;
                let subxform: Matrix4<f32> = step.xform;
                let geom: OpenMesh = step.geom;

                mesh = mesh.connect_single(&geom);
                
                let (mut submesh, subnodes) = rule_to_mesh(
                    &subrule, iters_left - 1);

                submesh = submesh.transform(subxform);

                nodes += subnodes;

                mesh = mesh.connect_single(&submesh);
            }
        }
        Rule::EmptyRule => {
            // do nothing
        }
    }
    (mesh, nodes)
}

fn main() {

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
    
    let try_save = |m: &OpenMesh, fname: &str| {
        let m_trimesh = m.to_trimesh().unwrap();
        std::fs::write(fname, m_trimesh.parse_as_obj()).unwrap();                                    
    };

    try_save(&m, "openmesh_cube.obj");
    try_save(&m2, "openmesh_cube2.obj");
    try_save(&m3, "openmesh_cube3.obj");

    {
        let count = 10;
        let mut mesh = m.clone();
        let mut inc = m.clone();
        for _ in 0..count {
            inc = inc.transform(xform);
            mesh = mesh.connect_single(&inc);
        }
        println!("mesh = {:?}", mesh);
        try_save(&mesh, "openmesh_cube_several.obj");
    }

    let r = Rule::Recurse(cube_thing_rule);

    let max_iters = 2;
    println!("Running rules...");
    let (cubemesh_, nodes) = rule_to_mesh(&r, max_iters);
    let cubemesh = cubemesh_.to_trimesh().unwrap();
    println!("Collected {} nodes, produced {} faces, {} vertices",
             nodes, cubemesh.no_faces(), cubemesh.no_vertices());
    println!("Writing OBJ...");
    std::fs::write("cubemesh.obj", cubemesh.parse_as_obj()).unwrap();

    /*
    let r2 = Rule::Recurse(curve_horn_start);
    println!("Running rules...");
    // Seed:
    let seed = {
        let indices: Vec<u32> = vec![0, 1, 2,  2, 1, 3];
        let positions: Vec<f64> = vec![0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  1.0, 1.0, 0.0];
        let mut s = tm::MeshBuilder::new().with_indices(indices).with_positions(positions).build().unwrap();
        s.apply_transformation(tm::Matrix4::from_translation(vec3(-0.5, -0.5, 0.0)));
        s
    };
    */
    // TEMP (while I figure shit out)
    println!("DEBUG-------------------------------");

}
