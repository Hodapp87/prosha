//use std::io;
use tri_mesh::prelude::*;
use nalgebra::base::dimension::{U1, U4};
//use nalgebra::Matrix4;

/// A type for custom mesh vertices. Initialize with [vertex][self::vertex].
pub type Vertex = nalgebra::Matrix<f32, U4, U1, nalgebra::ArrayStorage<f32, U4, U1>>;
// TODO: Why am I not just using nalgebra's RowVector type here?

/// Initializes a vertex for a custom mesh.
pub fn vertex(x: f32, y: f32, z: f32) -> Vertex {
    Vertex::new(x, y, z, 1.0)
}

struct OpenMesh {
    // Vertices (in homogeneous coordinates).  These must be in a
    // specific order: 'Entrance' loops, then 'body' vertices, then
    // 'exit' loops.
    verts: Vec<Vertex>,
    // Triangles, taken as every 3 values, treated each as indices
    // into 'verts':
    faces: Vec<usize>,
    // A list of indices into verts, telling the index at which each
    // 'entrance' vertex loop begins.  The loop implicitly ends where the
    // next one begins, or if it is the last one, at idxs_body._1.
    // Thus, this has one element per vertex loop, and must go in
    // ascending order.
    idxs_entrance: Vec<usize>,
    // The same as idxs_entrance, but for 'exit' vertex loops.  The final
    // loop is taken as ending at the end of the list.
    idxs_exit:  Vec<usize>,
    // The start and end (non-inclusive) of the 'body' vertices -
    // those that are neither an entrance nor an exit loop.
    idxs_body: (usize, usize),
}
// TODO: What is proper for an index, u32 or usize?

impl OpenMesh {
    
    fn transform(&self, xfm: nalgebra::Matrix4<f32>) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            faces: self.faces.clone(), // TODO: Use Rc?
            idxs_entrance: self.idxs_entrance.clone(), // TODO: Use Rc?
            idxs_exit: self.idxs_exit.clone(), // TODO: Use Rc?
            idxs_body: self.idxs_body.clone(), // TODO: Use Rc?
        }
    }

    fn connect(&self, others: &Vec<OpenMesh>) -> OpenMesh {

        let mut v = self.verts.clone();
        let mut f = self.faces.clone();

        for other in others {
            // We are offsetting all vertices in 'other' by everything
            // else in 'v', so we need to account for this when we
            // copy 'faces' (which has vector indices):
            let offset = v.len();
            v.extend(other.verts.iter());
            f.extend(other.faces.iter().map(|f| *f + offset));
        }
        
        // - Append entrance loop & body vertices of all.
        // - Append all faces - but add offset to those from 'other'.
        // - Connect up so that each of self's exit loops is an
        // entrance loop from one of 'other'

        return OpenMesh {
            verts: v,
            faces: f,
            idxs_entrance: self.idxs_entrance.clone(), // TODO
            idxs_exit: self.idxs_exit.clone(), // TODO
            idxs_body: self.idxs_body.clone(), // TODO
        };
    }
}

// TODO: Does OpenMesh subsume both 'geom' and 'seeds' in RuleStep?
// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?

enum Rule {
    // Recurse further.  Input is "seeds" that further geometry should
    // *replace*.  Generated geometry must have the same outer
    // boundary as the seeds, and be in the same coordinate space as
    // the input.
    Recurse(fn (Vec<Mesh>) -> Vec<RuleStep>),
    // Stop recursing here:
    EmptyRule,
}
// TODO: Rename rules?

struct RuleStep {
    // The 'final' geometry generated at this step.
    geom: Mesh,
    // The 'seed' geometry from this step.  If recursion stops
    // (whether because rule is EmptyRule or because recursion depth
    // has been hit), this will be transformed with 'xform' and
    // appended with 'geom'. If recursion continues, this geometry is
    // passed as the input to the next rule.  (TODO: rule_to_mesh
    // needs to do the 'recursion stops' part.)
    //
    // This is in the coordinate space that 'rule' should run in -
    // thus, if it is transformed with 'xform', it will be in the same
    // coordinate space as 'geom'.
    seeds: Vec<Mesh>,
    // The next rule to run.  If EmptyRule, then stop here (and
    // 'xform' is irrelevant).
    rule: Box<Rule>,
    // The transformation which puts 'seeds' and any geometry from
    // 'rule' (if applicable) into the same coordinate space as
    // 'geom'.
    xform: Mat4,
}

// is there a better way to do this?
fn empty_mesh() -> Mesh {
    MeshBuilder::new().with_indices(vec![]).with_positions(vec![]).build().unwrap()
}

fn curve_horn_start(_v: Vec<Mesh>) -> Vec<RuleStep> {
    // Seed is a square in XY, sidelength 1, centered at (0,0,0):
    let seed = {
        let indices: Vec<u32> = vec![0, 1, 2,  0, 2, 3];
        let positions: Vec<f64> = vec![0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  1.0, 1.0, 0.0,  0.0, 1.0, 0.0];
        let mut s = MeshBuilder::new().with_indices(indices).with_positions(positions).build().unwrap();
        s.apply_transformation(Matrix4::from_translation(vec3(-0.5, -0.5, 0.0)));
        s
    };
    vec![
        // Since neither of the other two rules *start* with geometry:
        RuleStep { geom: seed.clone(),
                   rule: Box::new(Rule::EmptyRule),
                   xform: Matrix4::identity(),
                   seeds: vec![]
        },
        // Recurse in both directions:
        RuleStep { geom: empty_mesh(),
                   rule: Box::new(Rule::Recurse(curve_horn_thing_rule)),
                   xform: Matrix4::identity(),
                   seeds: vec![seed.clone()],
        },
        RuleStep { geom: empty_mesh(),
                   rule: Box::new(Rule::Recurse(curve_horn_thing_rule)),
                   xform: Matrix4::from_angle_y(Rad::turn_div_2()),
                   seeds: vec![seed.clone()],
        },
    ]
}

use std::convert::TryFrom;

fn curve_horn_thing_rule(v: Vec<Mesh>) -> Vec<RuleStep> {
    
    let gen_geom = |seed: &Mesh| -> RuleStep {
        let mut mesh = seed.clone();

        let m: Mat4 = Matrix4::from_angle_y(Rad(0.1)) *
            Matrix4::from_scale(0.95) *
            Matrix4::from_translation(vec3(0.0, 0.0, 0.2));

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
        
        let joined = match MeshBuilder::new().
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

    let _m: Mat4 = Matrix4::from_cols(
        (xn*s).extend(0.0),   // new X
        (yn*s).extend(0.0),   // new Y
        (zn*s).extend(0.0),   // new Z
        v0.to_homogeneous(), // translation
    );
    return _m;
}

fn cube_thing_rule(_v: Vec<Mesh>) -> Vec<RuleStep> {

    let mesh = MeshBuilder::new().cube().build().unwrap();

    // Quarter-turn in radians:
    let qtr = Rad::turn_div_4();
    
    // Each element of this turns to a branch for the recursion:
    let turns: Vec<Mat4> = vec![
        Matrix4::identity(),
        Matrix4::from_angle_y(qtr),
        Matrix4::from_angle_y(qtr * 2.0),
        Matrix4::from_angle_y(qtr * 3.0),
        Matrix4::from_angle_z(qtr),
        Matrix4::from_angle_z(-qtr),
    ];

    let gen_rulestep = |rot: &Mat4| -> RuleStep {
        let m: Mat4 = rot *
            Matrix4::from_scale(0.5) *
            Matrix4::from_translation(vec3(6.0, 0.0, 0.0));
        let r = Rule::Recurse(cube_thing_rule);
        let mut m2 = mesh.clone();
        m2.apply_transformation(m);
        RuleStep { geom: m2, rule: Box::new(r), xform: m, seeds: vec![] }
    };
    // TODO: Why is 'mesh' present in each RuleStep?  This is just
    // duplicate geometry! Either 'm' applies to 'mesh' (and the
    // definition of RuleStep changes) - or 'mesh' needs to already be
    // transformed.

    turns.iter().map(gen_rulestep).collect()
}

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
        let (v1, v2) = self.m.edge_vertices(self.cur);
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

//fn mesh_boundary(m: &Mesh) -> Vec<tri_mesh::HalfEdgeID> {
//}

// TODO: Do I want to make 'geom' shared somehow, maybe with Rc? I
// could end up having a lot of identical geometry that need not be
// duplicated until it is transformed into the global space.
//
// This might produce bigger gains if I rewrite rule_to_mesh so that
// rather than repeatedly transforming meshes, it stacks
// transformations and then applies them all at once.

fn rule_to_mesh(rule: &Rule, seed: Vec<Mesh>, iters_left: u32) -> (Mesh, u32) {

    let mut mesh = MeshBuilder::new().with_indices(vec![]).with_positions(vec![]).build().unwrap();

    let mut nodes: u32 = 1;
    
    if iters_left <= 0 {
        return (mesh, nodes);
    }

    match rule {
        Rule::Recurse(func) => {
            for step in func(seed) {
                let subrule: Rule = *step.rule;
                let subxform: Mat4 = step.xform;
                let geom: Mesh = step.geom;

                mesh.append(&geom);
                
                let (mut submesh, subnodes) = rule_to_mesh(
                    &subrule, step.seeds, iters_left - 1);
                submesh.apply_transformation(subxform);
                nodes += subnodes;

                mesh.append(&submesh);
            }
        }
        Rule::EmptyRule => {
            // do nothing
        }
    }
    (mesh, nodes)
}

fn print_vector(v: &Vec4) -> String {
    return format!("{},{},{},{}", v.x, v.y, v.z, v.w);
}

fn print_matrix(m: &Mat4) {
    let mt = m.transpose();
    println!("[{}]\n[{}]\n[{}]\n[{}]",
             print_vector(&mt.x), print_vector(&mt.y),
             print_vector(&mt.z), print_vector(&mt.w));
}

fn main() {
    // Construct any mesh, this time, we will construct a simple icosahedron
    let mesh = MeshBuilder::new().icosahedron().build().unwrap();

    // Compute the extreme coordinates which defines the axis aligned bounding box..
    let (_min_coordinates, _max_coordinates) = mesh.extreme_coordinates();
    
    // .. or construct an actual mesh representing the axis aligned bounding box
    let _aabb = mesh.axis_aligned_bounding_box();

    let xform = points_to_xform(
        Point3::new(0.5, 0.5, 0.0),
        Point3::new(-0.5, 0.5, 0.0),
        Point3::new(2.0, -4.0, 0.0),
    );
    println!("points_to_xform:");
    print_matrix(&xform);
    
    // Export the bounding box to an obj file
    std::fs::write("foo.obj", mesh.parse_as_obj()).unwrap();

    let r = Rule::Recurse(cube_thing_rule);

    let max_iters = 2;
    println!("Running rules...");
    let (cubemesh, nodes) = rule_to_mesh(&r, vec![], max_iters);
    println!("Collected {} nodes, produced {} faces, {} vertices",
             nodes, cubemesh.no_faces(), cubemesh.no_vertices());
    println!("Writing OBJ...");
    std::fs::write("cubemesh.obj", cubemesh.parse_as_obj()).unwrap();

    let r2 = Rule::Recurse(curve_horn_start);
    println!("Running rules...");
    // Seed:
    let seed = {
        let indices: Vec<u32> = vec![0, 1, 2,  2, 1, 3];
        let positions: Vec<f64> = vec![0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  1.0, 1.0, 0.0];
        let mut s = MeshBuilder::new().with_indices(indices).with_positions(positions).build().unwrap();
        s.apply_transformation(Matrix4::from_translation(vec3(-0.5, -0.5, 0.0)));
        s
    };
    // TEMP (while I figure shit out)
    struct VID { val: usize }
    fn vertex_id_to_usize(v: VertexID) -> usize {
        let v: VID = unsafe { std::mem::transmute(v) };
        v.val
    }
    println!("DEBUG-------------------------------");
    let mb = MeshBound::new(&seed).unwrap();
    let pos = seed.positions_buffer();
    for bound_edge in mb {
        let (v1, v2) = seed.edge_vertices(bound_edge);
        let v1idx = vertex_id_to_usize(v1);
        let v2idx = vertex_id_to_usize(v2);
        
        println!("Boundary edge {}, vertices = {},{}, {:?}",
                 bound_edge, v1, v2, seed.edge_positions(bound_edge));
        println!("v1idx={} pos[...]=[{},{},{}], v2idx={}, pos[...]=[{},{},{}]",
                 v1idx, pos[3*v1idx], pos[3*v1idx+1], pos[3*v1idx+2],
                 v2idx, pos[3*v2idx], pos[3*v2idx+1], pos[3*v2idx+2]);
    }
    println!("DEBUG-------------------------------");

    let (mut mesh, nodes) = rule_to_mesh(&r2, vec![seed], 75);
    println!("Collected {} nodes, produced {} faces, {} vertices",
             nodes, mesh.no_faces(), mesh.no_vertices());
    println!("Trying to merge...");
    match mesh.merge_overlapping_primitives() {
        Err(e) => {
            println!("Couldn't merge overlapping primitives!");
            println!("Error: {:?}", e);
        }
        Ok(v) => {
            println!("Merged to {} faces, {} vertices",
                     mesh.no_faces(), mesh.no_vertices());
        }
    }
    println!("Writing OBJ...");
    std::fs::write("curve_horn_thing.obj", mesh.parse_as_obj()).unwrap();
    // TODO: Can I make the seed geometry part of the rule itself?
}
