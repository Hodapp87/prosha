//use std::io;
use tri_mesh::prelude::*;

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
        //let seed2 = mesh.clone();

        // Put all vertices together:
        let mut pos = seed.positions_buffer();
        pos.append(&mut mesh.positions_buffer());
        let num_verts = seed.no_vertices();
        let mut indices: Vec<u32> = vec![0; 2 * num_verts * 3];
        let nv2 = u32::try_from(num_verts).unwrap();

        // TODO: don't I need to check if these are boundary edges or
        // something?  I have to be traversing in some kind of sane
        // order!  That is, this loop needs to be traversing the outer
        // boundary one vertex at a time.
        
        for i in 0..num_verts {
            let j = u32::try_from(i).unwrap();
            let k = u32::try_from(num_verts + i).unwrap();
            // First triangle:
            indices[6*i + 0] = j;
            indices[6*i + 1] = (j + 1) % nv2;
            indices[6*i + 2] = k;
            // Second triangle:
            indices[6*i + 3] = j; // k;
            indices[6*i + 4] = (j + 1) % nv2; // j + 1;
            indices[6*i + 5] = k; // (k + 1) % (2 * nv2);
        }
        // TODO: Above needs some clarity
        // (also, to be fixed)

        let joined = match MeshBuilder::new().with_positions(pos).with_indices(indices).build() {
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

    // TODO: If I transform a clone of the seed then this will leave
    // the vertices in the same order and preserve things like
    // boundaries.  I would then need to copy positions_buffer for
    // both the seed and the transformed seed (to give the vertices),
    // and use this in conjuction with my own indices in order to do
    // the zig-zag connection

    // DEBUG
    /*
    for seed in &v {

        let boundary = seed.edge_iter().filter(|e| seed.is_edge_on_boundary(*e));

        for halfedge_id in boundary {
            let (v1, v2) = seed.edge_vertices(halfedge_id);
            println!("Boundary half-edge {}, verts {} & {}",
                     halfedge_id, v1, v2);
        }

        // So, I have my boundary edges in no particular order.
        // I suppose I could use a Walker to give some order.
        //
        // But: how do I then connect these vertices up?  I can easily
        // get boundary vertices and make something based on those but
        // that then limits me to the 'cage' thing I was stuck with
        // prior.
    }
    */

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
    // TODO: Something is wrong here
    let mb = MeshBound::new(&seed).unwrap();
    for bound_edge in mb {
        println!("Boundary edge: {}", bound_edge);
    }

    let (mesh, nodes) = rule_to_mesh(&r2, vec![seed], 75);
    println!("Collected {} nodes, produced {} faces, {} vertices",
             nodes, mesh.no_faces(), mesh.no_vertices());
    println!("Writing OBJ...");
    std::fs::write("curve_horn_thing.obj", mesh.parse_as_obj()).unwrap();
    // TODO: Can I make the seed geometry part of the rule itself?
}
