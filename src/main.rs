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

fn curve_horn_thing_rule(v: Vec<Mesh>) -> Vec<RuleStep> {
 
    // TODO:
    // Accept 'input' seeds in v.
    // Draw everything relative to v.
    /*
    let gen_geom = |seed: &Mesh| -> RuleStep {
    }*/

    // TODO:
    // Of what is drawn here, for any seed that is passed to a future
    // rule, compute an xform which moves that seed to a sane
    // coordinate system.  Pass the transformed seed and xform
    // forward.

    for seed in v {
        for halfedge_id in seed.edge_iter() {
            let (v1, v2) = seed.edge_vertices(halfedge_id);
            println!("Half-edge {}, verts {} & {}: {}",
                     halfedge_id,
                     v1, v2,
                     if seed.is_edge_on_boundary(halfedge_id) {
                         "boundary"
                     } else {
                         "non-boundary"
                     }
                     );
        }
    }

    panic!("Not implemented");
    return vec![];
}

fn points_to_xform(v0: Point3<f64>, v1: Point3<f64>, v2: Point3<f64>) -> Mat4 {
    let x:  Vec3 = (v1 - v0).normalize();
    let z:  Vec3 = x.cross(v2 - v0).normalize();
    let y:  Vec3 = z.cross(x);

    let _m: Mat4 = Matrix4::from_cols(
        x.extend(0.0),       // new X
        y.extend(0.0),       // new Y
        z.extend(0.0),       // new Z
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

    // Is there a better way to do this?
    let _empty_mesh = MeshBuilder::new().with_indices(vec![]).with_positions(vec![]).build().unwrap();
    
    // Compute the extreme coordinates which defines the axis aligned bounding box..
    let (_min_coordinates, _max_coordinates) = mesh.extreme_coordinates();
    
    // .. or construct an actual mesh representing the axis aligned bounding box
    let _aabb = mesh.axis_aligned_bounding_box();

    let xform = points_to_xform(
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(4.0, 1.0, 0.0),
        Point3::new(2.0, 4.0, 0.0),
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

    let r2 = Rule::Recurse(curve_horn_thing_rule);
    println!("Running rules...");
    // Seed:
    let indices: Vec<u32> = vec![0, 1, 2,  2, 1, 3];
    let positions: Vec<f64> = vec![0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  1.0, 1.0, 0.0];
    let seed = MeshBuilder::new().with_indices(indices).with_positions(positions).build().unwrap();
    let (mesh, nodes) = rule_to_mesh(&r2, vec![seed], max_iters);
    println!("Collected {} nodes, produced {} faces, {} vertices",
             nodes, mesh.no_faces(), mesh.no_vertices());
    println!("Writing OBJ...");
    std::fs::write("curve_horn_thing.obj", mesh.parse_as_obj()).unwrap();
    // TODO: Can I make the seed geometry part of the rule itself?
}
