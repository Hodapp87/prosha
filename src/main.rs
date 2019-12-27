//use std::io;
use tri_mesh::prelude::*;

enum Rule {
    // Recurse further:
    Recurse(fn (Vec<Mesh>) -> Vec<RuleStep>),
    // Stop recursing here:
    EmptyRule,
}
// TODO: Rename rules?

struct RuleStep {
    // The geometry generated at this step
    geom: Mesh,
    // The next rule to run on this geometry
    rule: Box<Rule>,
    // The transformation to apply to any results of 'rule' (if
    // applicable)
    xform: Mat4,
}

fn test_rule(_v: Vec<Mesh>) -> Vec<RuleStep> {

    let r = Rule::Recurse(test_rule);
    
    let mesh = MeshBuilder::new().cube().build().unwrap();
    //let mtx = Matrix4::identity();
    let m: Mat4 = 
        Matrix4::from_translation(vec3(1.0, 0.0, 0.0)) *
        Matrix4::from_angle_z(Rad(0.1)) *
        Matrix4::from_scale(0.9);
    
    let ret: Vec<RuleStep> = vec![
        RuleStep { geom: mesh, rule: Box::new(r), xform: m },
    ];

    ret
}

fn rule_to_mesh(rule: &Rule, xform: Mat4, iter_num: u32) -> Mesh {

    let max_iters: u32 = 50;
    let mut mesh = MeshBuilder::new().with_indices(vec![]).with_positions(vec![]).build().unwrap();

    if iter_num >= max_iters {
        return mesh;
    }
    
    match rule {
        Rule::Recurse(func) => {
            for step in func(vec![]) {
                let subrule: Rule = *step.rule;
                let subxform: Mat4 = step.xform;
                let geom: Mesh = step.geom;

                mesh.append(&geom);
                
                let mut submesh: Mesh = rule_to_mesh(&subrule, subxform, iter_num + 1);
                submesh.apply_transformation(xform);

                mesh.append(&submesh);
            }
            mesh
        }
        Rule::EmptyRule => {
            mesh
        }
    }
}

// This isn't kosher:
//type Rule = fn (Vec<Mesh>) -> (Mesh, Vec<(Mesh, Box<Rule>)>);

fn mesh_builder_example() -> Result<Mesh, tri_mesh::mesh_builder::Error> {
    let indices: Vec<u32> = vec![0, 1, 2,
                                 0, 2, 3,
                                 0, 3, 1];
    let positions: Vec<f64> = vec![0.0, 0.0, 0.0,
                                   1.0, 0.0, -0.5,
                                   -1.0, 0.0, -0.5,
                                   0.0, 0.0, 1.0];
    let mesh = MeshBuilder::new().
        with_indices(indices).
        with_positions(positions).
        build()?;

    assert_eq!(mesh.no_faces(), 3);
    assert_eq!(mesh.no_vertices(), 4);

    Ok(mesh)
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
    
    // Export the bounding box to an obj file
    std::fs::write("foo.obj", mesh.parse_as_obj()).unwrap();

    // Try some vector stuff:

    let m: Mat4 = Matrix4::from_translation(vec3(5.0, 0.0, 0.0));
    println!("translation: ");
    print_matrix(&m);
    /*print_vector(&m.x);
    print_vector(&m.y);
    print_vector(&m.z);
    print_vector(&m.w);*/

    let m2: Mat4 = Matrix4::from_scale(2.0);
    println!("scale: ");
    print_matrix(&m2);
    let m3 = m * m2;
    println!("translation * scale: ");
    print_matrix(&m3);
    let m4 = m2 * m;
    println!("scale * translation: ");
    print_matrix(&m4);

    let mut mesh2 = mesh_builder_example().unwrap();
    std::fs::write("foo2.obj", mesh2.parse_as_obj()).unwrap();

    mesh2.apply_transformation(m);
    
    mesh2.append(&mesh);

    std::fs::write("foo3.obj", mesh2.parse_as_obj()).unwrap();

    let r = Rule::Recurse(test_rule);

    match r {
        Rule::Recurse(f) => {
            let _v = f(vec![]);
        }
        Rule::EmptyRule => {
            println!("Empty");
        }
    }

    let cubemesh: Mesh = rule_to_mesh(&r, Matrix4::identity(), 0);
    std::fs::write("cubemesh.obj", cubemesh.parse_as_obj()).unwrap();
}
