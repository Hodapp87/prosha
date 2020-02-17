use nalgebra::*;
//pub mod examples;

use crate::openmesh::{OpenMesh, Tag, Mat4, Vertex, vertex};
use crate::rule::{Rule, RuleStep};
use crate::prim;

fn curve_horn_start() -> RuleStep {
    let id = nalgebra::geometry::Transform3::identity().to_homogeneous();
    let flip180 = nalgebra::geometry::Rotation3::from_axis_angle(
        &nalgebra::Vector3::y_axis(),
        std::f32::consts::PI).to_homogeneous();
    RuleStep {
        geom: OpenMesh {
            verts: vec![],
            faces: vec![
                Tag::Exit(1, 0), Tag::Exit(1, 2), Tag::Exit(0, 1),
                Tag::Exit(1, 2), Tag::Exit(0, 3), Tag::Exit(0, 1), 
                Tag::Exit(0, 0), Tag::Exit(0, 2), Tag::Exit(1, 1), 
                Tag::Exit(0, 2), Tag::Exit(1, 3), Tag::Exit(1, 1), 
                Tag::Exit(0, 3), Tag::Exit(1, 2), Tag::Exit(0, 2), 
                Tag::Exit(1, 2), Tag::Exit(1, 3), Tag::Exit(0, 2), 
                Tag::Exit(1, 0), Tag::Exit(0, 1), Tag::Exit(0, 0),
                Tag::Exit(1, 1), Tag::Exit(1, 0), Tag::Exit(0, 0), 
                // The above is connecting group 0 to group 1,
                // straight across + with diagonal - but with group 1
                // being flipped 180, so we remap vertices (0,1,2,3)
                // to (1,0,3,2) and then flip winding order.
            ],
            exit_groups: vec![4, 4],
        },
        final_geom: prim::empty_mesh(),
        children: vec![
            (Rule::Recurse(curve_horn_thing_rule), id), // exit group 0
            (Rule::Recurse(curve_horn_thing_rule), flip180), // exit group 1
        ],
    }
    // TODO: The starting vertices above are duplicated because I
    // don't have any way for an exit vertex to stand in for multiple
    // child vertices that happen to share the same location.  I don't
    // yet know a good way around this, so I am duplicating vertices.
}

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
            // The below is just connecting two groups of 4 vertices
            // each, straight across and then to the next.  Note that
            // since 'verts' doesn't go in a circle, it will look a
            // little strange.
            Tag::Body(1), Tag::Exit(0, 3), Tag::Exit(0, 1),
            Tag::Body(1), Tag::Body(3), Tag::Exit(0, 3),
            Tag::Exit(0, 0), Tag::Body(2), Tag::Body(0),
            Tag::Exit(0, 0), Tag::Exit(0, 2), Tag::Body(2),
            Tag::Body(2), Tag::Exit(0, 3), Tag::Body(3),
            Tag::Body(2), Tag::Exit(0, 2), Tag::Exit(0, 3),
            Tag::Body(0), Tag::Body(1), Tag::Exit(0, 1),
            Tag::Body(0), Tag::Exit(0, 1), Tag::Exit(0, 0),
            // TODO: I should really generate these, not hard-code them.
        ],
        exit_groups: vec![4],
    };

    // TODO: This could be made slightly nicer by taking it to a peak
    // instead of just flattening it in XY, but this is a pretty minor
    // change.
    let final_geom = OpenMesh {
        verts: final_verts,
        faces: vec![
            Tag::Body(0), Tag::Body(1), Tag::Body(3),
            Tag::Body(0), Tag::Body(3), Tag::Body(2),
        ],
        exit_groups: vec![],
    };
    
    RuleStep{
        geom: geom,
        final_geom: final_geom,
        children: vec![
            (Rule::Recurse(curve_horn_thing_rule), m), // exit group 0
        ],
    }
}

fn cube_thing_rule() -> RuleStep {

    let mesh = prim::cube();

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
        final_geom: prim::empty_mesh(), // no exit groups
        children: turns.iter().map(gen_rulestep).collect(),
    }
}

pub fn main() {

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
    run_test(Rule::Recurse(curve_horn_start), 100, "curve_horn2");
}
