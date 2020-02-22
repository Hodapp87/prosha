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
            verts: vec![
                vertex(-0.5, -0.5, 0.0),
                vertex(-0.5,  0.5, 0.0),
                vertex( 0.5,  0.5, 0.0),
                vertex( 0.5, -0.5, 0.0),
            ],
            faces: vec![],
        },
        final_geom: prim::empty_mesh(),
        children: vec![
            (Rule::Recurse(curve_horn_thing_rule), id),
            (Rule::Recurse(curve_horn_thing_rule), flip180),
        ],
    }
    // TODO: Fix the consequences of the 180 flip
}

fn curve_horn_thing_rule() -> RuleStep {

    let y = &Vector3::y_axis();
    
    let m: Mat4 = geometry::Rotation3::from_axis_angle(y, 0.1).to_homogeneous() *
        Matrix4::new_scaling(0.95) *
        geometry::Translation3::new(0.0, 0.0, 0.2).to_homogeneous();

    let verts = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let next_verts: Vec<Vertex> = verts.iter().map(|v| m * v).collect();
    
    let geom = OpenMesh {
        verts: next_verts.clone(),
        faces: vec![
            // The below is just connecting two groups of 4 vertices
            // each, straight across and then to the next.
            Tag::Body(1), Tag::Parent(0), Tag::Body(0),
            Tag::Parent(1), Tag::Parent(0), Tag::Body(1),
            Tag::Body(2), Tag::Parent(1), Tag::Body(1),
            Tag::Parent(2), Tag::Parent(1), Tag::Body(2),
            Tag::Body(3), Tag::Parent(2), Tag::Body(2),
            Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
            Tag::Body(0), Tag::Parent(3), Tag::Body(3),
            Tag::Parent(0), Tag::Parent(3), Tag::Body(0),
            // TODO: I should really generate these, not hard-code them.
        ],
    };

    // TODO: This could be made slightly nicer by taking it to a peak
    // instead of just flattening it in XY, but this is a pretty minor
    // change.
    let final_geom = OpenMesh {
        verts: vec![],
        faces: vec![
            Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
            Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
        ],
    };
    
    RuleStep{
        geom: geom,
        final_geom: final_geom,
        children: vec![
            (Rule::Recurse(curve_horn_thing_rule), m),
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
        final_geom: prim::empty_mesh(),
        children: turns.iter().map(gen_rulestep).collect(),
    }
}

// Conversion from Python & automata_scratch
fn ram_horn_start() -> RuleStep {
    let id = nalgebra::geometry::Transform3::identity().to_homogeneous();
    let flip180 = nalgebra::geometry::Rotation3::from_axis_angle(
        &nalgebra::Vector3::y_axis(),
        std::f32::consts::PI).to_homogeneous();
    RuleStep {
        geom: OpenMesh {
            verts: vec![
                // 'Bottom' vertices:
                vertex(-0.5, -0.5, -0.5), // 0
                vertex(-0.5,  0.5, -0.5), // 1
                vertex( 0.5,  0.5, -0.5), // 2
                vertex( 0.5, -0.5, -0.5), // 3
                // 'Top' vertices:
                vertex(-0.5, -0.5, 0.5),  // 4 (above 0)
                vertex(-0.5,  0.5, 0.5),  // 5 (above 1)
                vertex( 0.5,  0.5, 0.5),  // 6 (above 2)
                vertex( 0.5, -0.5, 0.5),  // 7 (above 3)
                // Top edge midpoints:
                vertex(-0.5,  0.0, 0.5),  // 8 (connects 4-5)
                vertex( 0.0,  0.5, 0.5),  // 9 (connects 5-6)
                vertex( 0.5,  0.0, 0.5),  // 10 (connects 6-7)
                vertex( 0.0, -0.5, 0.5),  // 11 (connects 7-4)
                // Top middle:
                vertex( 0.0,  0.0, 0.5),
            ],
            faces: vec![
                // bottom face:
                Tag::Body(0), Tag::Body(1), Tag::Body(2),
                Tag::Body(0), Tag::Body(2), Tag::Body(3),
                // two faces straddling edge from vertex 0:
                Tag::Body(0), Tag::Body(4), Tag::Body(8),
                Tag::Body(0), Tag::Body(11), Tag::Body(4),
                // two faces straddling edge from vertex 1:
                Tag::Body(1), Tag::Body(5), Tag::Body(9),
                Tag::Body(1), Tag::Body(8), Tag::Body(5),
                // two faces straddling edge from vertex 2:
                Tag::Body(2), Tag::Body(6), Tag::Body(10),
                Tag::Body(2), Tag::Body(9), Tag::Body(6),
                // two faces straddling edge from vertex 3:
                Tag::Body(3), Tag::Body(7), Tag::Body(11),
                Tag::Body(3), Tag::Body(10), Tag::Body(7),
                // four faces from edge (0,1), (1,2), (2,3), (3,0):
                Tag::Body(0), Tag::Body(8), Tag::Body(1),
                Tag::Body(1), Tag::Body(9), Tag::Body(2),
                Tag::Body(2), Tag::Body(10), Tag::Body(3),
                Tag::Body(3), Tag::Body(11), Tag::Body(0),
            ],
        },
        final_geom: prim::empty_mesh(),
        children: vec![
            (Rule::Recurse(ram_horn), id), // exit group 0
            (Rule::Recurse(ram_horn), id), // exit group 1
            (Rule::Recurse(ram_horn), id), // exit group 2
            (Rule::Recurse(ram_horn), id), // exit group 3
        ],
    }
}

fn ram_horn() -> RuleStep {
    // Incomplete
    RuleStep {
        geom: OpenMesh {
            verts: vec![],
            faces: vec![],
        },
        final_geom: prim::empty_mesh(),
        children: vec![],
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
    //run_test(Rule::Recurse(curve_horn_thing_rule), 100, "curve_horn_thing");
    run_test(Rule::Recurse(curve_horn_start), 100, "curve_horn2");
    run_test(Rule::Recurse(ram_horn_start), 10, "ram_horn");
}
