use nalgebra::*;
//pub mod examples;

use crate::openmesh::{OpenMesh, Tag, Mat4, Vertex, vertex};
use crate::rule::{Rule, RuleStep, Child};
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
            Child {
                rule: Rule::Recurse(curve_horn_thing_rule),
                xf: id,
                vmap: vec![0,1,2,3],
            },
            Child {
                rule: Rule::Recurse(curve_horn_thing_rule),
                xf: flip180,
                vmap: vec![3,2,1,0],
            },
        ],
    }
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
            Child {
                rule: Rule::Recurse(curve_horn_thing_rule),
                xf: m,
                vmap: vec![0,1,2,3],
            },
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

    let gen_rulestep = |rot: &Mat4| -> Child {
        let m: Mat4 = rot *
            Matrix4::new_scaling(0.5) *
            geometry::Translation3::new(6.0, 0.0, 0.0).to_homogeneous();
        Child {
            rule: Rule::Recurse(cube_thing_rule),
            xf: m,
            vmap: vec![],
        }
    };

    RuleStep {
        geom: mesh,
        final_geom: prim::empty_mesh(),
        children: turns.iter().map(gen_rulestep).collect(),
    }
}

// Conversion from Python & automata_scratch
fn ram_horn_start() -> RuleStep {
    let opening_xform = |i| {
        let r = std::f32::consts::FRAC_PI_2 * i;
        ((geometry::Rotation3::from_axis_angle(
            &nalgebra::Vector3::z_axis(), r).to_homogeneous()) *
         geometry::Translation3::new(0.25, 0.25, 1.0).to_homogeneous() *
         Matrix4::new_scaling(0.5) *
         geometry::Translation3::new(0.0, 0.0, -1.0).to_homogeneous())
    };
    RuleStep {
        geom: OpenMesh {
            verts: vec![
                // 'Top' vertices:
                vertex(-0.5, -0.5, 1.0),  //  0 (above 9)
                vertex(-0.5,  0.5, 1.0),  //  1 (above 10)
                vertex( 0.5,  0.5, 1.0),  //  2 (above 11)
                vertex( 0.5, -0.5, 1.0),  //  3 (above 12)
                // Top edge midpoints:
                vertex(-0.5,  0.0, 1.0),  //  4 (connects 0-1)
                vertex( 0.0,  0.5, 1.0),  //  5 (connects 1-2)
                vertex( 0.5,  0.0, 1.0),  //  6 (connects 2-3)
                vertex( 0.0, -0.5, 1.0),  //  7 (connects 3-0)
                // Top middle:
                vertex( 0.0,  0.0, 1.0),  //  8
                // 'Bottom' vertices:
                vertex(-0.5, -0.5, 0.0),  //  9
                vertex(-0.5,  0.5, 0.0),  // 10
                vertex( 0.5,  0.5, 0.0),  // 11
                vertex( 0.5, -0.5, 0.0),  // 12
            ],
            faces: vec![
                // bottom face:
                Tag::Body(9), Tag::Body(10), Tag::Body(11),
                Tag::Body(9), Tag::Body(11), Tag::Body(12),
                // two faces straddling edge from vertex 0:
                Tag::Body(9), Tag::Body(0), Tag::Body(4),
                Tag::Body(9), Tag::Body(7), Tag::Body(0),
                // two faces straddling edge from vertex 1:
                Tag::Body(10), Tag::Body(1), Tag::Body(5),
                Tag::Body(10), Tag::Body(4), Tag::Body(1),
                // two faces straddling edge from vertex 2:
                Tag::Body(11), Tag::Body(2), Tag::Body(6),
                Tag::Body(11), Tag::Body(5), Tag::Body(2),
                // two faces straddling edge from vertex 3:
                Tag::Body(12), Tag::Body(3), Tag::Body(7),
                Tag::Body(12), Tag::Body(6), Tag::Body(3),
                // four faces from edge (0,1), (1,2), (2,3), (3,0):
                Tag::Body(9), Tag::Body(4), Tag::Body(10),
                Tag::Body(10), Tag::Body(5), Tag::Body(11),
                Tag::Body(11), Tag::Body(6), Tag::Body(12),
                Tag::Body(12), Tag::Body(7), Tag::Body(9),
            ],
        },
        final_geom: prim::empty_mesh(),
        children: vec![
            Child {
                rule: Rule::Recurse(ram_horn),
                xf: opening_xform(0.0),
                vmap: vec![5,2,6,8],
            },
            Child {
                rule: Rule::Recurse(ram_horn),
                xf: opening_xform(1.0),
                vmap: vec![4,1,5,8],
            },
            Child {
                rule: Rule::Recurse(ram_horn),
                xf: opening_xform(2.0),
                vmap: vec![7,0,4,8],
            },
            Child {
                rule: Rule::Recurse(ram_horn),
                xf: opening_xform(3.0),
                vmap: vec![6,3,7,8],
            },
            // TODO: These vertex mappings appear to be right.
            // Explain *why* they are right.
        ],
    }
}

fn ram_horn() -> RuleStep {
    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Mat4 = geometry::Translation3::new(0.0, 0.0, 0.8).to_homogeneous() *
        geometry::Rotation3::from_axis_angle(&v, 0.3).to_homogeneous() *
        Matrix4::new_scaling(0.9);
    let seed = vec![
        vertex(-0.5, -0.5, 1.0),
        vertex(-0.5,  0.5, 1.0),
        vertex( 0.5,  0.5, 1.0),
        vertex( 0.5, -0.5, 1.0),
    ];
    let next = seed.iter().map(|v| incr * v).collect();
    let geom = OpenMesh {
        verts: next,
        faces: vec![
            Tag::Body(1), Tag::Parent(0), Tag::Body(0),
            Tag::Parent(1), Tag::Parent(0), Tag::Body(1),
            Tag::Body(2), Tag::Parent(1), Tag::Body(1),
            Tag::Parent(2), Tag::Parent(1), Tag::Body(2),
            Tag::Body(3), Tag::Parent(2), Tag::Body(2),
            Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
            Tag::Body(0), Tag::Parent(3), Tag::Body(3),
            Tag::Parent(0), Tag::Parent(3), Tag::Body(0),
        ],
    };
    let final_geom = OpenMesh {
        verts: vec![],
        faces: vec![
            Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
            Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
        ],
    };
    RuleStep {
        geom: geom,
        final_geom: final_geom,
        children: vec![
            Child {
                rule: Rule::Recurse(ram_horn),
                xf: incr,
                vmap: vec![0,1,2,3],
            },
        ],
    }
}

/*
fn ram_horn_branch() -> RuleStep {
    
}
 */

// Meant to be a copy of twist_from_gen from Python & automata_scratch
pub fn twist_start() -> RuleStep {
    //let ang=0.1;
    let dz=0.2;
    let dx0=2.0;
    let dy=0.1;
    let count=4;
    // TODO: Factor these out (see twist)

    let seed = vec![
        vertex(-0.5,  0.0, -0.5),
        vertex(-0.5,  0.0,  0.5),
        vertex( 0.5,  0.0,  0.5),
        vertex( 0.5,  0.0, -0.5),        
    ];
    // TODO: Subdivide method
    
    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;
    let y = &Vector3::y_axis();
    let xform = |i| {
        (geometry::Rotation3::from_axis_angle(y, qtr * (i as f32)).to_homogeneous() *
         geometry::Translation3::new(dx0, 0.0, 0.0).to_homogeneous())
    };
    
    // First generate 'count' children, each one shifted/rotated
    // differently:
    let children: Vec<Child> = (0..count).map(|i| {
        let xf = xform(i);
        Child {
            rule: Rule::Recurse(twist),
            xf: xf,
            vmap: (4*i..4*(i+count)).collect(), // N.B.
        }
    }).collect();

    // Use byproducts of this to make 'count' copies of 'seed' with
    // this same transform:
    let mut verts = vec![];
    for child in &children {
        verts.extend(seed.iter().map(|v| child.xf * v));
    }
    
    RuleStep {
        geom: OpenMesh {
            verts: verts,
            faces: vec![],
            // TODO: Close these initial faces off
        },
        final_geom: prim::empty_mesh(),
        children: children,
    }
}

pub fn twist() -> RuleStep {
    let ang=0.1;
    let dz=0.2;
    let dx0=2.0;
    let dy=0.1;
    let count=4;
    // TODO: Factor these out (see twist_start)

    let y = &Vector3::y_axis();
    let incr = geometry::Rotation3::from_axis_angle(y, ang).to_homogeneous() *
        geometry::Translation3::new(0.0, dy, 0.0).to_homogeneous();
    
    let seed = vec![
        vertex(-0.5,  0.0, -0.5),
        vertex(-0.5,  0.0,  0.5),
        vertex( 0.5,  0.0,  0.5),
        vertex( 0.5,  0.0, -0.5),        
        // TODO: Likewise factor these out
    ].iter().map(|v| incr * v).collect();
    
    RuleStep {
        geom: OpenMesh {
            verts: seed,
            faces: vec![
                Tag::Body(1), Tag::Parent(0), Tag::Body(0),
                Tag::Parent(1), Tag::Parent(0), Tag::Body(1),
                Tag::Body(2), Tag::Parent(1), Tag::Body(1),
                Tag::Parent(2), Tag::Parent(1), Tag::Body(2),
                Tag::Body(3), Tag::Parent(2), Tag::Body(2),
                Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
                Tag::Body(0), Tag::Parent(3), Tag::Body(3),
                Tag::Parent(0), Tag::Parent(3), Tag::Body(0),
            ],
        },
        final_geom: prim::empty_mesh(), // TODO: Close properly
        children: vec![
            Child {
                rule: Rule::Recurse(twist),
                xf: incr,
                vmap: vec![0,1,2,3],
            },
        ],
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
    // this can't work on its own because the resultant OpenMesh still
    // has parent references:
    //run_test(Rule::Recurse(curve_horn_thing_rule), 100, "curve_horn_thing");
    run_test(Rule::Recurse(curve_horn_start), 100, "curve_horn2");
    run_test(Rule::Recurse(ram_horn_start), 200, "ram_horn");
    run_test(Rule::Recurse(twist_start), 200, "twist");
}
