use std::rc::Rc;
use nalgebra::*;
//pub mod examples;

use crate::openmesh::{OpenMesh, Tag};
use crate::xform::{Transform, vertex};
use crate::rule::{Rule, RuleFn, RuleEval, Child};
use crate::prim;
use crate::util;

use std::time::Instant;

fn cube_thing() -> Rule<()> {

    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;
    //let x = &Vector3::x_axis();
    let y = &Vector3::y_axis();
    let z = &Vector3::z_axis();
    
    // Each element of this turns to a branch for the recursion:
    let id = Transform::new();
    let turns: Vec<Transform> = vec![
        id.clone(),
        id.rotate(y, qtr),
        id.rotate(y, qtr * 2.0),
        id.rotate(y, qtr * 3.0),
        id.rotate(z, qtr),
        id.rotate(z, -qtr),
    ];

    let rec = move |self_: Rc<Rule<()>>| -> RuleEval<()> {

        let xforms = turns.iter().map(|xf| xf.scale(0.5).translate(6.0, 0.0, 0.0));
        RuleEval {
            geom: Rc::new(prim::cube()),
            final_geom: Rc::new(prim::empty_mesh()),
            children: xforms.map(move |xf| Child {
                rule: self_.clone(),
                xf: xf,
                vmap: vec![],
            }).collect(),
        }
    };
    
    Rule { eval: Rc::new(rec), ctxt: () }
}

// Meant to be a copy of twist_from_gen from Python & automata_scratch
fn twist(f: f32, subdiv: usize) -> Rule<()> {
    // TODO: Clean this code up.  It was a very naive conversion from
    // the non-closure version.
    let xf = Transform::new().rotate(&Vector3::x_axis(), -0.7);
    let seed = {
        let s = vec![vertex(-0.5,  0.0, -0.5),
                     vertex( 0.5,  0.0, -0.5),
                     vertex( 0.5,  0.0,  0.5),
                     vertex(-0.5,  0.0,  0.5)];
        util::subdivide_cycle(&xf.transform(&s), subdiv)
    };
    let n = seed.len();
    let dx0: f32 = 1.5;
    let dy: f32 = 0.1/f;
    let ang: f32 = 0.05/f;
    let count: usize = 4;
    
    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;
    let y = Vector3::y_axis();

    let incr_inner = Transform::new().translate(-dx0, 0.0, 0.0).rotate(&y, ang).translate(dx0, dy, 0.0);
    let incr_outer = Transform::new().translate(-dx0*2.0, 0.0, 0.0).rotate(&y, ang/2.0).translate(dx0*2.0, dy, 0.0);

    let seed2 = seed.clone();
    // TODO: Why do I need the above?
    let recur = move |incr: Transform| -> RuleFn<()> {

        let seed_next = incr.transform(&seed2);

        let geom = Rc::new(util::zigzag_to_parent(seed_next.clone(), n));
        // TODO: Cleanliness fix - why not just make these return meshes?
        let (vc, faces) = util::connect_convex(&seed_next, true);
        let final_geom = Rc::new(OpenMesh {
            verts: vec![vc],
            faces: faces,
        });
        
        let c = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
            RuleEval {
                geom: geom.clone(),
                final_geom: final_geom.clone(),
                children: vec![
                    Child {
                        rule: self_.clone(),
                        xf: incr,
                        vmap: (0..n).collect(),
                    },
                ],
            }
        };
        Rc::new(c)
    };
    // TODO: Can a macro do anything to clean up some of the
    // repetition with HOFs & closures?
    
    let start = move |_| -> RuleEval<()> {
        
        let child = |incr, dx, i, ang0, div| -> (OpenMesh, Child<()>) {
            let xform = Transform::new().
                rotate(&y, ang0 + (qtr / div * (i as f32))).
                translate(dx, 0.0, 0.0);
            
            let c = Child {
                rule: Rc::new(Rule { eval: (recur.clone())(incr), ctxt: () }),
                // TODO: Cleanliness fix - can macros clean up above?
                xf: xform,
                vmap: (0..(n+1)).collect(),
                // N.B. n+1, not n. the +1 is for the centroid below.
            };
            let mut vs = xform.transform(&seed);
            // and in the process, generate faces for these seeds:
            let (centroid, f) = util::connect_convex(&vs, false);
            vs.push(centroid);
            (OpenMesh { verts: vs, faces: f }, c)
        };
        
        // Generate 'count' children, shifted/rotated differently:
        let inner = (0..count).map(|i| child(incr_inner, dx0, i, 0.0, 1.0));
        let outer = (0..count).map(|i| child(incr_outer, dx0*2.0, i, qtr/2.0, 2.0));

        RuleEval::from_pairs(inner.chain(outer), prim::empty_mesh())
    };
    
    Rule { eval: Rc::new(start), ctxt: () }
}

fn ramhorn() -> Rule<()> {

    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 0.8).
        rotate(&v, 0.3).
        scale(0.9);

    let seed = vec![
        vertex(-0.5, -0.5, 1.0),
        vertex(-0.5,  0.5, 1.0),
        vertex( 0.5,  0.5, 1.0),
        vertex( 0.5, -0.5, 1.0),
    ];
    let next = incr.transform(&seed);
    let geom = Rc::new(OpenMesh {
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
    });
    let final_geom = Rc::new(OpenMesh {
        verts: vec![],
        faces: vec![
            Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
            Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
        ],
    });
    
    let recur = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
        RuleEval {
            geom: geom.clone(),
            final_geom: final_geom.clone(),
            children: vec![
                Child {
                    rule: self_.clone(),
                    xf: incr,
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    };
    
    let opening_xform = |i| {
        let r = std::f32::consts::FRAC_PI_2 * i;
        Transform::new().
            rotate(&nalgebra::Vector3::z_axis(), r).
            translate(0.25, 0.25, 1.0).
            scale(0.5).
            translate(0.0, 0.0, -1.0)
    };

    let start = move |_| -> RuleEval<()> {

        RuleEval {
            geom: Rc::new(OpenMesh {
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
            }),
            final_geom: Rc::new(prim::empty_mesh()),
            children: vec![
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(recur.clone()), ctxt: () }),
                    xf: opening_xform(0.0),
                    vmap: vec![5,2,6,8],
                },
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(recur.clone()), ctxt: () }),
                    xf: opening_xform(1.0),
                    vmap: vec![4,1,5,8],
                },
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(recur.clone()), ctxt: () }),
                    xf: opening_xform(2.0),
                    vmap: vec![7,0,4,8],
                },
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(recur.clone()), ctxt: () }),
                    xf: opening_xform(3.0),
                    vmap: vec![6,3,7,8],
                },
                // TODO: These vertex mappings appear to be right.
                // Explain *why* they are right.
                // TODO: Factor out the repetition here.
            ],
        }
    };

    Rule { eval: Rc::new(start), ctxt: () }
}

#[derive(Copy, Clone)]
struct RamHornCtxt {
    depth: usize,
}

fn ramhorn_branch(depth: usize) -> Rule<RamHornCtxt> {

    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 0.8).
        rotate(&v, 0.4).
        scale(0.95);

    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let next = incr.transform(&seed);
    let geom = Rc::new(OpenMesh {
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
    });
    let final_geom = Rc::new(OpenMesh {
        verts: vec![],
        faces: vec![
            Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
            Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
        ],
    });

    let opening_xform = |i| {
        let r = std::f32::consts::FRAC_PI_2 * i;
        Transform::new().
            rotate(&nalgebra::Vector3::z_axis(), r).
            translate(0.25, 0.25, 0.0).
            scale(0.5)
    };
    
    let trans_verts = vec![
        // 'Top' vertices:
        vertex(-0.5, -0.5, 0.0),  //  0 (above 9)
        vertex(-0.5,  0.5, 0.0),  //  1 (above 10)
        vertex( 0.5,  0.5, 0.0),  //  2 (above 11)
        vertex( 0.5, -0.5, 0.0),  //  3 (above 12)
        // Top edge midpoints:
        vertex(-0.5,  0.0, 0.0),  //  4 (connects 0-1)
        vertex( 0.0,  0.5, 0.0),  //  5 (connects 1-2)
        vertex( 0.5,  0.0, 0.0),  //  6 (connects 2-3)
        vertex( 0.0, -0.5, 0.0),  //  7 (connects 3-0)
        // Top middle:
        vertex( 0.0,  0.0, 0.0),  //  8
    ];
    let trans_faces = vec![
        // two faces straddling edge from vertex 0:
        Tag::Parent(0), Tag::Body(0), Tag::Body(4),
        Tag::Parent(0), Tag::Body(7), Tag::Body(0),
        // two faces straddling edge from vertex 1:
        Tag::Parent(1), Tag::Body(1), Tag::Body(5),
        Tag::Parent(1), Tag::Body(4), Tag::Body(1),
        // two faces straddling edge from vertex 2:
        Tag::Parent(2), Tag::Body(2), Tag::Body(6),
        Tag::Parent(2), Tag::Body(5), Tag::Body(2),
        // two faces straddling edge from vertex 3:
        Tag::Parent(3), Tag::Body(3), Tag::Body(7),
        Tag::Parent(3), Tag::Body(6), Tag::Body(3),
        // four faces from edge (0,1), (1,2), (2,3), (3,0):
        Tag::Parent(0), Tag::Body(4), Tag::Parent(1),
        Tag::Parent(1), Tag::Body(5), Tag::Parent(2),
        Tag::Parent(2), Tag::Body(6), Tag::Parent(3),
        Tag::Parent(3), Tag::Body(7), Tag::Parent(0),
    ];
    let trans_geom = Rc::new(OpenMesh {
        verts: trans_verts.clone(),
        faces: trans_faces.clone(),
    });
    let trans_children = move |recur: RuleFn<RamHornCtxt>, ctxt: RamHornCtxt| {
        vec![
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(0.0),
                vmap: vec![5,2,6,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(1.0),
                vmap: vec![4,1,5,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(2.0),
                vmap: vec![7,0,4,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(3.0),
                vmap: vec![6,3,7,8],
            },
            // TODO: These vertex mappings appear to be right.
            // Explain *why* they are right.
            // TODO: Factor out the repetition here.
        ]
    };
    
    let tg = trans_geom.clone();
    // TODO: Why is that necessary?
    let recur = move |self_: Rc<Rule<RamHornCtxt>>| -> RuleEval<RamHornCtxt> {
        if self_.ctxt.depth <= 0 {
            RuleEval {
                geom: tg.clone(),
                final_geom: final_geom.clone(),
                // This final_geom will leave midpoint/centroid
                // vertices, but stopping here means none are
                // connected anyway - so they can just be ignored.
                children: trans_children(self_.eval.clone(), RamHornCtxt { depth }),
            }
        } else {
            let next_rule = Rule {
                eval: self_.eval.clone(),
                ctxt: RamHornCtxt { depth: self_.ctxt.depth - 1 },
            };
            RuleEval {
                geom: geom.clone(),
                final_geom: final_geom.clone(),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: incr,
                        vmap: vec![0,1,2,3],
                    },
                ],
            }
        }
    };
    
    let trans = move |self_: Rc<Rule<RamHornCtxt>>| -> RuleEval<RamHornCtxt> {
        RuleEval {
            geom: trans_geom.clone(),
            final_geom: Rc::new(prim::empty_mesh()),
            children: trans_children(Rc::new(recur.clone()), self_.ctxt),
        }
    };

    let start = move |self_: Rc<Rule<RamHornCtxt>>| -> RuleEval<RamHornCtxt> {
        RuleEval {
            geom: Rc::new(OpenMesh {
                verts: vec![
                    vertex(-0.5, -0.5, -0.5),
                    vertex(-0.5,  0.5, -0.5),
                    vertex( 0.5,  0.5, -0.5),
                    vertex( 0.5, -0.5, -0.5),
                ],
                faces: vec![
                    Tag::Body(0), Tag::Body(1), Tag::Body(2),
                    Tag::Body(0), Tag::Body(2), Tag::Body(3),
                ],
            }),
            final_geom: Rc::new(prim::empty_mesh()),
            children: vec![
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(trans.clone()), ctxt: self_.ctxt }),
                    xf: Transform::new(),
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    };
    
    Rule { eval: Rc::new(start), ctxt: RamHornCtxt { depth } }
}

/*
#[derive(Copy, Clone)]
struct CurveHorn {
    seed: [Vertex; 4],
    id_xform: Mat4,
    flip180: Mat4,
    incr: Mat4,
}

impl CurveHorn {

    fn test_thing(&self) {
        let f: Box<dyn Fn() -> RuleEval> = Rc::new(move || self.do_nothing());
        println!("{:p}", f);
    }

    fn do_nothing(&self) -> RuleEval {
        RuleEval {
            geom: prim::empty_mesh(),
            final_geom: prim::empty_mesh(),
            children: vec![
                Child {
                    rule: Rule { eval: Rc::new(move || self.do_nothing()) },
                    xf: self.id_xform,
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    }
    
    fn init() -> Rule {
        let y = &Vector3::y_axis();
        let c = CurveHorn {
            seed: [
                vertex(-0.5, -0.5, 0.0),
                vertex(-0.5,  0.5, 0.0),
                vertex( 0.5,  0.5, 0.0),
                vertex( 0.5, -0.5, 0.0),
            ],
            id_xform: nalgebra::geometry::Transform3::identity().to_homogeneous(),
            flip180: nalgebra::geometry::Rotation3::from_axis_angle(
                &nalgebra::Vector3::y_axis(),
                std::f32::consts::PI).to_homogeneous(),
            incr: geometry::Rotation3::from_axis_angle(y, 0.1).to_homogeneous() *
                Matrix4::new_scaling(0.95) *
                geometry::Translation3::new(0.0, 0.0, 0.2).to_homogeneous(),
        };
        Rule { eval: Rc::new(move || c.do_nothing()) }
    }
}
    fn start(&self) -> RuleEval {
        RuleEval {
            geom: OpenMesh {
                verts: self.seed.to_vec(),
                faces: vec![],
            },
            final_geom: prim::empty_mesh(),
            children: vec![
                Child {
                    rule: Rule { eval: Rc::new(move || self.recur()) },
                    xf: self.id_xform,
                    vmap: vec![0,1,2,3],
                },
                Child {
                    rule: Rule { eval: Rc::new(move || self.recur()) },
                    xf: self.flip180,
                    vmap: vec![3,2,1,0],
                },
            ],
        }
    }

    fn recur(&self) -> RuleEval {

        let verts = self.seed.clone();
        let next_verts: Vec<Vertex> = transform(&verts, &self.incr);
        
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
       g         Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
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
        
        RuleEval{
            geom: geom,
            final_geom: final_geom,
            children: vec![
                Child {
                    rule: Rule { eval: Rc::new(move || self.recur()) },
                    xf: self.incr,
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    }
}
*/

pub fn main() {
    
    fn run_test<S>(r: &Rc<Rule<S>>, iters: usize, name: &str, use_old: bool) {
        println!("---------------------------------------------------");
        println!("Running {} with {}...",
                 name, if use_old { "to_mesh" } else { "to_mesh_iter" });
        if false {
            let start = Instant::now();
            let n = 5;
            for _ in 0..n {
                Rule::to_mesh_iter(r.clone(), iters);
            }
            let elapsed = start.elapsed();
            println!("DEBUG: {} ms per run", elapsed.as_millis() / n);
        }
        let mesh_fn = if use_old { Rule::to_mesh } else { Rule::to_mesh_iter };
        let (mesh, nodes) = mesh_fn(r.clone(), iters);
        println!("Evaluated {} rules to {} verts", nodes, mesh.verts.len());
        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        mesh.write_stl_file(&fname).unwrap();
    }
    
    /*
    run_test(CubeThing::init(), Rule { eval: CubeThing::rec }, 3, "cube_thing");
    // this can't work on its own because the resultant OpenMesh still
    // has parent references:
    //run_test(Rule { eval: recur }, 100, "curve_horn_thing");
    run_test(CurveHorn::init(), Rule { eval: CurveHorn::start }, 100, "curve_horn2");
    run_test(RamHorn::init(), Rule { eval: RamHorn::start }, 200, "ram_horn");
    run_test(Twist::init(), Rule { eval: Twist::start }, 200, "twist");
    */

    //run_test_iter(CubeThing::init(), 3, "cube_thing2");
    //run_test_iter(CurveHorn::init(), 100, "curve_horn2_iter");
    //run_test_iter(RamHorn::init(), 100, "ram_horn2");
    // TODO: If I increase the above from 100 to ~150, Blender reports
    // that the very tips are non-manifold.  I am wondering if this is
    // some sort of numerical precision issue.
    
    //run_test_iter(Twist::init(1.0, 2), 100, "twist");
    // This is a stress test:
    // let f = 20;
    // run_test_iter(Twist::init(f as f32, 32), 100*f, "twist2");

    run_test(&Rc::new(cube_thing()), 3, "cube_thing3", false);
    run_test(&Rc::new(twist(1.0, 2)), 200, "twist", false);
    run_test(&Rc::new(ramhorn()), 100, "ram_horn3", false);
    run_test(&Rc::new(ramhorn_branch(6)), 31/*42*/, "ram_horn_branch", false);
}
