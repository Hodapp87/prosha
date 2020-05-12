use std::rc::Rc;

use nalgebra::*;
use rand::Rng;

use crate::util;
use crate::util::VecExt;
use crate::mesh::{Mesh, MeshFunc, VertexUnion};
use crate::xform::{Transform, Vertex, vertex, Mat4};
use crate::rule::{Rule, RuleFn, RuleEval, Child};
use crate::prim;

/*
pub fn cube_thing() -> Rule<()> {

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
*/

pub fn barbs() -> Rule<()> {

    let (b0, b1);
    let base_verts: Vec<VertexUnion> = vec_indexed![
        @b0: VertexUnion::Vertex(vertex(-0.5, -0.5, 0.0)),
             VertexUnion::Vertex(vertex(-0.5,  0.5, 0.0)),
             VertexUnion::Vertex(vertex( 0.5,  0.5, 0.0)),
        @b1: VertexUnion::Vertex(vertex( 0.5, -0.5, 0.0)),
    ];
    let n = base_verts.len();

    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 1.0).
        rotate(&Vector3::z_axis(), 0.15).
        rotate(&Vector3::x_axis(), 0.1).
        scale(0.95);

    let incr2: Transform = Transform::new().
        translate(0.0, 0.0, 1.5).
        rotate(&Vector3::y_axis(), -0.2).
        scale(0.8);

    let b = base_verts.clone();
    let barb = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
        let mut next_verts = b.clone();
        let (a0, a1) = next_verts.append_indexed(
            &mut (0..4).map(|i| VertexUnion::Arg(i)).collect()
        );

        let geom = util::parallel_zigzag(next_verts.clone(), b0..b1+1, a0..a1);
        /*let (vc, faces) = util::connect_convex(&next_verts, true);
        let final_geom = Rc::new(OpenMesh {
            verts: vec![vc],
            alias_verts: vec![],
            faces: faces,
        });
         */

        RuleEval {
            geom: Rc::new(geom.transform(&incr)),
            final_geom: Rc::new(prim::empty_meshfunc()), // TODO
            children: vec![
                Child {
                    rule: self_.clone(),
                    xf: incr2,
                    vmap: (0..n).collect(),
                }
            ]
        }
    };

    let barb_ = Rc::new(barb);

    let b = base_verts.clone();
    let main = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
        let mut next_verts = b.clone();
        let (a0, a1) = next_verts.append_indexed(
            &mut (0..4).map(|i| VertexUnion::Arg(i)).collect()
        );

        // TODO: Once I start doing the barbs this will go away
        let geom = util::parallel_zigzag(next_verts.clone(), b0..b1+1, a0..a1);
        /*
        let (vc, faces) = util::connect_convex(&next_verts, true);
        let final_geom = Rc::new(MeshFunc {
            verts: vec![vc],
            faces: faces,
        });
         */

        RuleEval {
            geom: Rc::new(geom.transform(&incr)),
            final_geom: Rc::new(prim::empty_meshfunc()), // TODO
            children: vec![
                Child {
                    rule: self_.clone(),
                    xf: incr,
                    vmap: (0..n).collect(),
                },
                Child {
                    rule:  Rc::new(Rule { eval: barb_.clone(), ctxt: () }),
                    xf: Transform::new().
                        translate(0.0, 0.0, 1.4).
                        rotate(&Vector3::y_axis(), -std::f32::consts::FRAC_PI_2).
                        scale(0.9),
                    vmap: vec![b0, b0 + 1, a0 + 1, a0],
                },
            ],
        }
    };

    let main_ = Rc::new(main);
    let base = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
        RuleEval {
            geom: Rc::new(MeshFunc {
                verts: base_verts.clone(),
                faces: vec![
                    b0, b0 + 1, b0 + 2,
                    b0, b0 + 2, b0 + 3,
                ],
            }),
            final_geom: Rc::new(prim::empty_meshfunc()),
            children: vec![
                Child {
                    rule: Rc::new(Rule { eval: main_.clone(), ctxt: () }),
                    xf: Transform::new(),
                    vmap: (0..n).collect(),
                },
            ],
        }
    };

    Rule { eval: Rc::new(base), ctxt: () }
}

/*
// Meant to be a copy of twist_from_gen from Python &
// automata_scratch, but has since acquired a sort of life of its own
pub fn twist(f: f32, subdiv: usize) -> Rule<()> {
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
    let dx0: f32 = 2.0;
    let dy: f32 = 0.1/f;
    let ang: f32 = 0.1/f;
    let count: usize = 4;
    
    // Quarter-turn in radians:
    let qtr = std::f32::consts::FRAC_PI_2;
    let y = Vector3::y_axis();

    let incr_inner = Transform::new().translate(-dx0, 0.0, 0.0).rotate(&y, ang).translate(dx0, dy, 0.0);
    let incr_outer = Transform::new().translate(-dx0*2.0, 0.0, 0.0).rotate(&y, ang/2.0).translate(dx0*2.0, dy, 0.0);

    let seed2 = seed.clone();
    // TODO: Why do I need the above?
    // TODO: Could a macro get rid of some of this or would it just be
    // equally cumbersome because I'd have to sort of pass 'seed'
    // explicitly?
    let recur = move |incr: Transform| -> RuleFn<()> {

        let seed_next = incr.transform(&seed2);

        let geom = Rc::new(util::zigzag_to_parent(seed_next.clone(), n));
        // TODO: Cleanliness fix - why not just make these return meshes?
        let (vc, faces) = util::connect_convex(&seed_next, true);
        let final_geom = Rc::new(OpenMesh {
            verts: vec![vc],
            alias_verts: vec![],
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
            (OpenMesh { verts: vs, faces: f, alias_verts: vec![] }, c)
        };
        
        // Generate 'count' children, shifted/rotated differently:
        let inner = (0..count).map(|i| child(incr_inner, dx0, i, 0.0, 1.0));
        //let outer = (0..count).map(|i| child(incr_outer, dx0*2.0, i, qtr/2.0, 2.0));
        let outer = (0..0).map(|i| child(incr_outer, dx0*2.0, i, qtr/2.0, 2.0));

        RuleEval::from_pairs(inner.chain(outer), prim::empty_mesh())
    };
    
    Rule { eval: Rc::new(start), ctxt: () }
}

#[derive(Copy, Clone)]
pub struct NestSpiral2Ctxt {
    init: bool,
    stack: [Transform; 2],
}

pub fn nest_spiral_2() -> Rule<NestSpiral2Ctxt> {
    let subdiv = 8;
    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let seed = util::subdivide_cycle(&seed, subdiv);

    let n = seed.len();
    let geom = Rc::new(util::zigzag_to_parent(seed.clone(), n));
    let (vc, faces) = util::connect_convex(&seed, true);
    let final_geom = Rc::new(OpenMesh {
        verts: vec![vc],
        alias_verts: vec![],
        faces: faces,
    });

    let rad = 1.0;
    let dz = 0.1;
    let rz = 0.1;
    let rad2 = 4.0;
    let rz2 = 0.1;

    let recur = move |self_: Rc<Rule<NestSpiral2Ctxt>>| -> RuleEval<NestSpiral2Ctxt> {
        //let x = &Vector3::x_axis();
        let z = &Vector3::z_axis();
        let stack = self_.ctxt.stack;
        let next_rule = Rule {
            eval: self_.eval.clone(),
            ctxt: NestSpiral2Ctxt {
                init: false,
                stack: [
                    Transform::new().rotate(z, rz2) * stack[0],
                    Transform::new().translate(0.0, 0.0, dz).rotate(z, rz) * stack[1],
                ],
            },
        };
        let xf = stack.iter().fold(Transform::new(), |acc,m| acc * (*m));
        if self_.ctxt.init {
            let mut s2 = seed.clone();
            let (centroid, f) = util::connect_convex(&s2, false);
            s2.push(centroid);
            let n2 = s2.len();
            let g = OpenMesh { verts: s2, faces: f, alias_verts: vec![] };
            RuleEval {
                geom: Rc::new(g.transform(&xf)),
                final_geom: Rc::new(prim::empty_mesh()),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n2).collect(),
                    },
                ],
            }
        } else {
            RuleEval {
                geom: Rc::new(geom.transform(&xf)),
                final_geom: Rc::new(final_geom.transform(&xf)),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n).collect(),
                    },
                ],
            }
        }
    };

    let count = 3;

    let r = Rc::new(recur);
    let start = move |self_: Rc<Rule<NestSpiral2Ctxt>>| -> RuleEval<NestSpiral2Ctxt> {
        let z = &Vector3::z_axis();
        let child = |i: usize| -> Child<NestSpiral2Ctxt> {
            let ang = std::f32::consts::PI * 2.0 * (i as f32) / (count as f32);
            Child {
                rule: Rc::new(Rule {
                    eval: r.clone(),
                    ctxt: NestSpiral2Ctxt {
                        init: true,
                        stack: [
                            Transform::new().translate(rad2, 0.0, 0.0),
                            Transform::new().rotate(z, ang).translate(rad, 0.0, 0.0),
                        ],
                    },
                }),
                xf: Transform::new(),
                vmap: vec![], // no parent vertices
            }
        };

        RuleEval {
            geom: Rc::new(prim::empty_mesh()),
            final_geom: Rc::new(prim::empty_mesh()),
            children: (0..count).map(child).collect(),
        }
    };

    Rule {
        eval: Rc::new(start),
        ctxt: NestSpiral2Ctxt {
            init: true,
            stack: [ // doesn't matter
                Transform::new(),
                Transform::new(),
            ],
        },
    }
}


#[derive(Copy, Clone)]
pub struct TorusCtxt {
    init: bool,
    count: usize,
    stack: [Transform; 3],
}

pub fn twisty_torus() -> Rule<TorusCtxt> {
    let subdiv = 8;
    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let xf = Transform::new().rotate(&Vector3::x_axis(), -0.9);
    let seed = util::subdivide_cycle(&xf.transform(&seed), subdiv);
    
    let n = seed.len();
    let geom = Rc::new(util::zigzag_to_parent(seed.clone(), n));
    let (vc, faces) = util::connect_convex(&seed, true);
    let final_geom = Rc::new(OpenMesh {
        verts: vec![vc],
        alias_verts: (0..(n+1)).collect(),
        faces: faces,
    });

    let rad = 1.0;
    let rad2 = 8.0;
    let rad3 = 24.0;
    let rz3 = 0.0004;
    let dx = 0.00;
    let rx = 0.01;
    let rz = 0.30;
    let ang = 0.1;

    let recur = move |self_: Rc<Rule<TorusCtxt>>| -> RuleEval<TorusCtxt> {
        let x = &Vector3::x_axis();
        let z = &Vector3::z_axis();
        let stack = self_.ctxt.stack;
        let count = self_.ctxt.count;
        let next_rule = Rule {
            eval: self_.eval.clone(),
            ctxt: TorusCtxt {
                init: false,
                count: count + 1,
                stack: [
                    Transform::new().rotate(z, rz3) * stack[0],
                    Transform::new().translate(dx, 0.0, 0.0).rotate(x, rx) * stack[1],
                    Transform::new().rotate(z, rz) * stack[2],
                ],
            },
        };
        let xf = stack.iter().fold(Transform::new(), |acc,m| acc * (*m));
        if self_.ctxt.init {
            let mut s2 = seed.clone();
            let (centroid, f) = util::connect_convex(&s2, false);
            s2.push(centroid);
            let n2 = s2.len();
            let g = OpenMesh { verts: s2, faces: f, alias_verts: vec![] };
            RuleEval {
                geom: Rc::new(g.transform(&xf)),
                final_geom: Rc::new(prim::empty_mesh()),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n2).collect(),
                    },
                ],
            }
        } else {
            RuleEval {
                geom: Rc::new(geom.transform(&xf)),
                final_geom: Rc::new(final_geom.transform(&xf)),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n).collect(),
                    },
                ],
            }
        }
    };

    Rule {
        eval: Rc::new(recur),
        ctxt: TorusCtxt {
            init: true,
            count: 0,
            stack: [
                Transform::new().translate(0.0, rad3, 0.0),
                Transform::new().translate(0.0, rad2, 0.0),
                Transform::new().translate(rad, 0.0, 0.0),
            ],
        },
    }
}

pub fn twisty_torus_hardcode() -> Rule<()> {
    let subdiv = 8;
    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let xf = Transform::new().rotate(&Vector3::x_axis(), -0.9);
    let seed = util::subdivide_cycle(&xf.transform(&seed), subdiv);
    let incr = Transform { mtx: Mat4::from_vec(vec![
        0.955234,    0.29576725, -0.0070466697, 0.0,
       -0.29581502,  0.9552189,  -0.007100463,  0.0,
        0.004630968, 0.008867174, 0.99994993,   0.0,
       -0.034161568, 0.290308,    0.07295418,   0.9999999,
    ])};
    
    let n = seed.len();
    
    let next = incr.transform(&seed);
    let geom = Rc::new(util::zigzag_to_parent(next.clone(), n));
    let (vc, faces) = util::connect_convex(&next, true);
    let final_geom = Rc::new(OpenMesh {
        verts: vec![vc],
        alias_verts: (0..(n+1)).collect(), // TODO: Fix parent/connect_convex
        faces: faces,
    });

    let rad = 1.0;
    let rad2 = 8.0;
    let rad3 = 24.0;

    let start = Transform::new().translate(0.0, rad3, 0.0) * Transform::new().translate(0.0, rad2, 0.0) * Transform::new().translate(rad, 0.0, 0.0);
    
    let recur = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
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

    let start = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
        let mut s2 = seed.clone();
        let (centroid, f) = util::connect_convex(&s2, false);
        s2.push(centroid);
        let n2 = s2.len();
        let g = OpenMesh { verts: s2, faces: f, alias_verts: vec![] };
        RuleEval {
            geom: Rc::new(g.transform(&xf)),
            final_geom: Rc::new(prim::empty_mesh()),
            children: vec![
                Child {
                    rule: Rc::new(Rule { eval: Rc::new(recur.clone()), ctxt: () }),
                    xf: incr,
                    vmap: (0..n2).collect(),
                },
            ],
        }
    };

    Rule {
        eval: Rc::new(start),
        ctxt: (),
    }
}

// This was a mistake that I'd like to understand later:
#[derive(Copy, Clone)]
pub struct WindChimeCtxt {
    init: bool,
    count: usize,
    stack: [Transform; 3],
}

pub fn wind_chime_mistake_thing() -> Rule<WindChimeCtxt> {
    let subdiv = 8;
    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let seed = util::subdivide_cycle(&seed, subdiv);
    
    let n = seed.len();
    let geom = Rc::new(util::zigzag_to_parent(seed.clone(), n));
    let (vc, faces) = util::connect_convex(&seed, true);
    let final_geom = Rc::new(OpenMesh {
        verts: vec![vc],
        alias_verts: (0..(n + 1)).collect(), // TODO: Check with parents (zigzag/connect_convex)
        faces: faces,
    });

    let rad = 1.0;
    let rad2 = 8.0;
    let dx0 = 2.0;
    let ang = 0.1;

    let recur = move |self_: Rc<Rule<WindChimeCtxt>>| -> RuleEval<WindChimeCtxt> {
        let x = &Vector3::x_axis();
        let z = &Vector3::z_axis();
        let stack = self_.ctxt.stack;
        let count = self_.ctxt.count;
        let next_rule = Rule {
            eval: self_.eval.clone(),
            ctxt: WindChimeCtxt {
                init: false,
                count: count + 1,
                stack: [
                    Transform::new().rotate(x, 0.01) * stack[0],
                    // stack[0], //Transform::new().rotate(z, 0.05 * (count as f32)).translate(0.0, rad2, 0.0),
                    Transform::new().rotate(z, 0.30) * stack[1],
                    Transform::new().translate(0.1, 0.0, 0.0) * stack[2],
                ],
            },
        };
        let xf = stack.iter().fold(Transform::new(), |acc,m| acc * (*m));
        if self_.ctxt.init {
            let mut s2 = seed.clone();
            let (centroid, f) = util::connect_convex(&s2, false);
            s2.push(centroid);
            let n2 = s2.len();
            let g = OpenMesh { verts: s2, faces: f, alias_verts: vec![] };
            RuleEval {
                geom: Rc::new(g.transform(&xf)),
                final_geom: Rc::new(prim::empty_mesh()),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n2).collect(),
                    },
                ],
            }
        } else {
            RuleEval {
                geom: Rc::new(geom.transform(&xf)),
                final_geom: Rc::new(final_geom.transform(&xf)),
                children: vec![
                    Child {
                        rule: Rc::new(next_rule),
                        xf: Transform::new(),
                        vmap: (0..n).collect(),
                    },
                ],
            }
        }
    };

    Rule {
        eval: Rc::new(recur),
        ctxt: WindChimeCtxt {
            init: true,
            count: 0,
            stack: [
                Transform::new().translate(0.0, rad2, 0.0),
                Transform::new().translate(rad, 0.0, 0.0),
                Transform::new(), // .translate(dx0, 0.0, 0.0),
            ],
        },
    }
}

pub fn ramhorn() -> Rule<()> {

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
        alias_verts: vec![0, 1, 2, 3],
        verts: next,
        faces: vec![
            5, 0, 4,
            1, 0, 5,
            6, 1, 5,
            2, 1, 6,
            7, 2, 6,
            3, 2, 7,
            4, 3, 7,
            0, 3, 4,
        ],
    });
    let final_geom = Rc::new(OpenMesh {
        verts: vec![],
        alias_verts: vec![0, 1, 2, 3],
        faces: vec![
            0, 2, 1,
            0, 3, 2,
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
                alias_verts: vec![],
                faces: vec![
                    // bottom face:
                    9, 10, 11,
                    9, 11, 12,
                    // two faces straddling edge from vertex 0:
                    9, 0, 4,
                    9, 7, 0,
                    // two faces straddling edge from vertex 1:
                    10, 1, 5,
                    10, 4, 1,
                    // two faces straddling edge from vertex 2:
                    11, 2, 6,
                    11, 5, 2,
                    // two faces straddling edge from vertex 3:
                    12, 3, 7,
                    12, 6, 3,
                    // four faces from edge (0,1, (1,2, (2,3, (3,0):
                    9, 4, 10,
                    10, 5, 11,
                    11, 6, 12,
                    12, 7, 9,
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
pub struct RamHornCtxt {
    depth: usize,
}

pub fn ramhorn_branch(depth: usize, f: f32) -> Rule<RamHornCtxt> {

    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 0.8 * f).
        rotate(&v, 0.4 * f).
        scale(1.0 - (1.0 - 0.95)*f);

    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let next = incr.transform(&seed);
    let geom = Rc::new(OpenMesh {
        verts: next,
        alias_verts: vec![],
        faces: util::parallel_zigzag_faces(4),
        // TODO: Fix this (parallel_zigzag_faces has parents)
    });
    let final_geom = Rc::new(OpenMesh {
        verts: vec![],
        alias_verts: vec![0, 1, 2, 3],
        faces: vec![
            0, 2, 1,
            0, 3, 2,
        ],
    });

    let opening_xform = |i| {
        let r = std::f32::consts::FRAC_PI_2 * i;
        Transform::new().
            rotate(&nalgebra::Vector3::z_axis(), r).
            translate(0.25, 0.25, 0.0).
            scale(0.5)
    };

    // 'transition' geometry (when something splits):
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
        0, 4, 8,
        0, 11, 4,
        // two faces straddling edge from vertex 1:
        1, 5, 9,
        1, 8, 5,
        // two faces straddling edge from vertex 2:
        2, 6, 10,
        2, 9, 6,
        // two faces straddling edge from vertex 3:
        3, 7, 11,
        3, 10, 7,
        // four faces from edge (0,1), (1,2), (2,3), (3,0):
        0, 8, 1,
        1, 9, 2,
        2, 10, 3,
        3, 11, 0,
    ];
    let trans_geom = Rc::new(OpenMesh {
        alias_verts: vec![0, 1, 2, 3],
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
                verts: Transform::new().translate(0.0, 0.0, -0.5).transform(&seed),
                alias_verts: vec![],
                faces: vec![
                    0, 1, 2,
                    0, 2, 3,
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

#[derive(Copy, Clone)]
pub struct RamHornCtxt2 {
    depth: usize,
}

pub fn ramhorn_branch_random(depth: usize, f: f32) -> Rule<RamHornCtxt2> {

    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 0.8 * f).
        rotate(&v, 0.4 * f).
        scale(1.0 - (1.0 - 0.95)*f);

    let seed = vec![
        vertex(-0.5, -0.5, 0.0),
        vertex(-0.5,  0.5, 0.0),
        vertex( 0.5,  0.5, 0.0),
        vertex( 0.5, -0.5, 0.0),
    ];
    let next = incr.transform(&seed);
    let geom = Rc::new(OpenMesh {
        verts: next,
        faces: util::parallel_zigzag_faces(4),
        alias_verts: vec![],
        // TODO: Fix parents with parallel_zigzag
    });
    let final_geom = Rc::new(OpenMesh {
        verts: vec![],
        alias_verts: vec![0, 1, 2, 3],
        faces: vec![
            0, 2, 1,
            0, 3, 2,
        ],
    });

    let opening_xform = |i| {
        let r = std::f32::consts::FRAC_PI_2 * i;
        Transform::new().
            rotate(&nalgebra::Vector3::z_axis(), r).
            translate(0.25, 0.25, 0.0).
            scale(0.5)
    };

    // 'transition' geometry (when something splits):
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
        0, 4, 8,
        0, 11, 4,
        // two faces straddling edge from vertex 1:
        1, 5, 9,
        1, 8, 5,
        // two faces straddling edge from vertex 2:
        2, 6, 10,
        2, 9, 6,
        // two faces straddling edge from vertex 3:
        3, 7, 11,
        3, 10, 7,
        // four faces from edge (0,1), (1,2), (2,3), (3,0):
        0, 8, 1,
        1, 9, 2,
        2, 10, 3,
        3, 11, 0,
    ];
    let trans_geom = Rc::new(OpenMesh {
        alias_verts: vec![0, 1, 2, 3],
        verts: trans_verts.clone(),
        faces: trans_faces.clone(),
    });
    let trans_children = move |recur: RuleFn<RamHornCtxt2>, ctxt: RamHornCtxt2| {
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
    let recur = move |self_: Rc<Rule<RamHornCtxt2>>| -> RuleEval<RamHornCtxt2> {
        if self_.ctxt.depth <= 0 {
            let d2 = rand::thread_rng().gen_range(2, 60);
            RuleEval {
                geom: tg.clone(),
                final_geom: final_geom.clone(),
                // This final_geom will leave midpoint/centroid
                // vertices, but stopping here means none are
                // connected anyway - so they can just be ignored.
                children: trans_children(self_.eval.clone(), RamHornCtxt2 { depth: d2 }),
            }
        } else {
            let next_rule = Rule {
                eval: self_.eval.clone(),
                ctxt: RamHornCtxt2 { depth: self_.ctxt.depth - 1 },
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

    let trans = move |self_: Rc<Rule<RamHornCtxt2>>| -> RuleEval<RamHornCtxt2> {
        RuleEval {
            geom: trans_geom.clone(),
            final_geom: Rc::new(prim::empty_mesh()),
            children: trans_children(Rc::new(recur.clone()), self_.ctxt),
        }
    };

    let start = move |self_: Rc<Rule<RamHornCtxt2>>| -> RuleEval<RamHornCtxt2> {
        RuleEval {
            geom: Rc::new(OpenMesh {
                verts: Transform::new().translate(0.0, 0.0, -0.5).transform(&seed),
                alias_verts: vec![],
                faces: vec![
                    0, 1, 2,
                    0, 2, 3,
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

    Rule { eval: Rc::new(start), ctxt: RamHornCtxt2 { depth } }
}
 */

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
