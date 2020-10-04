use std::rc::Rc;
use std::f32::consts::{FRAC_PI_2, FRAC_PI_3};
use std::f32;

use nalgebra::*;
use rand::Rng;

use crate::util;
use crate::util::VecExt;
use crate::mesh::{Mesh, MeshFunc, VertexUnion, vert_args};
use crate::xform::{Transform, Vertex, vertex, id};
use crate::rule::{Rule, RuleEval, Child};
use crate::prim;
use crate::dcel;
use crate::dcel::{VertSpec};

pub fn cube_thing() -> Rule<()> {

    // Quarter-turn in radians:
    let qtr = FRAC_PI_2;
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
            geom: Rc::new(prim::cube().to_meshfunc()),
            final_geom: Rc::new(prim::empty_mesh().to_meshfunc()),
            children: xforms.map(move |xf| Child {
                rule: self_.clone(),
                xf: xf,
                arg_vals: vec![],
            }).collect(),
        }
    };

    Rule { eval: Rc::new(rec), ctxt: () }
}

pub fn barbs(random: bool) -> Rule<()> {

    let (b0, bn);
    let base_verts: Vec<VertexUnion> = vec_indexed![
        @b0 VertexUnion::Vertex(vertex(-0.5, -0.5, 0.0)),
        VertexUnion::Vertex(vertex(-0.5,  0.5, 0.0)),
        VertexUnion::Vertex(vertex( 0.5,  0.5, 0.0)),
        VertexUnion::Vertex(vertex( 0.5, -0.5, 0.0)),
        @bn,
    ];

    let barb_incr = |random| {
        if random {
            let t = rand::thread_rng().gen_range(0.45, 0.55);
            let s = rand::thread_rng().gen_range(0.7, 0.9);
            let ry = rand::thread_rng().gen_range(-0.3, -0.1);
            let rx = rand::thread_rng().gen_range(-0.04, 0.04);
            let rz = rand::thread_rng().gen_range(-0.04, 0.04);
            id().translate(0.0, 0.0, t).
                rotate(&Vector3::y_axis(), ry).
                rotate(&Vector3::x_axis(), rx).
                rotate(&Vector3::z_axis(), rz).
                scale(s)
        } else {
            id().translate(0.0, 0.0, 0.5).
                rotate(&Vector3::y_axis(), -0.2).
                scale(0.8)
        }
    };

    let barb = rule_fn!(() => |self_, base_verts| {
        let mut next_verts = base_verts;
        let (a0, a1) = next_verts.append_indexed(vert_args(0..4));

        let geom = util::parallel_zigzag(next_verts, b0..bn, a0..a1);
        let final_geom = MeshFunc {
            verts: vert_args(0..4),
            faces: vec![ 0, 2, 1,   0, 3, 2 ],
        };

        let b = barb_incr(random);

        RuleEval {
            geom: Rc::new(geom.transform(&b)),
            final_geom: Rc::new(final_geom), // no transform needed (no vertices)
            children: vec![ child_iter!(self_, b, b0..bn) ],
        }
    });

    let main_barb_xf = |i| {
        id().rotate(&Vector3::z_axis(), -FRAC_PI_2 * (i as f32)).
        rotate(&Vector3::y_axis(), -FRAC_PI_2).
        translate(0.5, 0.0, 0.5)
    };
    let main_incr = |random| {
        if random {
            //let t = rand::thread_rng().gen_range(0.75, 1.25);
            let s = rand::thread_rng().gen_range(0.85, 1.10);
            let rz = rand::thread_rng().gen_range(0.05, 0.25);
            let rx = rand::thread_rng().gen_range(0.08, 0.12);
            id().translate(0.0, 0.0, 1.0).
                rotate(&Vector3::z_axis(), rz).
                rotate(&Vector3::x_axis(), rx).
                scale(s)
        } else {
            id().translate(0.0, 0.0, 1.0).
                rotate(&Vector3::z_axis(), 0.15).
                rotate(&Vector3::x_axis(), 0.1).
                scale(0.95)
        }
    };

    let main = rule_fn!(() => |self_, base_verts| {
        let mut next_verts = base_verts;
        let (a0, _) = next_verts.append_indexed(vert_args(0..4));

        // This contributes no faces of its own - just vertices.
        let geom = MeshFunc { verts: next_verts.clone(), faces: vec![] };
        // (unless recursion ends here, of course)
        let final_geom = MeshFunc {
            verts: vert_args(0..4),
            faces: vec![ 0, 2, 1,   0, 3, 2 ],
        };

        RuleEval {
            geom: Rc::new(geom),
            final_geom: Rc::new(final_geom),
            children: vec![
                child_iter!(self_, main_incr(random), b0..bn),
                child!(rule!(barb, ()), main_barb_xf(0), b0 + 0, b0 + 1, a0 + 1, a0 + 0),
                child!(rule!(barb, ()), main_barb_xf(1), b0 + 1, b0 + 2, a0 + 2, a0 + 1),
                child!(rule!(barb, ()), main_barb_xf(2), b0 + 2, b0 + 3, a0 + 3, a0 + 2),
                child!(rule!(barb, ()), main_barb_xf(3), b0 + 3, b0 + 0, a0 + 0, a0 + 3),
                // TODO: Factor out repetition?
            ],
        }
    });

    let base = rule_fn!(() => |_s, base_verts| {
        RuleEval {
            geom: Rc::new(MeshFunc {
                verts: base_verts,
                faces: vec![ b0, b0 + 1, b0 + 2,   b0, b0 + 2, b0 + 3 ],
            }),
            // TODO: This might be buggy and leave some vertices lying around
            final_geom: Rc::new(prim::empty_meshfunc()),
            children: vec![ child_iter!(rule!(main, ()), id(), b0..bn) ],
        }
    });

    //rule!(Rc::new(base), ())
    Rule { eval: base, ctxt: () }
}

pub fn sierpinski() -> Rule<()> {

    // Initial height step:
    let dz = 0.10;
    // 'Extra' z rotation (0.0 for normal Sierpinski)
    let dr = 0.1;
    // Scale factor (0.5 for normal Sierpinski)
    let s = 0.51;

    let rt3 = (3.0).sqrt();

    // Indices:
    // b+0,b+1,b+2 = base vertices
    // t+0,t+1,t+2 = 'top' vertices above base
    // tm01, tm12, tm20 = midpoints of (t0,t1), (t1,t2), (t2,t0).
    let (b, t, tm01, tm12, tm20, n);
    let base_verts: Vec<VertexUnion> = {
        let v0 = vertex(rt3/3.0, 0.0, 0.0);
        let v1 = vertex(-rt3/6.0, 1.0/2.0, 0.0);
        let v2 = vertex(-rt3/6.0, -1.0/2.0, 0.0);
        let v0b = v0 + vertex(0.0, 0.0, dz);
        let v1b = v1 + vertex(0.0, 0.0, dz);
        let v2b = v2 + vertex(0.0, 0.0, dz);
        vec_indexed![
            @b  VertexUnion::Vertex(v0),
            VertexUnion::Vertex(v1),
            VertexUnion::Vertex(v2),
            @t  VertexUnion::Vertex(v0b),
            VertexUnion::Vertex(v1b),
            VertexUnion::Vertex(v2b),
            @tm01 VertexUnion::Vertex((v0b+v1b)/2.0),
            @tm12 VertexUnion::Vertex((v1b+v2b)/2.0),
            @tm20 VertexUnion::Vertex((v2b+v0b)/2.0),
            @n,
        ]
    };

    let tri_split = move |i| {
        let rt3 = (3.0).sqrt();
        let angle = 2.0 * FRAC_PI_3 * (i as f32) + dr;
        id().
            rotate(&Vector3::z_axis(), angle).
            translate(rt3/12.0, 0.0, 0.0).
            scale(s).
            translate(0.0, 0.0, dz)
    };

    let split = rule_fn!(() => |_s, base_verts| {

        let mut next_verts = base_verts.clone();
        let (a0, _) = next_verts.append_indexed(vert_args(0..3));

        RuleEval {
            geom: Rc::new(MeshFunc {
                verts: next_verts,
                faces: vec![
                    //a0, a0+1, a0+2,
                    // Outer:
                    tm01, a0+1,  t+1,
                    tm01, t+0,  a0+0,
                    tm01, a0+0,  a0+1,
                    tm12, a0+2,  t+2,
                    tm12, t+1,  a0+1,
                    tm12, a0+1,  a0+2,
                    tm20, a0+0,  t+0,
                    tm20, t+2,  a0+2,
                    tm20, a0+2,  a0+0,
                    // Inner:
                    tm01, tm12, tm20,
                ],
            }),
            final_geom: Rc::new(MeshFunc {
                verts: vert_args(0..n), // just duplicate same verts
                faces: vec![
                    t+0, tm01, tm20,
                    t+1, tm12, tm01,
                    t+2, tm20, tm12,
                ],
            }),
            children: vec![
                child!(_s, tri_split(0), t+0, tm01, tm20),
                child!(_s, tri_split(1), t+1, tm12, tm01),
                child!(_s, tri_split(2), t+2, tm20, tm12),
            ],
        }
    });

    let base = rule_fn!(() => |_s, base_verts| {
        RuleEval {
            geom: Rc::new(MeshFunc {
                verts: base_verts,
                faces: vec![
                    // Outer:
                    tm01, b+1,  t+1,
                    tm01, t+0,  b+0,
                    tm01, b+0,  b+1,
                    tm12, b+2,  t+2,
                    tm12, t+1,  b+1,
                    tm12, b+1,  b+2,
                    tm20, b+0,  t+0,
                    tm20, t+2,  b+2,
                    tm20, b+2,  b+0,
                    // Inner:
                    tm01, tm12, tm20,
                    // Bottom:
                    b+2, b+1, b+0,
                ],
            }),
            final_geom: Rc::new(MeshFunc {
                verts: vec![],
                faces: vec![],
            }),
            children: vec![
                child!(rule!(split, ()), tri_split(0), t+0, tm01, tm20),
                child!(rule!(split, ()), tri_split(1), t+1, tm12, tm01),
                child!(rule!(split, ()), tri_split(2), t+2, tm20, tm12),
            ],
        }
    });

    Rule { eval: base, ctxt: () }
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
    let qtr = FRAC_PI_2;
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

        //let vc = util::centroid(&seed_next);
        //let faces = util::connect_convex(0..n, n, true);
        let geom = util::parallel_zigzag(seed_next, 0..n, 0..n);
        let final_geom = MeshFunc {
            verts: vec![],
            faces: vec![],
            // TODO: get actual verts here
        };
        
        let c = move |self_: Rc<Rule<()>>| -> RuleEval<()> {
            RuleEval {
                geom: geom.clone(),
                final_geom: final_geom.clone(),
                children: vec![
                    Child {
                        rule: self_.clone(),
                        xf: incr,
                        arg_vals: (0..n).collect(),
                    },
                ],
            }
        };
        Rc::new(c)
    };
    // TODO: Can a macro do anything to clean up some of the
    // repetition with HOFs & closures?
    
    let start = move |_| -> RuleEval<()> {
        
        let child = |incr, dx, i, ang0, div| -> (MeshFunc, Child<()>) {
            let xform = Transform::new().
                rotate(&y, ang0 + (qtr / div * (i as f32))).
                translate(dx, 0.0, 0.0);
            
            let c = Child {
                rule: Rc::new(Rule { eval: (recur.clone())(incr), ctxt: () }),
                // TODO: Cleanliness fix - can macros clean up above?
                xf: xform,
                arg_vals: (0..(n+1)).collect(),
                // N.B. n+1, not n. the +1 is for the centroid below.
            };
            let mut vs = xform.transform(&seed);
            // and in the process, generate faces for these seeds:
            let (centroid, f) = util::connect_convex(&vs, false);
            vs.push(centroid);
            (MeshFunc { verts: vs, faces: f }, c)
        };
        
        // Generate 'count' children, shifted/rotated differently:
        let inner = (0..count).map(|i| child(incr_inner, dx0, i, 0.0, 1.0));
        //let outer = (0..count).map(|i| child(incr_outer, dx0*2.0, i, qtr/2.0, 2.0));
        let outer = (0..0).map(|i| child(incr_outer, dx0*2.0, i, qtr/2.0, 2.0));

        RuleEval::from_pairs(inner.chain(outer), prim::empty_mesh())
    };
    
    Rule { eval: Rc::new(start), ctxt: () }
}
*/

/*
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
                        arg_vals: (0..n2).collect(),
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
                        arg_vals: (0..n).collect(),
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
            let ang = PI * 2.0 * (i as f32) / (count as f32);
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
                arg_vals: vec![], // no parent vertices
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
    let geom = util::parallel_zigzag(seed, 0..n, n..(2*n));
    // TODO: where are parent Args?
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
                        arg_vals: (0..n2).collect(),
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
                        arg_vals: (0..n).collect(),
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
 */

/*
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
                    arg_vals: (0..n).collect(),
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
                    arg_vals: (0..n2).collect(),
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
                        arg_vals: (0..n2).collect(),
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
                        arg_vals: (0..n).collect(),
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
*/

pub fn ramhorn() -> Rule<()> {

    let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
    let incr: Transform = Transform::new().
        translate(0.0, 0.0, 0.8).
        rotate(&v, 0.3).
        scale(0.9);

    let (a0, a1, a2, a3, s4, s5, s6, s7);
    let seed = vec_indexed![
        @a0 VertexUnion::Arg(0),
        @a1 VertexUnion::Arg(1),
        @a2 VertexUnion::Arg(2),
        @a3 VertexUnion::Arg(3),
        @s4 VertexUnion::Vertex(vertex(-0.5, -0.5, 1.0)),
        @s5 VertexUnion::Vertex(vertex(-0.5,  0.5, 1.0)),
        @s6 VertexUnion::Vertex(vertex( 0.5,  0.5, 1.0)),
        @s7 VertexUnion::Vertex(vertex( 0.5, -0.5, 1.0)),
    ];
    let geom = MeshFunc {
        verts: seed,
        faces: vec![
            s5, a0, s4,
            a1, a0, s5,
            s6, a1, s5,
            a2, a1, s6,
            s7, a2, s6,
            a3, a2, s7,
            s4, a3, s7,
            a0, a3, s4,
        ],
    };
    let final_geom = MeshFunc {
        verts: vec![
            VertexUnion::Arg(s4),
            VertexUnion::Arg(s5),
            VertexUnion::Arg(s6),
            VertexUnion::Arg(s7),
        ],
        // TODO: Factor out this repetition
        faces: vec![
            0, 2, 1,
            0, 3, 2,
        ],
    };

    let geom2 = Rc::new(geom.transform(&incr));
    let fgeom2 = Rc::new(final_geom);
    let recur = rule_fn!(() => |self_, geom2, fgeom2| {
        RuleEval {
            geom: geom2,
            final_geom: fgeom2,
            children: vec![
                child!(self_, incr, s4, s5, s6, s7),
            ],
        }
    });
    
    let opening_xform = |i| {
        let r = FRAC_PI_2 * (i as f32);
        Transform::new().
            rotate(&nalgebra::Vector3::z_axis(), r).
            translate(0.25, 0.25, 1.0).
            scale(0.5).
            translate(0.0, 0.0, -1.0)
    };

    let start = move |_| -> RuleEval<()> {

        RuleEval {
            geom: Rc::new(MeshFunc {
                verts: vec![
                    // 'Top' vertices:
                    VertexUnion::Vertex(vertex(-0.5, -0.5, 1.0)),  //  0 (above 9)
                    VertexUnion::Vertex(vertex(-0.5,  0.5, 1.0)),  //  1 (above 10)
                    VertexUnion::Vertex(vertex( 0.5,  0.5, 1.0)),  //  2 (above 11)
                    VertexUnion::Vertex(vertex( 0.5, -0.5, 1.0)),  //  3 (above 12)
                    // Top edge midpoints:
                    VertexUnion::Vertex(vertex(-0.5,  0.0, 1.0)),  //  4 (connects 0-1)
                    VertexUnion::Vertex(vertex( 0.0,  0.5, 1.0)),  //  5 (connects 1-2)
                    VertexUnion::Vertex(vertex( 0.5,  0.0, 1.0)),  //  6 (connects 2-3)
                    VertexUnion::Vertex(vertex( 0.0, -0.5, 1.0)),  //  7 (connects 3-0)
                    // Top middle:
                    VertexUnion::Vertex(vertex( 0.0,  0.0, 1.0)),  //  8
                    // 'Bottom' vertices:
                    VertexUnion::Vertex(vertex(-0.5, -0.5, 0.0)),  //  9
                    VertexUnion::Vertex(vertex(-0.5,  0.5, 0.0)),  // 10
                    VertexUnion::Vertex(vertex( 0.5,  0.5, 0.0)),  // 11
                    VertexUnion::Vertex(vertex( 0.5, -0.5, 0.0)),  // 12
                ],
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
                    // four faces from edge (0,1), (1,2), (2,3), (3,0):
                    9, 4, 10,
                    10, 5, 11,
                    11, 6, 12,
                    12, 7, 9,
                ],
            }),
            final_geom: Rc::new(prim::empty_mesh().to_meshfunc()),
            children: vec![
                child!(rule!(recur, ()), opening_xform(0), 5, 2, 6, 8),
                child!(rule!(recur, ()), opening_xform(1), 4, 1, 5, 8),
                child!(rule!(recur, ()), opening_xform(2), 7, 0, 4, 8),
                child!(rule!(recur, ()), opening_xform(3), 6, 3, 7, 8),
                // TODO: These vertex mappings appear to be right.
                // Explain *why* they are right.
            ],
        }
    };

    Rule { eval: Rc::new(start), ctxt: () }
}

/*

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
        let r = FRAC_PI_2 * i;
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
                arg_vals: vec![5,2,6,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(1.0),
                arg_vals: vec![4,1,5,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(2.0),
                arg_vals: vec![7,0,4,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(3.0),
                arg_vals: vec![6,3,7,8],
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
                        arg_vals: vec![0,1,2,3],
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
                    arg_vals: vec![0,1,2,3],
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
        let r = FRAC_PI_2 * i;
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
                arg_vals: vec![5,2,6,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(1.0),
                arg_vals: vec![4,1,5,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(2.0),
                arg_vals: vec![7,0,4,8],
            },
            Child {
                rule: Rc::new(Rule { eval: recur.clone(), ctxt }),
                xf: opening_xform(3.0),
                arg_vals: vec![6,3,7,8],
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
                        arg_vals: vec![0,1,2,3],
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
                    arg_vals: vec![0,1,2,3],
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
                    arg_vals: vec![0,1,2,3],
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
                PI).to_homogeneous(),
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
                    arg_vals: vec![0,1,2,3],
                },
                Child {
                    rule: Rule { eval: Rc::new(move || self.recur()) },
                    xf: self.flip180,
                    arg_vals: vec![3,2,1,0],
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
                    arg_vals: vec![0,1,2,3],
                },
            ],
        }
    }
}
*/

pub fn test_parametric() -> Mesh {

    let base_verts: Vec<Vertex> = vec![
        vertex(-1.0, -1.0, 0.0),
        vertex(-1.0,  1.0, 0.0),
        vertex( 1.0,  1.0, 0.0),
        vertex( 1.0, -1.0, 0.0),
    ];
    let base_verts = util::subdivide_cycle(&base_verts, 2);
    //let base_verts = util::subdivide_cycle(&base_verts, 16);

    let t0 = 0.0;
    let t1 = 15.0;
    let xform = |t: f32| -> Transform {
        id().
            translate(0.0, 0.0, t/5.0).
            rotate(&Vector3::z_axis(), -t/2.0).
            scale((0.8).powf(t))
    };

    crate::rule::parametric_mesh(base_verts, xform, t0, t1, 0.01)
}

pub fn test_dcel(fname: &str) {
    let mut mesh: dcel::DCELMesh<Vertex> = dcel::DCELMesh::new();
    let (f1, _) = mesh.add_face([
        VertSpec::New(vertex(-0.5, -0.5, 0.0)),
        VertSpec::New(vertex(-0.5,  0.5, 0.0)),
        VertSpec::New(vertex( 0.5,  0.5, 0.0)),
    ]);
    mesh.check();
    let (f2, edges) = mesh.add_face_twin1(mesh.faces[f1].halfedge, vertex(0.0, 0.0, 1.0));
    mesh.check();

    // From add_face_twin1, edges[0] is always the 'shared' edge:
    let edge = edges[0];
    let twin = {
        let he = &mesh.halfedges[edge];
        if he.has_twin {
            he.twin_halfedge
        } else {
            panic!("Can't find shared edge!");
        }
    };
    println!("Shared edges = {},{}", edge, twin);

    let ep = mesh.halfedges[edge].prev_halfedge;
    let en = mesh.halfedges[edge].next_halfedge;
    let tp = mesh.halfedges[twin].prev_halfedge;
    let tn = mesh.halfedges[twin].next_halfedge;
    println!("Connecting halfedges: {} and {}, {} and {}", en, tp, tn, ep);

    println!("DCEL mesh = {}", mesh);
    // As we're making *twin* halfedges, we go against the edge
    // direction:
    let (f3, _) = mesh.add_face_twin2(en, tp);
    mesh.check();
    let (f4, _) = mesh.add_face_twin2(tn, ep);
    mesh.check();

    println!("f1 verts: {:?}", mesh.face_to_verts(f1));
    println!("f2 verts: {:?}", mesh.face_to_verts(f2));
    println!("f3 verts: {:?}", mesh.face_to_verts(f3));
    println!("f4 verts: {:?}", mesh.face_to_verts(f4));

    //println!("DCEL mesh: ");
    //mesh.print();

    let faces = mesh.full_subdiv_face(f1, vec![
        vertex(-0.5, 0.0, 0.0),
        vertex(0.0, 0.5, 0.0),
        vertex(0.0, 0.0, 0.0),
    ]);
    println!("full_subdiv_face returned: {:?}", faces);

    //println!("DCEL mesh after subdiv");
    //mesh.check();
    //mesh.print();

    let mesh_conv = mesh.convert_mesh(|i| i);

    println!("Mesh = {:?}", mesh_conv);

    mesh_conv.write_stl_file(fname).unwrap();
}