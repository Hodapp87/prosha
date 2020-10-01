use std::borrow::Borrow;
use std::rc::Rc;
use std::f32;

use crate::mesh::{Mesh, MeshFunc, VertexUnion};
use crate::xform::{Transform, Vertex};

pub type RuleFn<S> = Rc<dyn Fn(Rc<Rule<S>>) -> RuleEval<S>>;

/// Definition of a rule.  In general, a `Rule`:
///
/// - produces geometry when it is evaluated
/// - tells what other rules to invoke, and what to do with their
/// geometry
pub struct Rule<S> {
    pub eval: RuleFn<S>,
    pub ctxt: S,
}

// TODO: It may be possible to have just a 'static' rule that requires
// no function call.
// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?

// TODO: Why *can't* I make this FnOnce?

// The above looks like it is going to require a lifetime parameter
// regardless, in which case I don't really need Box.

/// `RuleEval` supplies the results of evaluating some `Rule` for one
/// iteration: it contains the geometry produced at this step
/// (`geom`), and it tells what to do next depending on whether
/// recursion continues further, or is stopped here (due to hitting
/// some limit of iterations or some lower limit on overall scale).
///
/// That is:
/// - if recursion stops, `final_geom` is connected with `geom`.
/// - if recursion continues, the rules of `children` are evaluated,
/// and the resultant geometry is transformed and then connected with
/// `geom`.
pub struct RuleEval<S> {
    /// The geometry generated at just this iteration
    pub geom: Rc<MeshFunc>,

    /// The "final" geometry that is merged with `geom` via
    /// `connect()` in the event that recursion stops.  This must be
    /// in the same coordinate space as `geom`.
    ///
    /// Parent vertex references will be resolved directly to `geom`
    /// with no mapping.
    /// (TODO: Does this make sense? Nowhere else do I treat Arg(n) as
    /// an index - it's always a positional argument.)
    pub final_geom: Rc<MeshFunc>,

    /// The child invocations (used if recursion continues).  The
    /// 'parent' mesh, from the perspective of all geometry produced
    /// by `children`, is `geom`.
    pub children: Vec<Child<S>>,
}

/// `Child` evaluations, pairing another `Rule` with the
/// transformations and parent vertex mappings that should be applied
/// to it.
pub struct Child<S> {

    /// Rule to evaluate to produce geometry
    pub rule: Rc<Rule<S>>,

    /// The transform to apply to all geometry produced by `rule`
    /// (including its own `geom` and `final_geom` if needed, as well
    /// as all sub-geometry produced recursively).
    pub xf: Transform,

    /// The 'argument values' to apply to vertex arguments of a `MeshFunc`
    /// from `geom` and `final_geom` that `rule` produces when evaluated.
    /// The values of this are treated as indices into the parent
    /// `RuleEval` that produced this `Child`.
    ///
    /// In specific: if `arg_vals[i] = j` and `rule` produces some `geom` or
    /// `final_geom`, then any vertex of `VertexUnion::Arg(i)` will be mapped
    /// to `geom.verts[j]` in the *parent* geometry.
    pub arg_vals: Vec<usize>,
}

#[macro_export]
macro_rules! child {
    ( $Rule:expr, $Xform:expr, $( $Arg:expr ),* ) => {
        Child {
            rule: /*std::rc::Rc::new*/($Rule),
            xf: $Xform,
            arg_vals: vec![$($Arg,)*],
        }
    }
}

#[macro_export]
macro_rules! child_iter {
    ( $Rule:expr, $Xform:expr, $Args:expr ) => {
        Child {
            rule: /*std::rc::Rc::new*/($Rule),
            xf: $Xform,
            arg_vals: $Args.collect(), // does this even need a macro?
        }
    }
}

#[macro_export]
macro_rules! rule {
    ( $RuleFn:expr, $Ctxt:expr ) => {
        std::rc::Rc::new(Rule {
            eval: $RuleFn.clone(),
            ctxt: $Ctxt,
        })
    }
}

#[macro_export]
macro_rules! rule_fn {
    ( $Ty:ty => |$Self:ident $(,$x:ident)*| $Body:expr ) => {
        {
            $(let $x = $x.clone();)*
            std::rc::Rc::new(move |$Self: std::rc::Rc<Rule<$Ty>>| -> RuleEval<$Ty> {
                $(let $x = $x.clone();)*
                let $Self = $Self.clone();
                $Body
            })
        }
    }
}
// TODO: Shouldn't I fully-qualify Rule & RuleEval?
// TODO: Document all of the above macros
// TODO: Why must I clone twice?

impl<S> Rule<S> {

    /// Convert this `Rule` to mesh data, recursively (depth first).
    /// `iters_left` sets the maximum recursion depth.  This returns
    /// (geometry, number of rule evaluations).
    pub fn to_mesh(s: Rc<Rule<S>>, iters_left: usize) -> (MeshFunc, usize) {

        let mut evals = 1;

        let rs: RuleEval<S> = (s.eval)(s.clone());
        if iters_left <= 0 {
            return ((*rs.final_geom).clone(), 1);
            // TODO: This is probably wrong because of the way that
            // sub.arg_vals is used below.  final_geom is not supposed to
            // have any vertex mapping applied.
        }

        // TODO: This logic is more or less right, but it
        // could perhaps use some un-tupling or something.

        let subgeom: Vec<(MeshFunc, Vec<usize>)> = rs.children.iter().map(|sub| {
            // Get sub-geometry (still un-transformed):
            let (submesh, eval) = Rule::to_mesh(sub.rule.clone(), iters_left - 1);
            // Tally up eval count:
            evals += eval;
            
            let m2 = submesh.transform(&sub.xf);
            
            (m2, sub.arg_vals.clone())
                // TODO: Fix clone?
        }).collect();
        
        // Connect geometry from this rule (not child rules):
        return (rs.geom.connect(subgeom).0, evals);
    }

    /// This should be identical to to_mesh, but implemented
    /// iteratively with an explicit stack rather than with recursive
    /// function calls.
    pub fn to_mesh_iter(s: Rc<Rule<S>>, max_depth: usize) -> (MeshFunc, usize) {

        struct State<S> {
            // The set of rules we're currently handling:
            rules: Vec<Child<S>>,
            // The next element of 'children' to handle:
            next: usize,
            // World transform of the *parent* of 'rules', that is,
            // not including any transform of any element of 'rules'.
            xf: Transform,
            // How many levels 'deeper' can we recurse?
            depth: usize,
        }

        // 'stack' stores at its last element our "current" State in
        // terms of a current world transform and which Child should
        // be processed next.  Every element prior to this is previous
        // states which must be kept around for further backtracking
        // (usually because they involve multiple rules).
        //
        // We evaluate our own rule to initialize the stack:
        let eval = (s.eval)(s.clone());
        let mut stack: Vec<State<S>> = vec![State {
            rules: eval.children,
            next: 0,
            xf: Transform::new(),
            depth: max_depth,
        }];
        let mut geom = (*eval.geom).clone();

        // Number of times we've evaluated a Rule:
        let mut eval_count = 1;

        // Stack depth (update at every push & pop):
        let mut n = stack.len();

        while !stack.is_empty() {

            // s = the 'current' state:
            let s = &mut stack[n-1];
            let depth = s.depth;

            if s.rules.is_empty() {
                stack.pop();
                n -= 1;
                continue;
            }
            
            // Evaluate the rule:
            let child = &s.rules[s.next];
            let mut eval = (child.rule.eval)(child.rule.clone());
            eval_count += 1;

            // Make an updated world transform:
            let xf = s.xf * child.xf;

            // This rule produced some geometry which we'll
            // combine with the 'global' geometry:
            let new_geom = eval.geom.transform(&xf);
            
            // See if we can still recurse further:
            if depth <= 0 {
                // As we're stopping recursion, we need to connect
                // final_geom with all else in order to actually close
                // geometry properly:
                let final_geom = eval.final_geom.transform(&xf);
                // TODO: Fix the awful hack below.  I do this only to
                // generate an identity mapping for arg_vals when I don't
                // actually need arg_vals.
                let m = {
                    let mut m_ = 0;
                    for v in &final_geom.verts {
                        match *v {
                            VertexUnion::Arg(a) => {
                                if a > m_ {
                                    m_ = a;
                                }
                            },
                            VertexUnion::Vertex(_) => (),
                        }
                    }
                    m_ + 1
                };
                let arg_vals: Vec<usize> = (0..m).collect();
                let (geom2, _) = new_geom.connect(vec![(final_geom, arg_vals)]);
                
                geom = geom.connect(vec![(geom2, child.arg_vals.clone())]).0;
                // TODO: Fix clone?

                // If we end recursion on one child, we must end it
                // similarly on every sibling (i.e. get its geometry &
                // final geometry, and merge it in) - so we increment
                // s.next and let the loop re-run.
                s.next += 1;
                if s.next >= s.rules.len() {
                    // Backtrack only at the last child:
                    stack.pop();
                    n -= 1;
                }
                continue;
            }

            let (g, offsets) = geom.connect(vec![(new_geom, child.arg_vals.clone())]);
            geom = g;

            // 'eval.children' may contain (via 'arg_vals') references to
            // indices of 'new_geom'. However, we don't connect() to
            // 'new_geom', but to the global geometry we just merged it
            // into.  To account for this, we must shift 'arg_vals' by
            // the offset that 'geom.connect' gave us.
            let off = offsets[0];
            // (We pass a one-element vector to geom.connect() above
            // so offsets always has just one element.)
            for child in eval.children.iter_mut() {
                child.arg_vals = child.arg_vals.iter().map(|n| n + off).collect();
            }

            // We're done evaluating this rule, so increment 'next'.
            // If that was the last rule at this level (i.e. ignoring
            // eval.children), remove it - we're done with it.
            s.next += 1;
            if s.next >= s.rules.len() {
                stack.pop();
                n -= 1;
            }

            if !eval.children.is_empty() {
                // Recurse further (i.e. put more onto stack):
                stack.push(State {
                    rules: eval.children,
                    next: 0,
                    xf: xf,
                    depth: depth - 1,
                });
                n += 1;
            }
        }
        
        return (geom, eval_count); 
    }
    
}

impl<S> RuleEval<S> {
    /// Turn an iterator of (MeshFunc, Child) into a single RuleEval.
    /// All meshes are merged, and the `arg_vals` in each child has the
    /// correct offsets applied to account for this merge.
    ///
    /// (`final_geom` is passed through to the RuleEval unmodified.)
    pub fn from_pairs<T, U>(m: T, final_geom: MeshFunc) -> RuleEval<S>
        where U: Borrow<MeshFunc>,
              T: IntoIterator<Item = (U, Child<S>)>
    {
        let (meshes, children): (Vec<_>, Vec<_>) = m.into_iter().unzip();
        let (mesh, offsets) = MeshFunc::append(meshes);

        // Patch up arg_vals in each child, and copy everything else:
        let children2: Vec<Child<S>> = children.iter().zip(offsets.iter()).map(|(c,off)| {
            Child {
                rule: c.rule.clone(),
                xf: c.xf.clone(),
                // simply add offset:
                arg_vals: c.arg_vals.iter().map(|i| i + off).collect(),
            }
        }).collect();

        RuleEval {
            geom: Rc::new(mesh),
            final_geom: Rc::new(final_geom),
            children: children2,
        }
    }

}

// fa001f47d40de989da6963e442f31c278c88abc8

/// Produce a mesh from a starting frame, and a function `f` which produces
/// transformations that change continuously over its argument (the range
/// of which is given by `t0` and `t1`).  By convention, `f(t0)` should
/// always produce an identity transformation.
///
/// Facetization is guided by the given error, `max_err`, which is treated
/// as a distance in 3D space.
pub fn parametric_mesh<F>(frame: Vec<Vertex>, f: F, t0: f32, t1: f32, max_err: f32) -> Mesh
    where F: Fn(f32) -> Transform
{
    let n = frame.len();

    // Sanity checks:
    if t1 <= t0 {
        panic!("t1 must be greater than t0");
    }
    if n < 3 {
        panic!("frame must have at least 3 vertices");
    }

    struct frontierVert {
        vert: Vertex, // Vertex position
        t: f32, // Parameter value; f(t) should equal vert
        frame_idx: usize, // Index of 'frame' this sits in the trajectory of
        mesh_idx: usize, // Index of this vertex in the mesh
        neighbor1: usize, // Index of 'frontier' of one neighbor
        neighbor2: usize, // Index of 'frontier' of other neighbor
    };

    // Init 'frontier' with each 'frame' vertex, and start it at t=t0.
    let mut frontier: Vec<frontierVert> = frame.iter().enumerate().map(|(i,v)| frontierVert {
        vert: *v,
        t: t0,
        frame_idx: i,
        mesh_idx: i,
        neighbor1: (i - 1) % n,
        neighbor2: (i + 1) % n,
    }).collect();
    // Every vertex in 'frontier' has a trajectory it follows - which is
    // simply the position as we transform the original vertex by f(t),
    // and increment t through the range [t0, t1].
    //
    // The goal is to advance the vertices, one at a time, building up
    // new triangles every time we advance one, until each vertex
    // reaches t=t1 - in a way that forms the mesh we want.

    // That mesh will be built up here, starting with frame vertices:
    // (note initial value of mesh_idx)
    let mut verts: Vec<Vertex> = frame.clone();
    let mut faces: Vec<usize> = vec![];

    while !frontier.is_empty() {

        // Pick a vertex to advance.
        //
        // Heuristic for now: pick the 'furthest back' (lowest t)
        let (i,v) = frontier.iter().enumerate().min_by(|(i,f), (j, g)|
            f.t.partial_cmp(&g.t).unwrap_or(std::cmp::Ordering::Equal)).unwrap();
        // TODO: Make this less ugly?

        if v.t >= t1 {
            break;
        }

        println!("DEBUG: Moving vertex {}, {:?} (t={}, frame_idx={})", i, v.vert, v.t, v.frame_idx);

        let mut dt = (t1 - t0) / 100.0;
        let vf = frame[v.frame_idx];
        for iter in 0..100 {
            // Consider an edge from f(v.t)*vf to f(v.t + dt)*vf.
            // These two endpoints have zero error from the trajectory
            // (because they are directly on it).
            //
            // If we assume some continuity in f, then we can guess that
            // the worst error occurs at the midpoint of the edge:
            let edge_mid = 0.5*(f(v.t).mtx + f(v.t + dt).mtx)*vf;
            // ...relative to the trajectory midpoint:
            let traj_mid = f(v.t + dt/2.0).mtx * vf;
            let err = (edge_mid - traj_mid).norm();

            println!("DEBUG iter {}: dt={}, edge_mid={:?}, traj_mid={:?}, err={}", iter, dt, edge_mid, traj_mid, err);

            let r = (err - max_err).abs() / max_err;
            if r < 0.10 {
                println!("err close enough");
                break;
            } else if err > max_err {
                dt = dt / 2.0;
                println!("err > max_err, reducing dt to {}", dt);
            } else {
                dt = dt * 1.2;
                println!("err < max_err, increasing dt to {}", dt);
            }
        }

        let t = v.t + dt;
        let v_next = f(t).mtx * vf;

        // Add this vertex to our mesh:
        let pos = verts.len();
        verts.push(v_next);
        // There are 3 other vertices of interest: the one we started
        // from (v) and its two neighbors. We make two edges - one on
        // each side of the edge (v, v_next).
        faces.append(&mut vec![
            v.mesh_idx, pos, frontier[v.neighbor1].mesh_idx,
            pos, v.mesh_idx, frontier[v.neighbor2].mesh_idx,
        ]);

        // Replace this vertex in the frontier:
        frontier[i] = frontierVert {
            vert: v_next,
            frame_idx: v.frame_idx,
            mesh_idx: pos,
            t: t,
            neighbor1: v.neighbor1,
            neighbor2: v.neighbor2,
        }
    }

    // Move this vertex further along, i.e. t + dt.  (dt is set by
    // the furthest we can go while remaining within 'err', i.e. when we
    // make our connections we look at how far points on the *edges*
    // diverge from the trajectory of the  continuous transformation).

    // Add this vertex to the mesh, and connect it to: the vertex we
    // started with, and the two neighbors of that vertex.

    // Repeat at "Pick a vertex...".

    // Don't move t + dt past t1.  Once a frontier vertex is placed at
    // that value of t, remove it.


    // Missing: Anything about when to subdivide an edge.
    // If I assume a good criterion of "when" to subdivide an edge, the
    // "how" is straightforward: find the edge's two neighbors in the
    // frontier. Trace them back to their 'original' vertices at t=t0
    // (these should just be stored alongside each frontier member),
    // produce an interpolated vertex. Produce an interpolated t from
    // respective t of the two neighbors in the frontier; use that t
    // to move the 'interpolated' vertex along its trajectory.
    //
    // Add new vertex to mesh (and make the necessary connections)
    // and to frontier.

    // But still missing from that: When do I collapse a subdivision
    // back down?
    return Mesh { verts, faces };
}
