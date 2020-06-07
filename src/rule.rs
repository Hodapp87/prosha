use std::borrow::Borrow;
use std::rc::Rc;
use std::f32;

use crate::mesh::{Mesh, MeshFunc, VertexUnion};
use crate::xform::{Transform, Vertex};
use crate::dcel;
use crate::dcel::{DCELMesh, VertSpec};

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

    #[derive(Copy, Clone, Debug)]
    struct VertexTrajectory {
        // Vertex position
        vert: Vertex,
        // Parameter value; f(t) should equal vert
        t: f32,
        // "Starting" vertex position, i.e. at f(t0). Always either a frame
        // vertex, or a linear combination of two neighboring ones.
        frame_vert: Vertex,
    };

    let mut mesh: DCELMesh<VertexTrajectory> = DCELMesh::new();

    #[derive(Clone, Debug)]
    struct frontierVert {
        traj: VertexTrajectory,
        // If the boundaries on either side of this vertex lie on a face
        // (which is the case for all vertices *except* the initial ones),
        // then this gives the halfedges of those boundaries. halfedges[0]
        // connects the 'prior' vertex on the frontier to this, and
        // halfedges[1] connect this to the 'next' vertex on the fronter.
        // (Direction matters. If halfedges[0] is given, it must *end* at
        // 'traj.vert'. If halfedges[1] is given, it must *begin* at 'traj.vert'.)
        halfedges: [Option<usize>; 2],
        // If this vertex is already in 'mesh', its vertex index there:
        vert_idx: Option<usize>,
    };

    // Init 'frontier' with each 'frame' vertex, and start it at t=t0.
    let mut frontier: Vec<frontierVert> = frame.iter().enumerate().map(|(i,v)| frontierVert {
        traj: VertexTrajectory {
            vert: *v,
            t: t0,
            frame_vert: *v,
        },
        halfedges: [None; 2],
        vert_idx: None,
    }).collect();
    // Every vertex in 'frontier' has a trajectory it follows - which is
    // simply the position as we transform the original vertex by f(t),
    // and increment t through the range [t0, t1].
    //
    // The goal is to advance the vertices, one at a time, building up
    // new triangles every time we advance one, until each vertex
    // reaches t=t1 - in a way that forms the mesh we want.

    while !frontier.is_empty() {
        /*
        for (i, f) in frontier.iter().enumerate() {
            println!("DEBUG: frontier[{}]: vert={},{},{} vert_idx={:?} t={} halfedges={:?}", i, f.vert.x, f.vert.y, f.vert.z, f.vert_idx, f.t, f.halfedges);
            match f.vert_idx {
                Some(vert) => {
                    match f.halfedges[1] {
                        Some(idx) => {
                            if mesh.halfedges[idx].vert != vert {
                                println!("  Error: halfedge[1]={} starts at {}, should start at {}", idx, mesh.halfedges[idx].vert, vert);
                            }
                        },
                        None => {},
                    };
                    match f.halfedges[0] {
                        Some(idx) => {
                            let v2 = mesh.halfedges[mesh.halfedges[idx].next_halfedge].vert;
                            if v2 != vert {
                                println!("  Error: halfedge[0]={} ends at {}, should start at {}", idx, v2, vert);
                            }
                        },
                        None => {},
                    };
                },
                None => {},
            };
        }
         */

        // Pick a vertex to advance.  For now, this is an easy
        // heuristic: pick the 'furthest back' one (lowest t value).
        let (i, v) = {

            let mut t_min = t1;
            let mut i_min = 0;
            let mut all_done = true;

            for (i, v) in frontier.iter().enumerate() {
                if v.traj.t < t_min {
                    i_min = i;
                    t_min = v.traj.t;
                }
                all_done &= v.traj.t >= t1;
            }

            // If no vertex can be advanced, we're done:
            if all_done {
                println!("All past t1");
                break;
            }
            (i_min, frontier[i_min].clone())
        };

        // Indices (into 'frontier') of previous and next:
        let i_prev = (i + n - 1) % n;
        let i_next = (i + 1) % n;

        println!("DEBUG: Moving frontier vertex {}, {:?} (t={}, frame_vert={:?})", i, v.traj.vert, v.traj.t, v.traj.frame_vert);

        // Move this vertex further along, i.e. t + dt.  (dt is set by
        // the furthest we can go while remaining within 'err', i.e. when we
        // make our connections we look at how far points on the *edges*
        // diverge from the trajectory of the  continuous transformation).
        let vf = v.traj.frame_vert;
        let vt = v.traj.t;
        let dt_max = t1 - vt;
        let mut dt = ((t1 - t0) / 100.0).min(dt_max);
        for iter in 0..100 {
            // Consider an edge from f(v.t)*vf to f(v.t + dt)*vf.
            // These two endpoints have zero error from the trajectory
            // (because they are directly on it).
            //
            // If we assume some continuity in f, then we can guess that
            // the worst error occurs at the midpoint of the edge:
            let edge_mid = 0.5 * (f(vt).mtx + f(vt + dt).mtx) * vf;
            // ...relative to the trajectory midpoint:
            let traj_mid = f(vt + dt / 2.0).mtx * vf;
            let err = (edge_mid - traj_mid).norm();

            //println!("DEBUG iter {}: dt={}, edge_mid={:?}, traj_mid={:?}, err={}", iter, dt, edge_mid, traj_mid, err);

            let r = (err - max_err).abs() / max_err;
            if r < 0.10 {
                //println!("err close enough");
                println!("Error under threshold in {} iters, dt={}", iter, dt);
                break;
            } else if err > max_err {
                dt = dt / 2.0;
                //println!("err > max_err, reducing dt to {}", dt);
            } else {
                dt = (dt * 1.2).min(dt_max);
                //println!("err < max_err, increasing dt to {}", dt);
            }
        }
        // (note that this is just a crappy numerical approximation of
        // dt given a desired error and it would probably be possible
        // to do this directly given an analytical form of the
        // curvature of f at some starting point)

        let t = (vt + dt).min(t1);
        let v_next = VertexTrajectory {
            vert: f(t).mtx * vf,
            t: t,
            frame_vert: vf,
        };

        // DEBUG
        /*
        let l1 = (v.vert - v_next).norm();
        let l2 = (frontier[v.neighbor1].vert - v.vert).norm();
        let l3 = (frontier[v.neighbor2].vert - v.vert).norm();
        println!("l1={} l2={} l3={}", l1, l2, l3);
         */

        // Add two faces to our mesh. They share two vertices, and thus
        // the boundary in between those vertices.

        //println!("DEBUG: Adding 1st face");

        // First face: connect prior frontier vertex to 'v' & 'v_next'.
        // f1 is that face,
        // edge1 is: prior frontier vertex -> v_next.
        // edge_v_next is: v_next -> v
        let (f1, edge1, edge_v_next) = match v.halfedges[0] {
            // However, the way we add the face depends on whether we are
            // adding to an existing boundary or not:
            None => {
                let neighbor = &frontier[i_prev];
                //println!("DEBUG: add_face()");
                let (f1, edges) = mesh.add_face([
                    VertSpec::New(v_next), // edges[0]: v_next -> v
                    match v.vert_idx {
                        None => VertSpec::New(v.traj),
                        Some(idx) => VertSpec::Idx(idx),
                    }, // edges[1]: v -> neighbor
                    match neighbor.vert_idx {
                        None => VertSpec::New(neighbor.traj),
                        Some(idx) => VertSpec::Idx(idx),
                    }, // edges[2]: neighbor -> v_next
                ]);

                if neighbor.vert_idx.is_none() {
                    // If neighbor.vert_idx is None, then we had to
                    // add its vertex to the mesh for the face we just
                    // made - so mark it in the frontier:
                    frontier[i_prev].vert_idx = Some(mesh.halfedges[edges[2]].vert);
                }

                (f1, edges[2], edges[0])
            },
            Some(edge_idx) => {
                // edge_idx is half-edge from prior frontier vertex to 'v'
                //println!("DEBUG: add_face_twin1({}, ({},{},{}))", edge_idx, v_next.x, v_next.y, v_next.z);
                let (f1, edges) = mesh.add_face_twin1(edge_idx, v_next);
                // edges[0] is twin of edge_idx, thus v -> neighbor
                // edges[1] is neighbor -> v_next
                // edges[2] is v_next -> v
                (f1, edges[1], edges[2])
            },
        };
        //println!("DEBUG: edge1={} edge_v_next={}", edge1, edge_v_next);

        // DEBUG
        //mesh.check();
        //mesh.print();

        //println!("DEBUG: Adding 2nd face");
        // Second face: connect next frontier vertex to 'v' & 'v_next'.
        // f2 is that face.
        // edge2 is: v_next -> next frontier vertex
        let (f2, edge2) = match v.halfedges[1] {
            // Likewise, the way we add the second face depends on
            // the same (but for the other side).  Regardless,
            // they share the boundary between v_next and v.vert - which
            // is edges1[0].
            None => {
                let neighbor = &frontier[i_next];

                let (f2, edges) = match neighbor.vert_idx {
                    None => {
                        //let v = neighbor.traj.vert;
                        //println!("DEBUG: add_face_twin1({}, ({},{},{}))", edge_v_next, v.x, v.y, v.z);
                        let (f2, edges) = mesh.add_face_twin1(edge_v_next, neighbor.traj);
                        // Reasoning here is identical to "If neighbor.vert_idx
                        // is None..." above:
                        frontier[i_next].vert_idx = Some(mesh.halfedges[edges[2]].vert);

                        (f2, edges)
                    },
                    Some(vert_idx) => {
                        //println!("DEBUG: add_face_twin1_ref({}, {})", edge_v_next, vert_idx);
                        mesh.add_face_twin1_ref(edge_v_next, vert_idx)
                    }
                };
                // For either branch:
                // edges[0] is twin of edge_v_next, thus: v -> v_next
                // edges[1] is: v_next -> neighbor
                // egdes[2] is: neighbor -> v

                (f2, edges[1])
            },
            Some(edge_idx) => {
                //println!("DEBUG: add_face_twin2({}, {})", edge_idx, edge_v_next);
                let (f2, edges) = mesh.add_face_twin2(edge_idx, edge_v_next);
                // edges[0] is twin of edge_idx, thus: neighbor -> v
                // edges[1] is twin of edge_v_next, thus: v -> v_next
                // edges[2] is: v_next -> neighbor
                //println!("DEBUG: add_face_twin2 returned {:?}", edges);
                (f2, edges[2])
            },
        };
        //println!("DEBUG: edge2={}", edge2);

        // DEBUG
        //println!("DEBUG: Added 2nd face");
        //mesh.check();
        //mesh.print();

        // The 'shared' half-edge should start at v.vert - hence edges[1].

        // and add two faces:
        /*
        faces.append(&mut vec![
            v_next_idx, v.mesh_idx, frontier[v.neighbor[0]].mesh_idx, // face_idx
            v.mesh_idx, v_next_idx, frontier[v.neighbor[1]].mesh_idx, // face_idx + 1
        ]);
         */

        // Replace this vertex in the frontier:
        frontier[i] = frontierVert {
            traj: VertexTrajectory {
                vert: v_next.vert,
                frame_vert: vf,
                t: t,
            },
            halfedges: [Some(edge1), Some(edge2)],
            vert_idx: Some(mesh.halfedges[edge_v_next].vert),
        };
        // Since 'halfedges' forms more or less a doubly-linked list, we
        // must go to previous & next vertices and update them:
        frontier[i_prev].halfedges[1] = Some(edge1);
        frontier[i_next].halfedges[0] = Some(edge2);
    }

    // A stack of face indices for faces in 'mesh' that are still
    // undergoing subdivision
    let mut stack: Vec<usize> = (0..mesh.num_faces).collect();

    while !stack.is_empty() {
        let face = stack.pop().unwrap();
        println!("DEBUG: Examining face: {:?}", face);

        let v_idx = mesh.face_to_verts(face);
        if v_idx.len() != 3 {
            panic!(format!("Face {} has {} vertices, not 3?", face, v_idx.len()));
        }
        let tr = [mesh.verts[v_idx[0]].v,
                         mesh.verts[v_idx[1]].v,
                         mesh.verts[v_idx[2]].v];
        let (v0, v1, v2) = (tr[0].vert, tr[1].vert, tr[2].vert);
        let d01 = (v0 - v1).norm();
        let d02 = (v0 - v2).norm();
        let d12 = (v1 - v2).norm();
        // Law of cosines:
        let cos0 = (d01*d01 + d02*d02 - d12*d12) / (2.0 * d01 * d02);
        let cos1 = (d01*d01 + d12*d12 - d02*d02) / (2.0 * d01 * d12);
        let cos2 = (d02*d02 + d12*d12 - d01*d01) / (2.0 * d02 * d12);
        // TODO: Perhaps use https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates
        // to go by circumradius instead. Or - find it by law of sines?
        println!("DEBUG: d01={} d02={} d12={} cos0={} cos1={} cos2={}", d01, d02, d12, cos0, cos1, cos2);
        if (cos0 < 0.0 || cos0 > 0.7 ||
            cos1 < 0.0 || cos1 > 0.7 ||
            cos2 < 0.0 || cos2 > 0.7) {
            println!("DEBUG: Angles out of range!");
        } else {
            println!("DEBUG: Angles OK");
        }
        // TODO: Figure out how to subdivide in this case

        // The triangle forms a plane.  Get this plane's normal vector.
        let a = (v0 - v1).xyz();
        let b = (v0 - v2).xyz();
        let normal = a.cross(&b).normalize();
        // Make a new point that is on the surface, but roughly a
        // midpoint in parameter space. Exact location isn't crucial.
        let t_mid = (tr[0].t + tr[1].t + tr[2].t) / 3.0;
        let v_mid = (tr[0].frame_vert + tr[1].frame_vert + tr[2].frame_vert) / 3.0;
        let p = f(t_mid).mtx * v_mid;

        let d = p.xyz().dot(&normal);

        println!("DEBUG: t_mid={} v_mid={},{},{} p={},{},{}", t_mid, v_mid.x, v_mid.y, v_mid.z, p.x, p.y, p.z);
        println!("DEBUG: d={}", d);

        // DEBUG
        /*
        let n = verts.len();
        verts.push(p);
        faces.push(face.verts[0]);
        faces.push(face.verts[1]);
        faces.push(n);
         */
        if (d <= max_err) {
            // This triangle is already in the mesh, and already popped
            // off of the stack. We're done.
            continue;
        }
        println!("DEBUG: d > err, splitting this triangle");

        // The face has 3 edges.  We split each of them (making a new
        // vertex in the middle), and then make 3 new edges - one
        // between each pair of new vertices - to replace the face
        // with four smaller faces.

        // This split is done in 'parameter' space:
        let pairs = [(0,1), (1,2), (0,2)];
        let mut mids: Vec<Vertex> = pairs.iter().map(|(i,j)| {
            let t = (tr[*i].t + tr[*j].t) / 2.0;
            let v = (tr[*i].frame_vert + tr[*j].frame_vert) / 2.0;
            f(t).mtx * v
        }).collect();

        // DEBUG
        //let n = verts.len();
        // Index n+0 is (0,1), n+1 is (1,2), n+2 is (0,2)
        /*
        verts.append(&mut mids);
        faces[face.face] = n + 0;
        faces[face.face + 1] = n + 1;
        faces[face.face + 2] = n + 2;
        faces.extend_from_slice(&[
            face.verts[0], n + 0, n + 2,
            face.verts[1], n + 1, n + 0,
            face.verts[2], n + 2, n + 1,
            //n + 0,         n + 1, n + 2,
        ]);
        */
   }

    return mesh.convert_mesh(|i| i.vert );
}
