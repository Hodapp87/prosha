use crate::openmesh::{OpenMesh, Tag, Mat4};
//use crate::prim;
use std::rc::Rc;

/// Definition of a rule.  In general, a `Rule`:
///
/// - produces geometry when it is evaluated
/// - tells what other rules to invoke, and what to do with their
/// geometry
pub struct Rule {
    pub eval: Box<dyn Fn(Rc<Rule>) -> RuleEval>,
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
pub struct RuleEval {
    /// The geometry generated at just this iteration
    pub geom: OpenMesh,

    /// The "final" geometry that is merged with `geom` via
    /// `connect()` in the event that recursion stops.  This must be
    /// in the same coordinate space as `geom`.
    ///
    /// Parent vertex references will be resolved directly to `geom`
    /// with no mapping.
    pub final_geom: OpenMesh,

    /// The child invocations (used if recursion continues).  The
    /// 'parent' mesh, from the perspective of all geometry produced
    /// by `children`, is `geom`.
    pub children: Vec<Child>,
}

/// `Child` evaluations, pairing another `Rule` with the
/// transformations and parent vertex mappings that should be applied
/// to it.
pub struct Child {

    /// Rule to evaluate to produce geometry
    pub rule: Rc<Rule>,

    /// The transform to apply to all geometry produced by `rule`
    /// (including its own `geom` and `final_geom` if needed, as well
    /// as all sub-geometry produced recursively).
    pub xf: Mat4,

    /// The parent vertex mapping: a mapping to apply to turn a
    /// Tag::Parent vertex reference into a vertex index of the parent
    /// mesh.  That is, if `rule` produces an `OpenMesh` with a face
    /// of `Tag::Parent(n)`, this will correspond to index `vmap[n]`
    /// in the parent mesh.
    pub vmap: Vec<usize>,
}

impl Rule {

    // TODO: Do I want to make 'geom' shared somehow, maybe with Rc? I
    // could end up having a lot of identical geometry that need not be
    // duplicated until it is transformed into the global space.
    //
    // This might produce bigger gains if I rewrite to_mesh so that
    // rather than repeatedly transforming meshes, it stacks
    // transformations and then applies them all at once.

    /// Convert this `Rule` to mesh data, recursively (depth first).
    /// `iters_left` sets the maximum recursion depth.  This returns
    /// (geometry, number of rule evaluations).
    pub fn to_mesh(s: Rc<Self>, iters_left: u32) -> (OpenMesh, usize) {

        let mut evals = 1;

        let rs: RuleEval = (s.eval)(s.clone());
        if iters_left <= 0 {
            return (rs.final_geom, 1);
            // TODO: This is probably wrong because of the way that
            // sub.vmap is used below.  final_geom is not supposed to
            // have any vertex mapping applied.
        }

        // TODO: This logic is more or less right, but it
        // could perhaps use some un-tupling or something.

        let subgeom: Vec<(OpenMesh, &Vec<usize>)> = rs.children.iter().map(|sub| {
            // Get sub-geometry (still un-transformed):
            let (submesh, eval) = Rule::to_mesh(sub.rule.clone(), iters_left - 1);
            // Tally up eval count:
            evals += eval;
            
            let m2 = submesh.transform(&sub.xf);
            
            (m2, &sub.vmap)
        }).collect();
        
        // Connect geometry from this rule (not child rules):
        return (rs.geom.connect(&subgeom).0, evals);
    }

    /// This should be identical to to_mesh, but implemented
    /// iteratively with an explicit stack rather than with recursive
    /// function calls.
    pub fn to_mesh_iter(s: Rc<Self>, max_depth: usize) -> (OpenMesh, usize) {

        struct State {
            // The set of rules we're currently handling:
            rules: Vec<Child>,
            // The next element of 'children' to handle:
            next: usize,
            // World transform of the *parent* of 'rules', that is,
            // not including any transform of any element of 'rules'.
            xf: Mat4,
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
        let mut stack: Vec<State> = vec![State {
            rules: eval.children,
            next: 0,
            xf: nalgebra::geometry::Transform3::identity().to_homogeneous(),
            depth: max_depth,
        }];
        let mut geom = eval.geom;

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
                // generate an identity mapping for vmap when I don't
                // actually need vmap.
                let m = {
                    let mut m_ = 0;
                    for f in &final_geom.faces {
                        match f {
                            Tag::Parent(i) => {
                                if *i > m_ {
                                    m_ = *i;
                                }
                            }
                            _ => {}
                        }
                    }
                    m_ + 1
                };
                let vmap: Vec<usize> = (0..m).collect();
                let (geom2, _) = new_geom.connect(&vec![(final_geom, &vmap)]);
                
                geom = geom.connect(&vec![(geom2, &child.vmap)]).0;
                
                // and backtrack:
                stack.pop();
                n -= 1;
                continue;
            }

            let (g, offsets) = geom.connect(&vec![(new_geom, &child.vmap)]);
            geom = g;

            // 'new_geom' may itself be parent geometry for
            // 'eval.children' (via Tag::Parent), and vmap is there to
            // resolve Tag::Parent references to the right vertices in
            // 'new_geom'.
            //
            // However, we connect() on the global geometry which we
            // merged 'new_geom' into, not 'new_geom' directly.  To
            // account for this, we must shift vmap by the offset that
            // 'geom.connect' gave us:
            for (offset, child) in offsets.iter().zip(eval.children.iter_mut()) {
                child.vmap = child.vmap.iter().map(|n| {
                    n + offset
                }).collect();
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
