use crate::openmesh::{OpenMesh, Mat4};
use crate::prim;

/// Definition of a rule.  In general, a `Rule`:
///
/// - produces geometry when it is evaluated
/// - tells what other rules to invoke, and what to do with their
/// geometry
pub enum Rule<A> {
    /// Produce some geometry, and possibly recurse further.
    Recurse(fn (&A) -> RuleEval<A>),
    /// Produce nothing and recurse no further.
    EmptyRule,
}
// TODO: Rename rules?
// TODO: It may be possible to have just a 'static' rule that requires
// no function call.
// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?

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
pub struct RuleEval<A> {
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
    pub children: Vec<Child<A>>,
}

/// `Child` evaluations, pairing another `Rule` with the
/// transformations and parent vertex mappings that should be applied
/// to it.
pub struct Child<A> {

    /// Rule to evaluate to produce geometry
    pub rule: Rule<A>,

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

impl<A> Rule<A> {

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
    pub fn to_mesh(&self, arg: &A, iters_left: u32) -> (OpenMesh, u32) {

        let mut evals: u32 = 1;

        if iters_left <= 0 {
            match self {
                Rule::Recurse(f) => {
                    let rs: RuleEval<A> = f(arg);
                    return (rs.final_geom, 1);
                }
                Rule::EmptyRule => {
                    return (prim::empty_mesh(), evals);
                }
            }
        }

        match self {
            Rule::Recurse(f) => {
                let rs: RuleEval<A> = f(arg);
                // TODO: This logic is more or less right, but it
                // could perhaps use some un-tupling or something.

                let subgeom: Vec<(OpenMesh, &Vec<usize>)> = rs.children.iter().map(|sub| {
                    // Get sub-geometry (still un-transformed):
                    let (submesh, eval) = sub.rule.to_mesh(arg, iters_left - 1);
                    // Tally up eval count:
                    evals += eval;
                    
                    let m2 = submesh.transform(&sub.xf);
                    
                    (m2, &sub.vmap)
                }).collect();
                
                // Connect geometry from this rule (not child rules):
                return (rs.geom.connect(&subgeom), evals);
            }
            Rule::EmptyRule => {
                return (prim::empty_mesh(), evals);
            }
        }
    }

    pub fn to_mesh_iter(&self, arg: &A, max_depth: usize) -> (OpenMesh, u32) {

        let mut geom = prim::empty_mesh();
        let mut stack: Vec<(RuleEval<A>, usize)> = vec![];

        match self {
            Rule::Recurse(f) => stack.push((f(arg), 0)),
            Rule::EmptyRule => {},
        }

        loop {
            let n = stack.len(); // TODO: Just keep a running total.
            // We can increment/decrement as we push/pop.
            let (eval, idx) = stack[n-1];
            // note that:
            // stack[n-2].children[idx] = eval    (so to speak)

            // I can't do the above. I can either...
            // - Implement Copy for RuleEval.
            // - push/pop all over the place and deal with the Option
            // unpacking and the extra state-changes.
            // - use something nicer than RuleEval?

            // I can't borrow it because a mutable borrow is already
            // done with the pop?
            
            // I don't need to share geometry. I use geometry only
            // once (though I may need to be careful on the rules with
            // final_geom), though that's not yet implemented.

            // Deriving automatically puts the Copy constraint on A,
            // and I am not sure I want to deal with that - but I have
            // to be able to copy Child regardless, thus Rule.

            // Function pointers support Copy, so Rule is fine.
            // Vectors by design do *not*.

            // Seems a little bizarre that none of this affects
            // recursive to_mesh... what am I doing differently?
            
            // See if it is time to backtrack:
            if n > max_depth || eval.children.is_empty() {
                // This has no parents:
                if n < 2 {
                    break;
                }

                // Backtrack:
                stack.pop();
                // TODO: Pop transform off of stack

                // If possible, step to the next sibling:
                let (parent, _) = &stack[n-2];
                if (idx + 1) < parent.children.len() {
                    let sib = parent.children[idx + 1];
                    match sib.rule {
                        Rule::Recurse(f) => {
                            let eval_sib = f(arg);
                            stack.push((eval_sib, idx + 1));
                            // TODO: Push transform onto stack
                            // TODO: Append geometry
                        },
                        Rule::EmptyRule => {
                            // Nowhere to recurse further
                        } 
                    }
                }
                continue;
            } else {
                // Otherwise, try to recurse to first child:
                let child = eval.children[0];
                match child.rule {
                    Rule::Recurse(f) => {
                        let eval_child = f(arg);
                        stack.push((eval_child, 0));
                        // TODO: Push transform onto stack
                        // TODO: Append geometry  
                    }
                    Rule::EmptyRule => {
                        // Do nothing.
                    }
                }
            }
        }

        // TODO: Return right number
        return (geom, 0); 

    }

    
}
