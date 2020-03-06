use crate::openmesh::{OpenMesh, Mat4};
use crate::prim;

/// Definition of a rule.  In general, a `Rule`:
///
/// - produces geometry when it is evaluated
/// - tells what other rules to invoke, and what to do with their
/// geometry
pub struct Rule<A> {
    pub eval: fn (&A) -> RuleEval<A>,
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
    pub fn to_mesh(&self, arg: &A, iters_left: u32) -> (OpenMesh, usize) {

        let mut evals = 1;

        let rs: RuleEval<A> = (self.eval)(arg);
        if iters_left <= 0 {
            return (rs.final_geom, 1);
        }

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
        return (rs.geom.connect(&subgeom).0, evals);
    }

    /// This should be identical to to_mesh, but implemented
    /// iteratively with an explicit stack rather than with recursive
    /// function calls.
    pub fn to_mesh_iter(&self, arg: &A, max_depth: usize) -> (OpenMesh, usize) {

        struct State<A> {
            // The set of rules we're currently handling:
            rules: Vec<Child<A>>,
            // The next element of 'children' to handle:
            next: usize,
            // The world transform of the *parent* of 'rules', that
            // is, not including any transform of any element of
            // 'rules'.
            xf: Mat4,
        }

        // 'stack' stores at its last element our "current" State in
        // terms of a current world transform and which Child should
        // be processed next.  Every element prior to this is previous
        // states which must be kept around for further backtracking
        // (usually because they involve multiple rules).
        let mut stack: Vec<State<A>> = vec![];
        let mut geom = prim::empty_mesh();

        // Set up the stack's initial state - evaluate our own rule
        let eval = (self.eval)(arg);
        stack.push(State {
            rules: eval.children,
            next: 0,
            xf: nalgebra::geometry::Transform3::identity().to_homogeneous(),
        });
        geom = eval.geom;

        // Number of times we've evaluated a Rule:
        let mut eval_count = 1;

        // Stack depth (update at every push & pop):
        let mut n = stack.len();

        while !stack.is_empty() {

            // TODO: This, more elegantly?
            if eval_count > max_depth {
                break;
            }
            
            println!("DEBUG: stack has len {}", n);
            let s = &mut stack[n-1];

            if s.next >= s.rules.len() {
                // If we've run out of child rules, have the *parent* node (if one) move on:
                if n >= 2 {
                    stack[n-2].next += 1;
                }
                // and backtrack:
                stack.pop();
                n -= 1;
                // (if there isn't one, it makes no difference,
                // because the loop will end)
                continue;
            }

            let child = &s.rules[s.next];
            // Evaluate the rule:
            let mut eval = (child.rule.eval)(arg);
            eval_count += 1;

            // Make an updated world transform:
            let xf = s.xf * child.xf; // TODO: Check order on this

            // This rule produced some geometry which we'll
            // combine with the 'global' geometry:
            let new_geom = eval.geom.transform(&xf);
            println!("DEBUG: Connecting {} faces, vmap={:?}, faces={:?}",
                     new_geom.verts.len(), child.vmap, new_geom.faces);
            let (g, offsets) = geom.connect(&vec![(new_geom, &child.vmap)]);
            geom = g;

            // 'new_geom' may itself be parent geometry for
            // something in 'eval.children' (via Tag::Parent),
            // and vmap is there to resolve those Tag::Parent
            // references to the right vertices in 'new_geom'.
            //
            // However, we connect() on the global geometry
            // which we merged 'new_geom' into, not 'new_geom'
            // directly.  To account for this, we must shift
            // vmap by the offset that 'geom.connect' gave us:
            for (offset, child) in offsets.iter().zip(eval.children.iter_mut()) {
                child.vmap = child.vmap.iter().map(|n| {
                    n + offset
                }).collect();
            }

            // TODO: Why does below work?
            if (s.next + 1) >= s.rules.len() {
                let m = stack.len();
                if m >= 2 {
                    stack[m-2].next += 1;
                }
                stack.pop();
                n -= 1;
            }
            // I guess we are "done" with the rule after we've
            // evaluated it, and it is then safe to increment
            // s.next.
            
            // Recurse further (i.e. put more onto stack):                    
            stack.push(State {
                rules: eval.children,
                next: 0,
                xf: xf,
            });
            n += 1;
        }
        // TODO: Recursion depth? What does that even mean here?
        // Maybe something more like 'branch depth'?

        // TODO: Handle final_geom

        return (geom, eval_count); 
    }
    
}
