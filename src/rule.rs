use crate::mesh::{MeshFunc, VertexUnion};
use crate::xform::{Transform};
//use crate::prim;
use std::borrow::Borrow;
use std::rc::Rc;

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
