use crate::openmesh::{OpenMesh, Mat4};
use crate::prim;

/// Definition of a rule.  In general, a `Rule`:
///
/// - produces geometry when it is evaluated
/// - tells what other rules to invoke, and what to do with their
/// geometry
pub enum Rule {
    /// Produce some geometry, and possibly recurse further.
    Recurse(fn () -> RuleStep),
    /// Produce nothing and recurse no further.
    EmptyRule,
}
// TODO: Rename rules?
// TODO: It may be possible to have just a 'static' rule that requires
// no function call.
// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?

/// `RuleStep` supplies the results of evaluating some `Rule` for one
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
pub struct RuleStep {
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
    pub rule: Rule,

    /// The transform to apply to all geometry produced by `rule`
    /// (including its own `geom` and `final_geom` if needed, as well
    /// as all sub-geometry produced recursively).
    pub xf: Mat4,

    /// The mapping to apply to turn a Tag::Parent vertex reference
    /// into a vertex index of the parent mesh.  That is, if `rule`
    /// produces an `OpenMesh` with a face of `Tag::Parent(n)`, this
    /// will correspond to index `vmap[n]` in the parent mesh.
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

    /// Convert this `Rule` to mesh data, recursively.  `iters_left`
    /// sets the maximum recursion depth.  This returns (geometry,
    /// number of rule evaluations).
    pub fn to_mesh(&self, iters_left: u32) -> (OpenMesh, u32) {

        let mut evals: u32 = 1;

        if iters_left <= 0 {
            match self {
                Rule::Recurse(f) => {
                    let rs: RuleStep = f();
                    return (rs.final_geom, 1);
                }
                Rule::EmptyRule => {
                    return (prim::empty_mesh(), evals);
                }
            }
        }

        match self {
            Rule::Recurse(f) => {
                let rs: RuleStep = f();
                // TODO: This logic is more or less right, but it
                // could perhaps use some un-tupling or something.

                let subgeom: Vec<(OpenMesh, &Vec<usize>)> = rs.children.iter().map(|sub| {
                    // Get sub-geometry (still un-transformed):
                    let (submesh, eval) = sub.rule.to_mesh(iters_left - 1);
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
}
