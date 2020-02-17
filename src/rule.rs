use crate::openmesh::{OpenMesh, Mat4};
use crate::prim;

// TODO: Do I benefit with Rc<Rule> below so Rule can be shared?
pub enum Rule {
    // Produce geometry, and possibly recurse further:
    Recurse(fn () -> RuleStep),
    // Stop recursing here:
    EmptyRule,
}
// TODO: Rename rules?
// TODO: It may be possible to have just a 'static' rule that requires
// no function call.

pub struct RuleStep {
    // The geometry generated by this rule on its own (not by any of
    // the child rules).
    pub geom: OpenMesh,

    // The "final" geometry, used only if recursion must be stopped.
    // This should be in the same coordinate space as 'geom', and
    // properly close any exit groups that it may have (and have no
    // exit groups of its own).
    pub final_geom: OpenMesh,

    // Child rules, paired with the transform that will be applied to
    // all of their geometry
    pub children: Vec<(Rule, Mat4)>,
}

impl Rule {

    // TODO: Do I want to make 'geom' shared somehow, maybe with Rc? I
    // could end up having a lot of identical geometry that need not be
    // duplicated until it is transformed into the global space.
    //
    // This might produce bigger gains if I rewrite rule_to_mesh so that
    // rather than repeatedly transforming meshes, it stacks
    // transformations and then applies them all at once.

    pub fn to_mesh(&self, iters_left: u32) -> (OpenMesh, u32) {

        
        
        let mut nodes: u32 = 1;

        if iters_left <= 0 {
            match self {
                Rule::Recurse(f) => {
                    let rs: RuleStep = f();
                    return (rs.final_geom, 1);
                }
                Rule::EmptyRule => {
                    return (prim::empty_mesh(), nodes);
                }
            }
        }

        match self {
            Rule::Recurse(f) => {
                let rs: RuleStep = f();

                // Get sub-geometry (from child rules) and transform it:
                let subgeom: Vec<(OpenMesh, Mat4, u32)> = rs.children.iter().map(|(subrule, subxform)| {
                    let (m,n) = subrule.to_mesh(iters_left - 1);
                    (m, *subxform, n)
                }).collect();

                // Tally up node count:
                subgeom.iter().for_each(|(_,_,n)| nodes += n);

                let g: Vec<OpenMesh> = subgeom.iter().map(|(m,x,_)| m.transform(*x)).collect();

                // Connect geometry from this rule (not child rules):
                return (rs.geom.connect(&g), nodes);
            }
            Rule::EmptyRule => {
                return (prim::empty_mesh(), nodes);
            }
        }
    }
}
