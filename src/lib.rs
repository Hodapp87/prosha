pub mod examples;
pub mod openmesh;
pub mod rule;
pub mod prim;
pub mod util;
pub mod xform;

//pub use crate::examples;
//pub use crate::openmesh::test_thing;

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;
    use std::time::Instant;
    use rule::Rule;

    fn run_test<S>(rule: Rule<S>, iters: usize, name: &str, use_old: bool) {
        let r = Rc::new(rule);
        println!("---------------------------------------------------");
        println!("Running {} with {}...",
                 name, if use_old { "to_mesh" } else { "to_mesh_iter" });
        if false {
            let start = Instant::now();
            let n = 5;
            for _ in 0..n {
                Rule::to_mesh_iter(r.clone(), iters);
            }
            let elapsed = start.elapsed();
            println!("DEBUG: {} ms per run", elapsed.as_millis() / n);
        }
        let mesh_fn = if use_old { Rule::to_mesh } else { Rule::to_mesh_iter };
        let (mesh, nodes) = mesh_fn(r.clone(), iters);
        println!("Evaluated {} rules to {} verts", nodes, mesh.verts.len());
        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        mesh.write_stl_file(&fname).unwrap();
    }

    // TODO: These tests don't test any conditions, so this is useful
    // short-hand to run, but not very meaningful as a test.
    
    #[test]
    fn cube_thing() {
        run_test(examples::cube_thing(), 3, "cube_thing3", false);
    }

    #[test]
    fn twist() {
        run_test(examples::twist(1.0, 2), 200, "screw", false);
    }

    // This one is very time-consuming to run:
    #[test]
    #[ignore]
    fn twist_full() {
        let f = 40;
        run_test(examples::twist(f as f32, 128), 100*f, "screw_full", false);
    }
    
    #[test]
    fn ramhorn() {
        run_test(examples::ramhorn(), 100, "ram_horn3", false);
    }

    #[test]
    fn ramhorn_branch() {
        run_test(examples::ramhorn_branch(24, 0.25), 32, "ram_horn_branch", false);
    }
}
// need this for now:
// cargo test -- --nocapture
// or: cargo test cube_thing -- --nocapture
