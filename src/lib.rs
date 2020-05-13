pub mod mesh;
pub mod rule;
pub mod prim;
#[macro_use]
pub mod util;
pub mod xform;
pub mod examples;

//pub use crate::examples;
//pub use crate::openmesh::test_thing;

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;
    use std::time::Instant;
    use rule::Rule;
    use nalgebra::*;

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
        mesh.to_mesh().write_stl_file(&fname).unwrap();
    }

    #[test]
    fn xform_order() {
        let geom = prim::cube();

        let y = &Vector3::y_axis();

        let dx = 4.0;
        let r = -0.5;
        
        let trans = xform::Transform::new().translate(dx, 0.0, 0.0);
        let rot = xform::Transform::new().rotate(y, r);

        let xf1 = trans.rotate(y, r);
        let xf2 = rot.translate(dx, 0.0, 0.0);

        // Rotate entire space, *then* translate in that rotated plane:
        geom.transform(&trans).transform(&rot).write_stl_file("xform_apply_trans_rot.stl").unwrap();
        geom.transform(&(rot * trans)).write_stl_file("xform_mul_rot_trans.stl").unwrap();
        geom.transform(&xf2).write_stl_file("xform_rot_trans.stl").unwrap();
        // Translate cube, *then* rotate it:
        geom.transform(&rot).transform(&trans).write_stl_file("xform_apply_rot_trans.stl").unwrap();
        geom.transform(&(trans * rot)).write_stl_file("xform_mul_trans_rot.stl").unwrap();
        geom.transform(&xf1).write_stl_file("xform_trans_rot.stl").unwrap();
    }

    /*
    // TODO: These tests don't test any conditions, so this is useful
    // short-hand to run, but not very meaningful as a test.
    #[test]
    fn cube_thing() {
        run_test(examples::cube_thing(), 3, "cube_thing3", false);
    }
    */

    #[test]
    fn barbs() { run_test(examples::barbs(), 50, "barbs", false); }

    /*
    #[test]
    fn twist() {
        run_test(examples::twist(1.0, 2), 200, "screw", false);
    }

    #[test]
    fn twisty_torus() {
        run_test(examples::twisty_torus(), 3000, "twisty_torus", false);
    }

    #[test]
    fn twisty_torus_hardcode() {
        run_test(examples::twisty_torus_hardcode(), 1000, "twisty_torus_hardcode", false);
    }
    
    #[test]
    #[ignore]
    fn twisty_torus_full() {
        run_test(examples::twisty_torus(), 40000, "twisty_torus_full", false);
    }

    #[test]
    #[ignore]
    fn wind_chime_mistake_thing() {
        run_test(examples::wind_chime_mistake_thing(), 400, "wind_chime_mistake_thing", false);
    }

    #[test]
    fn nest_spiral_2() {
        run_test(examples::nest_spiral_2(), 200, "nest_spiral_2", false);
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

    #[test]
    fn ramhorn_branch_random() {
        run_test(examples::ramhorn_branch_random(24, 0.25), 32, "ram_horn_branch_random", false);
    }
     */
}
// need this for now:
// cargo test -- --nocapture
// or: cargo test cube_thing -- --nocapture
