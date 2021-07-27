pub mod mesh;
#[macro_use]
pub mod prim;
#[macro_use]
pub mod util;
pub mod examples;
pub mod xform;

//pub use crate::examples;
//pub use crate::openmesh::test_thing;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::*;
    use std::rc::Rc;
    use std::time::Instant;

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
        geom.transform(&trans)
            .transform(&rot)
            .write_stl_file("xform_apply_trans_rot.stl")
            .unwrap();
        geom.transform(&(rot * trans))
            .write_stl_file("xform_mul_rot_trans.stl")
            .unwrap();
        geom.transform(&xf2)
            .write_stl_file("xform_rot_trans.stl")
            .unwrap();
        // Translate cube, *then* rotate it:
        geom.transform(&rot)
            .transform(&trans)
            .write_stl_file("xform_apply_rot_trans.stl")
            .unwrap();
        geom.transform(&(trans * rot))
            .write_stl_file("xform_mul_trans_rot.stl")
            .unwrap();
        geom.transform(&xf1)
            .write_stl_file("xform_trans_rot.stl")
            .unwrap();
    }

    // TODO: These tests don't test any conditions, so this is useful
    // short-hand to run, but not very meaningful as a test.
    #[test]
    fn barbs() {
        let name = "barbs";
        println!("---------------------------------------------------");
        let b = examples::Barbs::new();
        let m = b.run(100);

        println!("Got {} verts...", m.verts.len());

        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        m.write_stl_file(&fname).unwrap();
    }

    #[test]
    fn tree_thing1() {
        let name = "tree_thing";
        println!("---------------------------------------------------");
        let b = examples::TreeThing::new(0.6, 10);
        let m = b.run();

        println!("Got {} verts...", m.verts.len());

        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        m.write_stl_file(&fname).unwrap();
    }

    #[test]
    fn tree_thing2() {
        let name = "tree_thing2";
        println!("---------------------------------------------------");
        let b = examples::TreeThing::new(0.6, 10);
        let m = b.run2();

        println!("Got {} verts...", m.verts.len());

        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        m.write_stl_file(&fname).unwrap();
    }

    #[test]
    fn sierpinski() {
        let name = "sierpinski";
        println!("---------------------------------------------------");
        let b = examples::Sierpinski::new(0.50, 0.10, 0.0);
        //let b = examples::Sierpinski::new(0.51, 0.10, 0.1);
        let m = b.run();

        println!("Got {} verts...", m.verts.len());

        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        m.write_stl_file(&fname).unwrap();
    }

    #[test]
    fn nested_spiral() {
        let name = "nested_spiral";
        println!("---------------------------------------------------");
        let b = examples::NestedSpiral::new();
        //let b = examples::Sierpinski::new(0.51, 0.10, 0.1);
        let m = b.run();

        println!("Got {} verts...", m.verts.len());

        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        m.write_stl_file(&fname).unwrap();
    }
}
// need this for now:
// cargo test -- --nocapture
// or: cargo test cube_thing -- --nocapture
