use std::f32::consts::{FRAC_PI_2, FRAC_PI_3, PI};
use std::f32;

use nalgebra::*;
use rand::Rng;

use crate::util;
use crate::util::VecExt;
use crate::mesh::{Mesh};
use crate::xform::{Transform, Vertex, vertex, id};

pub struct Barbs {
    base_incr: Transform,
    barb_incr: Transform,
    sides: [Transform; 4],
    base: Vec<Vertex>,
    verts: Vec<Vertex>,
    faces: Vec<usize>,
}

impl Barbs {
    pub fn new() -> Barbs {
        // Incremental transform from each stage to the next:
        let base_incr = id().translate(0.0, 0.0, 1.0).
            rotate(&Vector3::z_axis(), 0.15).
            rotate(&Vector3::x_axis(), 0.1).
            scale(0.95);
        let barb_incr = id().translate(0.0, 0.0, 0.5).
            rotate(&Vector3::y_axis(), -0.2).
            scale(0.8);
        // 'Base' vertices, used throughout:
        let base = vec![
            vertex(-0.5, -0.5, 0.0),
            vertex(-0.5,  0.5, 0.0),
            vertex( 0.5,  0.5, 0.0),
            vertex( 0.5, -0.5, 0.0),
        ];
        // sides[i] gives transformation from a 'base' layer to the
        // i'th side (0 to 3):
        let mut sides: [Transform; 4] = [id(); 4];
        for i in 0..4 {
            sides[i] = id().
                rotate(&Vector3::z_axis(), -FRAC_PI_2 * (i as f32)).
                rotate(&Vector3::y_axis(), -FRAC_PI_2).
                translate(0.5, 0.0, 0.5);// * barb_incr;
        }
        Barbs {
            base_incr: base_incr,
            barb_incr: barb_incr,
            base: base,
            sides: sides,
            verts: vec![],
            faces: vec![],
        }
    }

    pub fn run(mut self, iters: usize) -> Mesh {
        // Make seed vertices, use them for 'bottom' face, and recurse:
        self.verts.append(&mut self.base.clone());
        self.faces.extend_from_slice(&[0, 1, 2,   0, 2, 3]);
        self.main(iters, id(), [0,1,2,3]);
        return Mesh {
            verts: self.verts,
            faces: self.faces,
        }
    }

    fn limit_check(&self, xform: &Transform, _: usize) -> bool {
        // Assume all scales are the same (for now)
        let (s, _, _) = xform.get_scale();
        return s < 0.005;
    }

    fn main(&mut self, iters: usize, xform: Transform, bound: [usize; 4]) {

        if self.limit_check(&xform, iters) {
            self.faces.extend_from_slice(&[bound[0], bound[2], bound[1],
                bound[0], bound[3], bound[2]]);
            return;
        }

        let xform2 = xform * self.base_incr;
        let g = xform2.transform(&self.base);
        let (a0, _) = self.verts.append_indexed(g);

        // TODO: Isn't there some cleaner way?
        self.main(iters - 1, xform2, [a0, a0+1, a0+2, a0+3]);
        self.barb(iters - 1, xform * self.sides[0], [bound[0], bound[1], a0+1, a0+0]);
        self.barb(iters - 1, xform * self.sides[1], [bound[1], bound[2], a0+2, a0+1]);
        self.barb(iters - 1, xform * self.sides[2], [bound[2], bound[3], a0+3, a0+2]);
        self.barb(iters - 1, xform * self.sides[3], [bound[3], bound[0], a0+0, a0+3]);
    }

    fn barb(&mut self, iters: usize, xform: Transform, bound: [usize; 4]) {

        if self.limit_check(&xform, iters) {
            self.faces.extend_from_slice(&[bound[0], bound[2], bound[1],
                                           bound[0], bound[3], bound[2]]);
            return;
        }

        let xform2 = xform * self.barb_incr;
        let g = xform2.transform(&self.base);
        let (a0, a1) = self.verts.append_indexed(g);
        self.faces.append(&mut util::parallel_zigzag2(bound.to_vec(), a0..a1));

        self.barb(iters - 1, xform2, [a0, a0+1, a0+2, a0+3]);
    }
}

pub struct TreeThing {
    incr: Transform,
    splits: [Transform; 4],
    base: Vec<Vertex>,
    trans: Vec<Vertex>,
    verts: Vec<Vertex>,
    faces: Vec<usize>,
    depth: usize,
}

impl TreeThing {
    pub fn new(f: f32, depth: usize) -> TreeThing {
        // Incremental transform for each layer:
        let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
        let incr: Transform = Transform::new().
            translate(0.0, 0.0, 0.9 * f).
            rotate(&v, 0.4 * f).
            scale(1.0 - (1.0 - 0.95) * f);
        // 'Base' vertices, used throughout:
        let base = vec![
            vertex(-0.5, -0.5, 0.0),
            vertex(-0.5,  0.5, 0.0),
            vertex( 0.5,  0.5, 0.0),
            vertex( 0.5, -0.5, 0.0),
        ];
        // 'Transition' vertices:
        let trans = vec![
            // Top edge midpoints:
            vertex(-0.5,  0.0, 0.0),  // 0 - connects b0-b1
            vertex( 0.0,  0.5, 0.0),  // 1 - connects b1-b2
            vertex( 0.5,  0.0, 0.0),  // 2 - connects b2-b3
            vertex( 0.0, -0.5, 0.0),  // 3 - connects b3-b0
            // Top middle:
            vertex( 0.0,  0.0, 0.0),  // 4 - midpoint/centroid of all
        ];
        // splits[i] gives transformation from a 'base' layer to the
        // i'th split (0 to 3):
        let mut splits: [Transform; 4] = [id(); 4];
        for i in 0..4 {
            let r = FRAC_PI_2 * (i as f32);
            splits[i] = id().
                rotate(&nalgebra::Vector3::z_axis(), r).
                translate(0.25, 0.25, 0.0).
                scale(0.5);
        }
        TreeThing {
            incr: incr,
            splits: splits,
            base: base,
            trans: trans,
            verts: vec![],
            faces: vec![],
            depth: depth,
        }
    }

    pub fn run(mut self) -> Mesh {
        // Make seed vertices, use them for 'bottom' face, and recurse:
        self.verts.append(&mut self.base.clone());
        self.faces.extend_from_slice(&[0, 1, 2,   0, 2, 3]);
        self.child(id(), self.depth, [0,1,2,3]);
        return Mesh {
            verts: self.verts,
            faces: self.faces,
        }
    }

    pub fn run2(mut self) -> Mesh {
        // Make seed vertices, use them for 'bottom' face, and recurse:
        self.verts.append(&mut self.base.clone());
        self.faces.extend_from_slice(&[0, 1, 2,   0, 2, 3]);
        self.trunk(id(), [0,1,2,3]);
        return Mesh {
            verts: self.verts,
            faces: self.faces,
        }
    }

    pub fn trunk(&mut self, xform: Transform, b: [usize; 4]) {

        if self.limit_check(&xform) {
            self.faces.extend_from_slice(&[b[0], b[2], b[1], b[0], b[3], b[2]]);
            return;
        }

        let incr = id().translate(0.0, 0.0, 1.0).
            rotate(&Vector3::z_axis(), 0.15).
            rotate(&Vector3::x_axis(), 0.1).
            scale(0.95);
        let mut sides: [Transform; 4] = [id(); 4];
        for i in 0..4 {
            sides[i] = id().
                rotate(&Vector3::z_axis(), -FRAC_PI_2 * (i as f32)).
                rotate(&Vector3::y_axis(), -FRAC_PI_2).
                translate(0.5, 0.0, 0.5);// * barb_incr;
        }

        let xform2 = xform * incr;
        let g = xform2.transform(&self.base);
        let (a0, _) = self.verts.append_indexed(g);

        self.trunk(xform2, [a0, a0+1, a0+2, a0+3]);
        self.child(xform * sides[0], self.depth, [b[0], b[1], a0 + 1, a0 + 0]);
        self.child(xform * sides[1], self.depth, [b[1], b[2], a0 + 2, a0 + 1]);
        self.child(xform * sides[2], self.depth, [b[2], b[3], a0 + 3, a0 + 2]);
        self.child(xform * sides[3], self.depth, [b[3], b[0], a0 + 0, a0 + 3]);
    }

    fn limit_check(&self, xform: &Transform) -> bool {
        // Assume all scales are the same (for now)
        let (s, _, _) = xform.get_scale();
        return s < 0.01;
    }

    fn child(&mut self, xform: Transform, depth: usize, b: [usize; 4]) {

        if self.limit_check(&xform) {
            self.faces.extend_from_slice(&[b[0], b[2], b[1], b[0], b[3], b[2]]);
            return;
        }

        if depth > 0 {
            // Just recurse on the current path:
            let xform2 = xform * self.incr;
            let (n0, n1) = self.verts.append_indexed(xform2.transform(&self.base));
            self.faces.append(&mut util::parallel_zigzag2(n0..n1, b.to_vec()));

            self.child(xform2, depth - 1, [n0, n0 + 1, n0 + 2, n0 + 3]);

        } else {
            // 'Transition' stage (splits to 4 parts):
            let (n, _) = self.verts.append_indexed(xform.transform(&self.base));
            let (m01, _) = self.verts.append_indexed(xform.transform(&self.trans));
            let (m12, m23, m30, c) = (m01 + 1, m01 + 2, m01 + 3, m01 + 4);
            self.faces.extend_from_slice(&[
                // two faces straddling edge from vertex 0:
                b[0], n+0, m01,
                b[0], m30, n+0,
                // two faces straddling edge from vertex 1:
                b[1], n+1, m12,
                b[1], m01, n+1,
                // two faces straddling edge from vertex 2:
                b[2], n+2, m23,
                b[2], m12, n+2,
                // two faces straddling edge from vertex 3:
                b[3], n+3, m30,
                b[3], m23, n+3,
                // four faces from edge (0,1), (1,2), (2,3), (3,0):
                b[0], m01, b[1],
                b[1], m12, b[2],
                b[2], m23, b[3],
                b[3], m30, b[0],
            ]);

            self.child(xform * self.splits[0], self.depth,[c, m12, n+2, m23]);
            self.child(xform * self.splits[1], self.depth,[c, m01, n+1, m12]);
            self.child(xform * self.splits[2], self.depth,[c, m30, n+0, m01]);
            self.child(xform * self.splits[3], self.depth,[c, m23, n+3, m30]);
        }
    }
}

pub struct Sierpinski {
    splits: [Transform; 3],
    base: Vec<Vertex>,
    verts: Vec<Vertex>,
    faces: Vec<usize>,
}

impl Sierpinski {
    // s = Scale factor (0.5 for normal Sierpinski)
    // dz = Initial height step:
    // dr = 'Extra' z rotation (0.0 for normal Sierpinski)
    pub fn new(s: f32, dz: f32, dr: f32) -> Sierpinski {
        let mut splits: [Transform; 3] = [id(); 3];
        let rt3 = (3.0).sqrt();
        for i in 0..3 {
            let angle = 2.0 * FRAC_PI_3 * (i as f32) + dr;
            splits[i] = id().
                rotate(&Vector3::z_axis(), angle).
                translate(rt3/12.0, 0.0, 0.0).
                scale(s).
                translate(0.0, 0.0, dz);
        }

        let v0 = vertex(rt3/3.0, 0.0, 0.0);
        let v1 = vertex(-rt3/6.0, 1.0/2.0, 0.0);
        let v2 = vertex(-rt3/6.0, -1.0/2.0, 0.0);
        let v0b = v0 + vertex(0.0, 0.0, dz);
        let v1b = v1 + vertex(0.0, 0.0, dz);
        let v2b = v2 + vertex(0.0, 0.0, dz);

        let base = vec![
            v0, v1, v2,
            v0b, v1b, v2b,
            (v0b+v1b)/2.0,
            (v1b+v2b)/2.0,
            (v2b+v0b)/2.0,
        ];
        // There might be some redundancy when these are used

        Sierpinski {
            splits: splits,
            base: base,
            verts: vec![],
            faces: vec![],
        }
    }


    pub fn run(mut self) -> Mesh {
        self.verts.append(&mut self.base.clone());
        self.faces.extend_from_slice(&[2, 1, 0]);

        self.tri_split(id(), [0, 1, 2]);
        return Mesh {
            verts: self.verts,
            faces: self.faces,
        }
    }

    pub fn tri_split(&mut self, xform: Transform, b: [usize; 3]) {

        if self.limit_check(&xform) {
            self.faces.extend_from_slice(&b);
            return;
        }

        let g = xform.transform(&self.base);
        let (a, _) = self.verts.append_indexed(g);
        let t = a+3;
        let (tm01, tm12, tm20) = (t+3, t+4, t+5);

        self.faces.extend_from_slice(&[
            // Outer:
            tm01, b[1],  t+1,
            tm01, t+0,   b[0],
            tm01, b[0],  b[1],
            tm12, b[2],  t+2,
            tm12, t+1,   b[1],
            tm12, b[1],  b[2],
            tm20, b[0],  t+0,
            tm20, t+2,   b[2],
            tm20, b[2],  b[0],
            // Inner:
            tm01, tm12, tm20,
        ]);

        self.tri_split(xform * self.splits[0], [t+0, tm01, tm20]);
        self.tri_split(xform * self.splits[1], [t+1, tm12, tm01]);
        self.tri_split(xform * self.splits[2], [t+2, tm20, tm12]);
    }

    fn limit_check(&self, xform: &Transform) -> bool {
        // Assume all scales are the same (for now)
        let (s, _, _) = xform.get_scale();
        return s < 0.01;
    }
}

pub struct NestedSpiral {
    base: Vec<Vertex>,
    verts: Vec<Vertex>,
    faces: Vec<usize>,
    seeds: Vec<Vec<Transform>>,
    incr: Vec<Transform>,
}

impl NestedSpiral {
    pub fn new() -> NestedSpiral {
        // 'Base' vertices, used throughout:
        let mut base = vec![
            (vertex(0.0, 0.0, 0.0)),
            (vertex(1.0, 0.0, 0.0)),
            (vertex(1.0, 0.0, 1.0)),
            (vertex(0.0, 0.0, 1.0)),
        ];
        base = id().scale(0.3).translate(-0.5, 0.0, -0.5).transform(&base);

        // TODO: Subdivide vertices down (Python version goes from 4
        // verts up to 16 by subdividing twice)

        let rot_y = |ang| id().rotate(&Vector3::y_axis(), ang);
        let rot_y_and_trans = |ang, dx| rot_y(ang).translate(dx, 0.0, 0.0);


        let get_xform = |i, j, k| {
            let i0 = PI*2.0*(i as f32) / 4.0;
            let j0 = PI*2.0*(j as f32) / 4.0;
            let k0 = PI*2.0*(k as f32) / 4.0;
            // Pairs of (seed transform, incremental transform):
            vec![
                id(),
                rot_y_and_trans(i0, 3.0),
                rot_y_and_trans(j0, 1.0),
                rot_y_and_trans(k0, 0.5),
            ]
            // TODO: Does it make sense that this should be the reverse of the Python ones?
        };

        let mut seeds = vec![];
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    seeds.push(get_xform(i, j, k));
                }
            }
        }

        let incr = vec![
            id().translate(0.0, 0.1, 0.0),
            rot_y(-0.03),
            rot_y(0.07),
            rot_y(-0.2),
        ];
        // TODO: This is kludgy (having it separate 'incremental' and 'seed'
        // transformations)

        NestedSpiral {
            base: base,
            verts: vec![],
            faces: vec![],
            incr:  incr,
            seeds: seeds,
        }
    }

    pub fn iter(&mut self, count: usize, b: [usize; 4], xforms: &Vec<Transform>) {

        if count <= 0 {
            self.faces.extend_from_slice(&[b[0], b[2], b[1], b[0], b[3], b[2]]);
            return;
        }

        // Compute global transform with product of all 'xforms':
        let global = xforms.iter().fold(id(), |acc, m| acc * (*m));

        // Generate geometry:
        let g = global.transform(&self.base);
        let (n0, n1) = self.verts.append_indexed(g);
        self.faces.append(&mut util::parallel_zigzag2(n0..n1, b.to_vec()));

        // Increment the individual transformations:
        let xforms_next: Vec<Transform> = xforms.iter()
            .zip(self.incr.iter())
            .map(|(m,incr)| ((*incr) * (*m)))
            .collect();

        self.iter(count - 1, [n0, n0 + 1, n0 + 2, n0 + 3], &xforms_next);
    }

    pub fn run(mut self) -> Mesh {

        let seeds = self.seeds.clone();
        let mut rng = rand::thread_rng();

        for seed in seeds {
            let global = seed.iter().fold(id(), |acc, m| acc * (*m));
            let g = global.transform(&self.base);
            let (n0, _) = self.verts.append_indexed(g);
            self.faces.extend_from_slice(&[n0, n0+1, n0+2,   n0, n0+2, n0+3]);

            let count: usize = rng.gen_range(100, 500);
            // Old effect: just set count to 500

            self.iter(count, [n0, n0 + 1, n0 + 2, n0 + 3], &seed);
        }

        return Mesh {
            verts: self.verts,
            faces: self.faces,
        };
    }
}
