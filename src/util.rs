use std::ops::Range;
use crate::mesh::{Mesh, MeshFunc, VertexUnion};
use crate::xform::{Vertex};

/// This is like `vec!`, but it can handle elements that are given
/// with `@var element` rather than `element`, e.g. like
/// `vec_indexed![foo, bar, @a baz, quux]`. The variable (which must
/// already be declared and a `usize`) is then assigned the index of the
/// element it precedes (2).  This can be used any number of times with
/// different elements and indices.
///
/// It can  also be used like `vec_indexed![foo, bar, baz, @b,]` in which
/// case `b` is the index after 'baz' (3) rather than before. This still
/// requires a trailing comma.
#[macro_export]
macro_rules! vec_indexed {
    // Thank you to GhostOfSteveJobs and Rantanen in the Rust discord.
    ($( $(@ $Index:ident)? $($Value:expr)?,)*) => {{
        let mut v = Vec::new();
        $( $($Index = v.len();)? $(v.push($Value);)? )*
        v
    }};
}

pub trait VecExt<T> {
    fn append_indexed(&mut self, other: Vec<T>) -> (usize, usize);
}

impl<T> VecExt<T> for Vec<T> {
    // Like `append`, but returning `(a, b)` which give the range of
    // elements just inserted.
    fn append_indexed(&mut self, mut other: Vec<T>) -> (usize, usize) {
        let a = self.len();
        self.append(&mut other);
        let b = self.len();
        (a, b)
    }
}

/// Linearly subdivides a list of points that are to be treated as a
/// cycle.  This produces 'count' points for every element of 'p'
/// (thus, the returned length will be `count*p.len()`).
pub fn subdivide_cycle(p: &Vec<Vertex>, count: usize) -> Vec<Vertex> {
    let n = p.len();
    (0..n).map(|i| {
        // The inner loop is interpolating between v1 and v2:
        let v1 = &p[i];
        let v2 = &p[(i+1) % n];
        (0..count).map(move |j| {
            // This is just lerping in count+1 equally spaced steps:
            let f = (j as f32) / (count as f32);
            v1*(1.0-f) + v2*f
        })
    }).flatten().collect()
}
// TODO: This can be generalized to an iterator or to IntoIterator
// trait bound

pub fn parallel_zigzag_faces(r1: Range<usize>, r2: Range<usize>) -> Vec<usize> {
    let count = r1.end - r1.start;
        if (r2.end - r2.start) != count {
        panic!("Ranges must have the same size");
    }
    if (r2.end > r1.start && r2.end < r1.end) ||
        (r1.end > r2.start && r1.end < r2.end) {
        panic!("Ranges cannot overlap");
    }

    (0..count).map(|i0| {
        // i0 is an *offset* for the 'current' index.
        // i1 is for the 'next' index, wrapping back to 0.
        let i1 = (i0 + 1) % count;
        vec![
            // Mind winding order!
            r1.start + i1, r2.start + i0, r1.start + i0,
            r2.start + i1, r2.start + i0, r1.start + i1,
        ]
    }).flatten().collect()
}

pub fn parallel_zigzag(verts: Vec<VertexUnion>, main: Range<usize>, parent: Range<usize>) -> MeshFunc {
    MeshFunc {
        verts: verts,
        faces: parallel_zigzag_faces(main, parent),
    }
}

pub fn parallel_zigzag_mesh(verts: Vec<Vertex>, main: Range<usize>, parent: Range<usize>) -> Mesh {
    Mesh {
        verts: verts,
        faces: parallel_zigzag_faces(main, parent),
    }
}

pub fn parallel_zigzag2<T,U>(main: T, parent: U) -> Vec<usize>
where T: IntoIterator<Item=usize>, U: IntoIterator<Item=usize>
{
    let m: Vec<usize> = main.into_iter().collect();
    let p: Vec<usize> = parent.into_iter().collect();
    let l = m.len();
    if l != p.len() {
        panic!("Vectors must be the same size!")
    }

    (0..l).map(|i0| {
        // i0 is an *offset* for the 'current' index.
        // i1 is for the 'next' index, wrapping back to 0.
        let i1 = (i0 + 1) % l;
        vec![
            // Mind winding order!
            m[i1], p[i0], m[i0],
            p[i1], p[i0], m[i1],
        ]
    }).flatten().collect()
}

pub fn centroid(verts: &Vec<Vertex>) -> Vertex {
    let n = verts.len();
    let mut centroid = Vertex::new(0.0, 0.0, 0.0, 0.0);
    for v in verts {
        centroid += v;
    }
    centroid /= n as f32;
    centroid
}

pub fn connect_convex(range: Range<usize>, target: usize, as_parent: bool) -> Vec<usize> {

    let count = range.end - range.start;
    if as_parent {
        (0..count).map(|i0| {
            let i1 = (i0 + 1) % count;
            vec![range.start + i1, range.start + i0, target]
        }).flatten().collect()
    } else {
        (0..count).map(|i0| {
            let i1 = (i0 + 1) % count;
            vec![range.start + i0, range.start + i1, target]
        }).flatten().collect()
    }
}
