use crate::openmesh::{Tag, Vertex};
//use crate::rule::{Rule, Child};

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

pub fn parallel_zigzag_faces(count: usize) -> Vec<Tag> {
    (0..count).map(|i0| {
        let i1 = (i0+1) % count;
        vec![
            Tag::Body(i1),   Tag::Parent(i0), Tag::Body(i0),
            Tag::Parent(i1), Tag::Parent(i0), Tag::Body(i1),
        ]
    }).flatten().collect()
}

pub fn connect_convex(verts: &Vec<Vertex>, as_parent: bool) -> (Vertex, Vec<Tag>) {
    let n = verts.len();
    let mut centroid = Vertex::new(0.0, 0.0, 0.0, 0.0);
    for v in verts {
        centroid += v;
    }
    centroid /= n as f32;

    let faces: Vec<Tag> = {
        if as_parent {
            (0..n).map(|f1| {
                let f2 = (f1 + 1) % n;
                vec![Tag::Parent(f2), Tag::Parent(f1), Tag::Body(0)]
            }).flatten().collect()
        } else {
            (0..n).map(|f1| {
                let f2 = (f1 + 1) % n;
                // n is used for new center vertex
                vec![Tag::Body(f1), Tag::Body(f2), Tag::Body(n)]
            }).flatten().collect()
        }
    };

    (centroid, faces)
}
