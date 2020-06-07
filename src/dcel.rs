use std::fmt;

use crate::mesh::{Mesh};
use crate::xform::{Vertex};

/// Doubly-connected edge list mesh (or a half-edge mesh),
/// parametrized over some vertex type.
#[derive(Clone, Debug)]
pub struct DCELMesh<V: Copy> {
    pub verts: Vec<DCELVertex<V>>,
    pub faces: Vec<DCELFace>,
    pub halfedges: Vec<DCELHalfEdge>,

    pub num_verts: usize,
    pub num_faces: usize,
    pub num_halfedges: usize,
}

/// A vertex of a mesh, combined with an arbitrary half-edge that has
/// this vertex as its origin.  This is always relative to some parent
/// Mesh<V>.
#[derive(Clone, Debug)]
pub struct DCELVertex<V> {
    /// The vertex itself.
    pub v: V,
    /// A half-edge (given as an index into 'halfedges');
    /// arbitrary, but `mesh.halfedges[halfedge] = v` must be true
    pub halfedge: usize,
}

/// A face, given as a half-edge that lies on its boundary (and must
/// traverse it counter-clockwise). This is always relative to some
/// parent Mesh<V>, as in Vertex.
#[derive(Clone, Debug)]
pub struct DCELFace {
    /// A boundary half-edge of this face (given as an index into
    /// 'halfedges').
    pub halfedge: usize,
}

/// A half-edge, given in terms of its origin vertex, the face that the
/// half-edge lies on the boundary of, its optional "twin" half-edge that
/// lies on an adjacent face, and previous and next half-edges (to
/// traverse the boundaries of the face). This is always relative to
/// some parent Mesh<V>, as in Vertex and Face.
#[derive(Clone, Debug)]
pub struct DCELHalfEdge {
    /// Origin vertex (given as an index into 'verts')
    pub vert: usize,
    /// Face this half-edge lies on the boundary of (given as an index
    /// into 'faces')
    pub face: usize,
    /// If false, ignore twin_halfedge.  (If this is true, then it must
    /// also be true for the twin.)
    pub has_twin: bool,
    /// The twin half-edge (given as an index into 'halfedges').
    /// The twin of the twin must point back to this HalfEdge.
    pub twin_halfedge: usize,
    /// The next half-edge on the boundary (given as an index into
    /// 'halfedges').  'prev_halfedge' of this half-edge must point
    /// back to this same HalfEdge.  Repeatedly following 'next_halfedge'
    /// must also lead back to this same HalfEdge.
    pub next_halfedge: usize,
    /// The previous half-edge on the boundary (given as an index into
    /// 'halfedges'). 'next_halfedge' of this half-edge must point
    /// back to this HalfEdge.  Repeatedly following 'prev_halfedge'
    /// must also lead back to this same HalfEdge.
    pub prev_halfedge: usize,
}

#[derive(Debug)]
pub enum VertSpec<V> {
    New(V),
    Idx(usize),
}

impl<V: Copy + std::fmt::Debug> fmt::Display for DCELMesh<V> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        let v_strs: Vec<String> = self.verts.iter().enumerate().map(|(i,v)| {
            format!("V{}=e{} {:?}", i, v.halfedge, v.v)
        }).collect();
        let v_str = v_strs.join(",");

        let f_strs: Vec<String> = self.faces.iter().enumerate().map(|(i,f)| {
            format!("F{}=e{}", i, f.halfedge)
        }).collect();
        let f_str = f_strs.join(", ");

        let e_strs: Vec<String> = self.halfedges.iter().enumerate().map(|(i,h)| {
            let twin = if h.has_twin {
                format!(" tw{}", h.twin_halfedge)
            } else {
                String::from("")
            };
            format!("E{}=v{} f{}{} n{} p{}", i, h.vert, h.face, twin, h.next_halfedge, h.prev_halfedge)
        }).collect();
        let e_str = e_strs.join(", ");

        write!(f, "DCELMesh({} verts, {}; {} faces, {}; {} halfedges, {})",
               self.num_verts, v_str,
               self.num_faces, f_str,
               self.num_halfedges, e_str)
    }
}

impl<V: Copy + std::fmt::Debug> DCELMesh<V> {

    pub fn new() -> DCELMesh<V> {
        DCELMesh {
            verts: vec![],
            faces: vec![],
            halfedges: vec![],
            num_verts: 0,
            num_faces: 0,
            num_halfedges: 0,
        }
    }

    pub fn print(&self) {

        println!("DCELMesh has {} verts, {} faces, {} halfedges:",
               self.num_verts,
               self.num_faces,
               self.num_halfedges);

        for (i,v) in self.verts.iter().enumerate() {
            println!("Vert {}: halfedge {}, {:?}", i, v.halfedge, v.v);
        }

        for (i,f) in self.faces.iter().enumerate() {
            println!("Face {}: halfedge {} (halfedges {:?}, verts {:?})",
                     i, f.halfedge, self.face_to_halfedges(i), self.face_to_verts(i));
        }

        for (i,h) in self.halfedges.iter().enumerate() {
            let twin = if h.has_twin {
                format!(", twin half-edge {}", h.twin_halfedge)
            } else {
                String::from("")
            };
            let v1 = self.verts[h.vert].v;
            let v2 = self.verts[self.halfedges[h.next_halfedge].vert].v;
            println!("Halfedge {}: vert {} (to {}), face {}, prev: {}, next: {}{}, ({:?} to {:?})",
                     i, h.vert, self.halfedges[h.next_halfedge].vert, h.face,
                     h.prev_halfedge, h.next_halfedge, twin,
                     v1, v2);
        }
    }

    /// Runs various checks on the mesh. This will return true if the mesh
    /// looks okay, and otherwise false. It will print messages as it
    /// runs.
    pub fn check(&self) -> bool {
        let mut pass = true;

        if self.num_halfedges != self.halfedges.len() {
            pass = false;
            println!("self.num_halfedges={} != self.halfedges.len()={}",
                     self.num_halfedges, self.halfedges.len());
        } else {
            println!("self.num_halfedges matches self.halfedges.len()");
        }

        if self.num_faces != self.faces.len() {
            pass = false;
            println!("self.faces={} != self.faces.len()={}",
                     self.num_faces, self.faces.len());
        } else {
            println!("self.num_faces matches self.faces.len()");
        }

        if self.num_verts != self.verts.len() {
            pass = false;
            println!("self.verts={} != self.verts.len()={}",
                     self.num_verts, self.verts.len());
        } else {
            println!("self.num_verts matches self.verts.len()");
        }

        for (i,v) in self.verts.iter().enumerate() {
            if v.halfedge >= self.halfedges.len() {
                println!("Vertex {}: halfedge index {} out of range",
                    i, v.halfedge);
                pass = false;
            }
            if self.halfedges[v.halfedge].vert != i {
                println!("Vertex {} names halfedge {}, which has a different origin vertex ({})",
                    i, v.halfedge, self.halfedges[v.halfedge].vert);
            }
        }

        for (i,edge) in self.halfedges.iter().enumerate() {
            if edge.vert >= self.verts.len() {
                println!("Halfedge {}: vertex index {} out of range", i, edge.vert);
                pass = false;
            }
            if edge.has_twin {
                let twin = &self.halfedges[edge.twin_halfedge];
                if !twin.has_twin {
                    println!("Halfedge {}: twin {} says it has no twin",
                        i, edge.twin_halfedge);
                    pass = false;
                } else if i != twin.twin_halfedge {
                    println!("Halfedge {} has twin {}, but reverse isn't true",
                             i, edge.twin_halfedge);
                    pass = false;
                } else if edge.vert != self.halfedges[twin.next_halfedge].vert {
                    println!("Halfedge {} starts at vertex {} but twin {} ends at vertex {}",
                             i, edge.vert, edge.twin_halfedge, self.halfedges[twin.next_halfedge].vert);
                    pass = false;
                }
            }
            let p = edge.prev_halfedge;
            if p >= self.halfedges.len() {
                println!("Halfedge {}: previous halfedge index {} out of range",
                    i, p);
                pass = false;
            }
            let n = edge.next_halfedge;
            if p >= self.halfedges.len() {
                println!("Halfedge {}: next halfedge index {} out of range",
                         i, n);
                pass = false;
            }
            let pn = self.halfedges[p].next_halfedge;
            if pn != i {
                println!("Halfedge {}: previous halfedge {} has next halfedge of {}, not {}",
                    i, p, pn, i);
                pass = false;
            }
            let np = self.halfedges[n].prev_halfedge;
            if np != i {
                println!("Halfedge {}: next halfedge {} has previous halfedge of {}, not {}",
                         i, n, np, i);
                pass = false;
            }
            // TODO: Check that following prev always leads back to start
            // likewise following next
        }

        for (i,face) in self.faces.iter().enumerate() {
            if face.halfedge >= self.halfedges.len() {
                println!("Face {}: halfedge index {} out of range",
                         i, face.halfedge);
                pass = false;
            }
            let face2 = self.halfedges[face.halfedge].face;
            if i != face2 {
                println!("Face {} gives boundary halfedge {}, which gives different face ({})",
                    i, face.halfedge, face2);
                pass = false;
            }

            // TODO: Check that face never visits same vertex twice?
            // This might belong in halfedge checking
        }

        if pass {
            println!("Mesh OK")
        } else {
            println!("Mesh has errors!")
        }

        pass
    }

    pub fn face_to_halfedges(&self, face_idx: usize) -> Vec<usize> {
        let mut edges: Vec<usize> = vec![];
        let start_idx = self.faces[face_idx].halfedge;
        edges.push(start_idx);

        let mut idx = self.halfedges[start_idx].next_halfedge;
        while start_idx != idx {
            edges.push(idx);
            idx = self.halfedges[idx].next_halfedge;
        }
        return edges;
    }

    pub fn face_to_verts(&self, face_idx: usize) -> Vec<usize> {
        self.face_to_halfedges(face_idx).iter().map(|e| {
            self.halfedges[*e].vert
        }).collect()
    }

    /// Adds a face that shares no edges with anything else  in the mesh.
    /// Returns (face index, half-edge indices); half-edge indices are
    /// given in the order of the vertices (i.e. the first half-edge's
    /// origin is verts[0], second is verts[1], third is verts[2]).
    pub fn add_face(&mut self, verts: [VertSpec<V>; 3]) -> (usize, [usize; 3]) {
        // *New* vertices will be at index v_n onward.
        let v_n = self.num_verts;
        // The face will be at index f_n:
        let f_n = self.num_faces;
        // The half-edges will be at indices e_n, e_n+1, e_n+2:
        let e_n = self.num_halfedges;

        // Half-edges and vertices can be inserted both at once:
        let mut new_verts: usize = 0;
        for i in 0..3 {
            let n = (i + 1) % 3;
            let p = (i + 2) % 3;
            // Either insert a new vertex, or use an existing one.
            // In both cases, 'v' is its index.
            let v = match verts[i] {
                VertSpec::New(v) => {
                    self.verts.push(DCELVertex {
                        v: v,
                        halfedge: e_n + i,
                    });
                    let idx = v_n + new_verts;
                    new_verts += 1;
                    idx
                },
                VertSpec::Idx(v) => v,
            };
            // Note that its half-edge is e_n + i, which technically
            // doesn't exist yet, but is inserted below:
            self.halfedges.push(DCELHalfEdge {
                vert: v,
                face: f_n,
                has_twin: false,
                twin_halfedge: 0,
                next_halfedge: e_n + n,
                prev_halfedge: e_n + p,
            });
        }
        self.num_halfedges += 3;
        self.num_verts += new_verts;

        // Finally, add the face (any halfedge is fine):
        self.faces.push(DCELFace { halfedge: e_n });
        self.num_faces += 1;

        (f_n, [e_n, e_n+1, e_n+2])
    }

    /// Add a face that lies on an existing boundary - i.e. one half-edge
    /// has a twin half-edge already on the mesh.  As this gives two
    /// vertices, only one other vertex needs specified.
    /// Returns (face index, halfedge indices).  Halfedge indices begin
    /// at the twin half-edge to the one specified.
    pub fn add_face_twin1(&mut self, twin: usize, vert: V) -> (usize, [usize; 3]) {
        // 'vert' will be at index v_n:
        let v_n = self.num_verts;
        // The half-edges will be at indices e_n, e_n+1, e_n+2:
        let e_n = self.num_halfedges;

        self.verts.push(DCELVertex {
            v: vert,
            halfedge: e_n + 2,
        });
        self.num_verts += 1;
        self.add_face_twin1_ref(twin, v_n)
    }

    /// Like `add_face_twin1`, but for a vertex already present in the
    /// mesh rather than a new one.  All else is identical.
    pub fn add_face_twin1_ref(&mut self, twin: usize, vert_idx: usize) -> (usize, [usize; 3]) {
        // The face will be at index f_n:
        let f_n = self.num_faces;
        // The half-edges will be at indices e_n, e_n+1, e_n+2:
        let e_n = self.num_halfedges;

        // Note the reversal of direction
        let twin_halfedge = &self.halfedges[twin];
        let v1 = self.halfedges[twin_halfedge.next_halfedge].vert;
        let v2 = twin_halfedge.vert;
        // twin is: v2 -> v1

        // Insert *its* twin, v1 -> v2, first:
        self.halfedges.push(DCELHalfEdge {
            vert: v1,
            face: f_n,
            has_twin: true,
            twin_halfedge: twin,
            next_halfedge: e_n + 1,
            prev_halfedge: e_n + 2,
        });
        // DEBUG
        if self.halfedges[twin].has_twin {
            panic!(format!("Trying to add twin to {}, which already has twin ({})",
                           twin, self.halfedges[twin].twin_halfedge));
        }
        self.halfedges[twin].has_twin = true;
        self.halfedges[twin].twin_halfedge = e_n;
        self.halfedges.push(DCELHalfEdge {
            vert: v2,
            face: f_n,
            has_twin: false,
            twin_halfedge: 0,
            next_halfedge: e_n + 2,
            prev_halfedge: e_n,
        });
        self.halfedges.push(DCELHalfEdge {
            vert: vert_idx,
            face: f_n,
            has_twin: false,
            twin_halfedge: 0,
            next_halfedge: e_n,
            prev_halfedge: e_n + 1,
        });

        self.num_halfedges += 3;

        // Finally, add the face (any halfedge is fine):
        self.faces.push(DCELFace { halfedge: e_n });
        self.num_faces += 1;

        (f_n, [e_n, e_n+1, e_n+2])
    }

    /// Add a face that lies on two connected boundaries - i.e. two of its
    /// half-edges have twins already on the mesh.
    ///
    /// Twin half-edges should be given in counter-clockwise order; that
    /// is, for the resultant face, one half-edge's twin will be twin1, and
    /// the next half-edge's twin will be twin2.
    /// Also: halfedge `twin2_idx` must end at the vertex that starts
    /// `twin1_idx`.
    ///
    /// Returns (face index, halfedge indices).  Halfedge indices begin
    /// at the twin halfedge to twin1, then twin halfedge of
    /// twin2, then the 'new' halfedge (which starts where twin2 starts,
    /// and ends where twin1 ends).
    pub fn add_face_twin2(&mut self, twin1_idx: usize, twin2_idx: usize) -> (usize, [usize; 3]) {
        // The face will be at index f_n:
        let f_n = self.num_faces;
        // The half-edges will be at indices e_n, e_n+1, e_n+2:
        let e_n = self.num_halfedges;

        // The origin vertex is 'between' the two edges, but because their
        // order is reversed (as twins), this is twin1's origin:
        let twin1 = &self.halfedges[twin1_idx];
        let twin2 = &self.halfedges[twin2_idx];
        let v1 = twin1.vert;
        let v2 = self.halfedges[twin1.next_halfedge].vert;
        // Final vertex is back around to twin2's origin:
        let v3 = twin2.vert;

        if v1 != self.halfedges[twin2.next_halfedge].vert {
            panic!(format!("twin2 ({}) must end where twin1 ({}) begins, but does not (vertex {} vs. {})",
                twin2_idx, twin1_idx, self.halfedges[twin2.next_halfedge].vert, v1));
        }

        // twin1 is: v1 -> v2, twin2 is: v3 -> v1.
        // so the twin of twin1 must be: v2 -> v1
        self.halfedges.push(DCELHalfEdge {
            vert: v2,
            face: f_n,
            has_twin: true,
            twin_halfedge: twin1_idx,
            next_halfedge: e_n + 1,
            prev_halfedge: e_n + 2,
        }); // index e_n
        self.halfedges[twin1_idx].has_twin = true;
        self.halfedges[twin1_idx].twin_halfedge = e_n;
        // and the twin of twin2 must be: v1 -> v3
        self.halfedges.push(DCELHalfEdge {
            vert: v1,
            face: f_n,
            has_twin: true,
            twin_halfedge: twin2_idx,
            next_halfedge: e_n + 2,
            prev_halfedge: e_n,
        }); // index e_n + 1
        self.halfedges[twin2_idx].has_twin = true;
        self.halfedges[twin2_idx].twin_halfedge = e_n + 1;
        // and final edge must be v3 -> v2:
        self.halfedges.push(DCELHalfEdge {
            vert: v3,
            face: f_n,
            has_twin: false,
            twin_halfedge: 0,
            next_halfedge: e_n,
            prev_halfedge: e_n + 1,
        }); // index e_n + 2
        self.num_halfedges += 3;

        // Finally, add the face (any halfedge is fine):
        self.faces.push(DCELFace { halfedge: e_n });
        self.num_faces += 1;

        (f_n, [e_n, e_n+1, e_n+2])
    }

    pub fn convert_mesh<F>(&self, f: F) -> Mesh
        where F: Fn(V) -> Vertex,
    {
        let n = self.faces.len();
        let mut faces: Vec<usize> = vec![0; 3 * n];

        for i in 0..n {

            let e0 = self.faces[i].halfedge;
            let h0 = &self.halfedges[e0];
            faces[3*i + 0] = h0.vert;
            let e1 = h0.next_halfedge;
            let h1 = &self.halfedges[e1];
            faces[3*i + 1] = h1.vert;
            let e2 = h1.next_halfedge;
            let h2 = &self.halfedges[e2];
            faces[3*i + 2] = h2.vert;
            if h2.next_halfedge != e0 {
                panic!(format!("Face {}: half-edges {},{},{} return to {}, not {}",
                               i, e0, e1, e2, h2.next_halfedge, e0));
            }
        }

        Mesh {
            verts: self.verts.iter().map(|e| f(e.v)).collect(),
            faces: faces,
        }
    }
}
