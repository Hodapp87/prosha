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
            panic!("Trying to add twin to {}, which already has twin ({})",
                   twin, self.halfedges[twin].twin_halfedge);
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
            panic!("twin2 ({}) must end where twin1 ({}) begins, but does not (vertex {} vs. {})",
                twin2_idx, twin1_idx, self.halfedges[twin2.next_halfedge].vert, v1);
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

    pub fn split_face(&mut self, face: usize, verts: Vec<V>) {
        // 'verts' maps 1:1 to vertices for 'face' (i.e. face_to_verts).

        let mut edge_idx = self.faces[face].halfedge;
        let n = verts.len();

        //
        let new_edges: Vec<(usize, usize)> = verts.iter().map(|v| {

            println!("DEBUG: halfedge 3: {:?}", self.halfedges[3]);
            println!("DEBUG: halfedge {}: {:?}", self.halfedges[3].next_halfedge, self.halfedges[self.halfedges[3].next_halfedge]);

            // As we iterate over every vertex, we walk the half-edges:
            let mut edge = self.halfedges[edge_idx].clone();
            if !edge.has_twin {
                panic!("Halfedge {} has no twin, and split_face (for now) requires twins", edge_idx);
            }
            // TODO: Remove the above limitation and just don't try to split
            // nonexistent twins. I think all logic works the same way.

            let next_idx = edge.next_halfedge;
            let twin_idx = edge.twin_halfedge;
            let mut twin = self.halfedges[twin_idx].clone();

            // This half-edge, and its twin, are both split (at this
            // vertex). Half-edges i & j will be the new ones created.
            let i_edge = self.num_halfedges;
            let j_edge = i_edge + 1;
            // Half of edge_idx is split into j_edge.
            // Half of twin_idx (its twin) is split into i_edge.

            println!("DEBUG: edge_idx={} next_idx={} twin_idx={} i_edge={} j_edge={}", edge_idx, next_idx, twin_idx, i_edge, j_edge);
            // This is where the vertex will be inserted:
            let v_idx = self.num_verts;
            println!("DEBUG: adding v_idx={}", v_idx);
            self.verts.push(DCELVertex {
                v: *v,
                halfedge: i_edge, // j_edge is also fine
            });
            self.num_verts += 1;

            edge.twin_halfedge = i_edge;
            let j_next = edge.next_halfedge;

            twin.twin_halfedge = j_edge;
            let i_next = twin.next_halfedge;
            println!("DEBUG: new twin of {} = {}; new twin of {} = {}", edge_idx, i_edge, twin_idx, j_edge);

            self.halfedges.push(DCELHalfEdge {
                vert: v_idx,
                face: 0, // This is set properly in the next loop
                has_twin: true,
                twin_halfedge: edge_idx,
                next_halfedge: i_next,
                prev_halfedge: twin_idx,
            }); // i_edge
            println!("DEBUG: edge {}: vert {} twin {} next {} prev {} (ends at vertex {})", i_edge, v_idx, edge_idx, i_next, twin_idx, self.halfedges[self.halfedges[twin_idx].next_halfedge].vert);

            self.halfedges.push(DCELHalfEdge {
                vert: v_idx,
                face: 0, // This is set properly in the next loop
                has_twin: true,
                twin_halfedge: twin_idx,
                next_halfedge: j_next,
                prev_halfedge: edge_idx,
            }); // j_edge
            println!("DEBUG: edge {}: vert {} twin {} next {} prev {} (ends at vertex {})", j_edge, v_idx, twin_idx, j_next, edge_idx, self.halfedges[self.halfedges[edge_idx].next_halfedge].vert);

            self.num_halfedges += 2;

            self.halfedges[edge_idx] = edge;
            self.halfedges[twin_idx] = twin;

            let r = (edge_idx, j_edge);

            edge_idx = next_idx;

            r
        }).collect();

        println!("DEBUG: {:?}", new_edges);

        // We then must connect some edges up 'across' vertices
        // in order to form the smaller face at each vertex.
        //
        // This is outside the loop because we need one iteration's
        // value (any iteration, doesn't matter) to reassign a face:
        let mut e_twin_idx = 0;
        let twin_edges: Vec<usize> = (0..n).map(|i0| {
            let i1 = (i0 + 1) % n;

            let (_, ep_idx) = new_edges[i0];
            let (en_idx, _) = new_edges[i1];

            // Halfedges will be inserted here:
            let e_cross_idx = self.num_halfedges;
            e_twin_idx = e_cross_idx + 1;
            // And the face here:
            let face_new = self.num_faces;

            println!("DEBUG: i0={} i1={} ep_idx={} en_idx={} e_cross_idx={} e_twin_idx={} face_new={}",
                i0, i1, ep_idx, en_idx, e_cross_idx, e_twin_idx, face_new);

            // So, the vertex for i0 had two halfedges (one pointing in,
            // one pointing out). We split both those half-edges earlier.
            // en_idx & ep_idx are the halves that are nearest the
            // vertex. The point of below is to form a smaller triangle
            // (which includes this vertex, and the point at which both
            // edges were split).  This requires a half-edge from ep_idx
            // (which *starts* where the other was split) to en_idx (which
            // *ends* where one edge was split).
            self.halfedges.push(DCELHalfEdge {
                vert: self.halfedges[self.halfedges[en_idx].twin_halfedge].vert,
                face: face_new,
                has_twin: true,
                twin_halfedge: e_twin_idx,
                next_halfedge: ep_idx,
                prev_halfedge: en_idx,
            }); // e_cross_idx
            println!("DEBUG: edge {}: vert {} twin {} next {} prev {} (ends at vertex {})", e_cross_idx, self.halfedges[en_idx].vert, e_twin_idx, ep_idx, en_idx, self.halfedges[self.halfedges[e_cross_idx].next_halfedge].vert);

            // It also requires a twin half-edge. These all form a single
            // central face with each edge sharing a boundary with the
            // 'cross' edge we just created.
            self.halfedges.push(DCELHalfEdge {
                vert: self.halfedges[ep_idx].vert,
                face: face, // Reuse index for the smaller *central* face
                has_twin: true,
                twin_halfedge: e_cross_idx,
                next_halfedge: 0, // TODO
                prev_halfedge: 0, // TODO
            }); // e_twin_idx
            println!("DEBUG: edge {}: vert {} twin {} next/prev incorrect", e_twin_idx, self.halfedges[ep_idx].vert, e_cross_idx);

            self.num_halfedges += 2;

            // en/ep also need directed to 'new' edges and each other
            self.halfedges[en_idx].prev_halfedge = ep_idx;
            self.halfedges[en_idx].next_halfedge = e_cross_idx;
            println!("DEBUG: edge {}: now next {} prev {} ends at vert {}", en_idx, e_cross_idx, ep_idx, self.halfedges[self.halfedges[en_idx].next_halfedge].vert);

            self.halfedges[ep_idx].next_halfedge = en_idx;
            self.halfedges[ep_idx].prev_halfedge = e_cross_idx;
            println!("DEBUG: edge {}: now next {} prev {} ends at vert {}", ep_idx, en_idx, e_cross_idx, self.halfedges[self.halfedges[ep_idx].next_halfedge].vert);

            self.halfedges[ep_idx].face = face_new;

            self.faces.push(DCELFace {
                halfedge: e_cross_idx, // en_idx or ep_idx is fine too
            });
            self.num_faces += 1;

            // We also need to split the opposite side to make the two
            // 'side' triangles (which means two new half-edges).
            // First, we have to find the halfedge that starts the
            // 'opposite' vertex (outer2):
            let base1 = self.halfedges[en_idx].twin_halfedge;
            let base2 = self.halfedges[base1].prev_halfedge;
            let outer1 = self.halfedges[base1].next_halfedge;
            let outer2 = self.halfedges[outer1].next_halfedge;
            let v_opp = self.halfedges[outer2].vert;
            // One face will reuse the old index:
            let face1 = self.halfedges[outer2].face;
            // Another will be inserted at this index:
            let face2 = self.num_faces;
            // Half-edges will be inserted here:
            let edge_side1 = self.num_halfedges;
            let edge_side2 = edge_side1 + 1;
            self.halfedges.push(DCELHalfEdge {
                vert: v_opp,
                face: face1,
                has_twin: true,
                twin_halfedge: edge_side2,
                next_halfedge: base1,
                prev_halfedge: outer1,
            }); // edge_side1
            self.halfedges.push(DCELHalfEdge {
                vert: self.halfedges[base1].vert,
                face: face2,
                has_twin: true,
                twin_halfedge: edge_side1,
                next_halfedge: outer2,
                prev_halfedge: base2,
            }); // edge_side2
            self.num_halfedges += 2;

            self.faces.push(DCELFace {
                halfedge: outer2, // base2 or edge_side2 is fine too
            });
            self.num_faces += 1;
            self.faces[face1].halfedge = outer1; // base1 or edge_side1 is fine too

            self.halfedges[outer1].next_halfedge = edge_side1;
            self.halfedges[outer1].prev_halfedge = base1;
            self.halfedges[outer1].face = face1;
            self.halfedges[base1].face = face1;
            self.halfedges[base1].prev_halfedge = edge_side1;
            self.halfedges[base1].next_halfedge = outer1;

            self.halfedges[outer2].next_halfedge = base2;
            self.halfedges[outer2].prev_halfedge = edge_side2;
            self.halfedges[outer2].face = face2;
            self.halfedges[base2].face = face2;
            self.halfedges[base2].prev_halfedge = outer2;
            self.halfedges[base2].next_halfedge = edge_side2;

            println!("DEBUG: base1={} base2={} outer1={} outer2={} face1={} face2={} edge_side1={} edge_side2={}",
                     base1, base2, outer1, outer2, face1, face2, edge_side1, edge_side2);
            e_twin_idx
        }).collect();

        println!("DEBUG: cross_edges={:?}", twin_edges);

        for i0 in 0..n {
            let i1 = (i0 + 1) % n;

            self.halfedges[twin_edges[i0]].next_halfedge = twin_edges[i1];
            self.halfedges[twin_edges[i1]].prev_halfedge = twin_edges[i0];
        }

        // We split up original 'face' completely and created four new faces,
        // We need something at this index, and the other three already have
        // indices, so reuse it for the smaller central face:
        self.faces[face].halfedge = e_twin_idx;
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
                panic!("Face {}: half-edges {},{},{} return to {}, not {}",
                       i, e0, e1, e2, h2.next_halfedge, e0);
            }
        }

        Mesh {
            verts: self.verts.iter().map(|e| f(e.v)).collect(),
            faces: faces,
        }
    }
}
