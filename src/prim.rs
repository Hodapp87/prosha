use nalgebra::*;
use crate::openmesh::{OpenMesh, Tag, vertex};

/// Returns an empty mesh (no vertices, no faces).
pub fn empty_mesh() -> OpenMesh {
    OpenMesh {
        verts: vec![],
        faces: vec![],
    }
}

/// Returns a cube of sidelength one centered at (0,0,0).
pub fn cube() -> OpenMesh {
    OpenMesh {
        verts: vec![
            vertex(0.0, 0.0, 0.0),
            vertex(1.0, 0.0, 0.0),
            vertex(0.0, 1.0, 0.0),
            vertex(1.0, 1.0, 0.0),
            vertex(0.0, 0.0, 1.0),
            vertex(1.0, 0.0, 1.0),
            vertex(0.0, 1.0, 1.0),
            vertex(1.0, 1.0, 1.0),
        ],
        faces: vec![
            Tag::Body(0), Tag::Body(3), Tag::Body(1),
            Tag::Body(0), Tag::Body(2), Tag::Body(3),
            Tag::Body(1), Tag::Body(7), Tag::Body(5),
            Tag::Body(1), Tag::Body(3), Tag::Body(7),
            Tag::Body(5), Tag::Body(6), Tag::Body(4),
            Tag::Body(5), Tag::Body(7), Tag::Body(6),
            Tag::Body(4), Tag::Body(2), Tag::Body(0),
            Tag::Body(4), Tag::Body(6), Tag::Body(2),
            Tag::Body(2), Tag::Body(7), Tag::Body(3),
            Tag::Body(2), Tag::Body(6), Tag::Body(7),
            Tag::Body(0), Tag::Body(1), Tag::Body(5),
            Tag::Body(0), Tag::Body(5), Tag::Body(4),
        ],
    }.transform(&geometry::Translation3::new(-0.5, -0.5, -0.5).to_homogeneous())
}
