use nalgebra::*;
use crate::openmesh::{OpenMesh, Tag, vertex};

// is there a better way to do this?
pub fn empty_mesh() -> OpenMesh {
    OpenMesh {
        verts: vec![],
        faces: vec![],
        exit_groups: vec![],
    }
}

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
        exit_groups: vec![],
    }.transform(geometry::Translation3::new(-0.5, -0.5, -0.5).to_homogeneous())
}
