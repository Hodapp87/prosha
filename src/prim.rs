use crate::mesh::Mesh;
use crate::xform::{vertex, Transform};

/// Returns an empty mesh (no vertices, no faces).
pub fn empty_mesh() -> Mesh {
    Mesh {
        verts: vec![],
        faces: vec![],
    }
}

/// Returns a cube of sidelength one centered at (0,0,0).
pub fn cube() -> Mesh {
    Mesh {
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
            0, 3, 1, 0, 2, 3, 1, 7, 5, 1, 3, 7, 5, 6, 4, 5, 7, 6, 4, 2, 0, 4, 6, 2, 2, 7, 3, 2, 6,
            7, 0, 1, 5, 0, 5, 4,
        ],
    }
    .transform(&Transform::new().translate(-0.5, -0.5, -0.5))
}
