use nalgebra::*;
use std::ops::Mul;

/// A type for mesh vertices. Initialize with [vertex][self::vertex].
pub type Vertex = Vector4<f32>;
/// A type for homogeneous transforms
pub type Mat4 = Matrix4<f32>;

#[derive(Clone, Copy)]
pub struct Transform {
    pub mtx: Mat4,
}

/// Initializes a vertex:
pub fn vertex(x: f32, y: f32, z: f32) -> Vertex {
    Vertex::new(x, y, z, 1.0)
}

impl Transform {

    pub fn new() -> Transform {
        Transform {
            mtx: geometry::Transform3::identity().to_homogeneous(),
        }
    }

    pub fn rotate(&self, axis: &Unit<Vector3<f32>>, angle: f32) -> Transform {
        Transform {
            mtx: self.mtx * Matrix4::from_axis_angle(axis, angle),
        }
    }

    pub fn translate(&self, x: f32, y: f32, z: f32) -> Transform {
        let xf = geometry::Translation3::new(x, y, z);
        Transform {
            mtx: self.mtx * xf.to_homogeneous(),
        }
    }

    pub fn scale(&self, f: f32) -> Transform {
        let xf = Matrix4::new_scaling(f);
        Transform {
            mtx: self.mtx * xf,
        }
    }

    /// Transforms a vector of vertices:
    pub fn transform(&self, verts: &Vec<Vertex>) -> Vec<Vertex> {
        verts.iter().map(|v| self.mtx * v).collect()
    }
}

impl Mul for Transform {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Transform { mtx: self.mtx * rhs.mtx }
    }
}

/// Convenience function for identity transformation
pub fn id() -> Transform {
    Transform::new()
}

/*
impl<'a> Mul for &'a Transform {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        &Transform { mtx: self.mtx * rhs.mtx }
    }
}
*/
