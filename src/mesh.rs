use std::fs::OpenOptions;
use std::io;

use crate::xform::{Vertex, Transform};

/// Basic face-vertex mesh.  `faces` contains indices of `verts` and is
/// taken in groups of 3 for each triangle.
#[derive(Clone, Debug)]
pub struct Mesh {
    pub verts: Vec<Vertex>,
    pub faces: Vec<usize>,
}

impl Mesh {
    /// Returns a new `Mesh` whose vertices have been transformed.
    pub fn transform(&self, xfm: &Transform) -> Mesh {
        Mesh {
            verts: xfm.transform(&self.verts),
            // TODO: Is the above faster if I pack vectors into a
            // bigger matrix, and transform that?
            faces: self.faces.clone(), // TODO: Use Rc?
        }
    }

    /// Write this mesh as an STL file.  This will fail if any element
    /// of `faces` is `Tag::Parent`.
    pub fn write_stl_file(&self, fname: &str) -> io::Result<()> {
        let mut file = OpenOptions::new().write(true).create(true).truncate(true).open(fname)?;
        self.write_stl(&mut file)
    }

    fn write_stl<W: std::io::Write>(&self, writer: &mut W) -> io::Result<()> {

        // Every group of 3 indices in self.faces is one triangle, so
        // pre-allocate in the format stl_io wants:
        let num_faces = self.faces.len() / 3;
        let mut triangles = vec![stl_io::Triangle {
            normal: [0.0; 3],
            vertices: [[0.0; 3]; 3],
        }; num_faces];

        // Turn every face into an stl_io::Triangle:
        for i in 0..num_faces {
            let v0 = self.verts[self.faces[3*i + 0]].xyz();
            let v1 = self.verts[self.faces[3*i + 1]].xyz();
            let v2 = self.verts[self.faces[3*i + 2]].xyz();

            let normal = (v1-v0).cross(&(v2-v0));

            triangles[i].normal.copy_from_slice(&normal.as_slice());
            triangles[i].vertices[0].copy_from_slice(v0.as_slice());
            triangles[i].vertices[1].copy_from_slice(v1.as_slice());
            triangles[i].vertices[2].copy_from_slice(v2.as_slice());
            // TODO: Is there a cleaner way to do the above?
        }

        // I could also solve this with something like
        // https://doc.rust-lang.org/std/primitive.slice.html#method.chunks_exact
        // however I don't know what performance difference may be.

        stl_io::write_stl(writer, triangles.iter())
    }
}
