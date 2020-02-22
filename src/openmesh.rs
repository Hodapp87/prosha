//pub mod openmesh;

use nalgebra::*;
use std::fs::OpenOptions;
use std::io;

/// A type for custom mesh vertices. Initialize with [vertex][self::vertex].
pub type Vertex = Vector4<f32>;
pub type Mat4 = Matrix4<f32>;

/// Initializes a vertex:
pub fn vertex(x: f32, y: f32, z: f32) -> Vertex {
    Vertex::new(x, y, z, 1.0)
}

#[derive(Clone, Debug)]
pub enum Tag {
    Body(usize),
    Parent(usize),
}

#[derive(Clone, Debug)]
pub struct OpenMesh {
    // Vertices (in homogeneous coordinates).
    pub verts: Vec<Vertex>,
    // Triangles, taken as every 3 values, treated each as indices
    // into 'verts':
    pub faces: Vec<Tag>,
}

impl OpenMesh {
    
    pub fn transform(&self, xfm: Mat4) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            // TODO: Is the above faster if I pack vectors into a
            // bigger matrix, and transform that?
            faces: self.faces.clone(), // TODO: Use Rc?
        }
    }

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

        let get_vert = |j| {
            match self.faces[j] {
                Tag::Body(n) => self.verts[n].xyz(),
                Tag::Parent(_) => panic!("Cannot write_stl() if mesh has parent references!"),
            }
        };
        // TODO: Handle this behavior
        
        // Turn every face into an stl_io::Triangle:
        for i in 0..num_faces {
            let v0 = get_vert(3*i + 0);
            let v1 = get_vert(3*i + 1);
            let v2 = get_vert(3*i + 2);
            
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

    pub fn connect(&self, others: &Vec<(OpenMesh, &Vec<usize>)>) -> OpenMesh {
        // TODO: Clean up Vec<usize> stuff

        // Copy body vertices & faces:
        let mut verts: Vec<Vertex> = self.verts.clone();
        let mut faces = self.faces.clone();

        for (other,mapping) in others {

            // body_offset corresponds to the position in 'verts' at
            // which we're appending everything in 'other.verts' -
            // thus, the offset we shift all indices in 'others' by.
            let body_offset = verts.len();
            
            // Copy all vertices from 'other':
            verts.append(&mut other.verts.clone());

            // Append its faces:
            faces.extend(other.faces.iter().map(|t| {
                match t {
                    // Apply aforementioned shift to its body vertices:
                    Tag::Body(n) => Tag::Body(n + body_offset),
                    // Since 'self' vertices are in the same order,
                    // parent vertex references retain same index:
                    Tag::Parent(n) => Tag::Body(mapping[*n]),
                }
            }));
        }

        OpenMesh {
            verts: verts,
            faces: faces,
        }
    }
}
