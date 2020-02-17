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
    Exit(usize, usize),  // (group, vertex)
}

#[derive(Clone, Debug)]
pub struct OpenMesh {
    // Vertices (in homogeneous coordinates).
    pub verts: Vec<Vertex>,
    // Triangles, taken as every 3 values, treated each as indices
    // into 'verts':
    pub faces: Vec<Tag>,
    pub exit_groups: Vec<usize>,
}

impl OpenMesh {
    
    pub fn transform(&self, xfm: Mat4) -> OpenMesh {
        OpenMesh {
            verts: self.verts.iter().map(|v| xfm * v).collect(),
            // TODO: Is the above faster if I pack vectors into a
            // bigger matrix, and transform that?
            faces: self.faces.clone(), // TODO: Use Rc?
            exit_groups: self.exit_groups.clone(),
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
                Tag::Exit(_, _) => panic!("Cannot write_stl() if mesh has exit groups!"),
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

    pub fn connect(&self, others: &Vec<OpenMesh>) -> OpenMesh {

        // Copy body vertices & faces:
        let mut verts: Vec<Vertex> = self.verts.clone();
        let mut faces = self.faces.clone();

        let mut exit_groups: Vec<usize> = vec![];
        
        let mut body_offset = self.verts.len();
        let mut exit_offset = 0;
        let mut offsets: Vec<usize> = vec![0; others.len()];
        for (i,other) in others.iter().enumerate() {

            // Append body vertices & exit vertices directly:
            verts.append(&mut other.verts.clone());
            
            // Append faces, shifting each kind by respective offset:
            faces.extend(other.faces.iter().map(|t| {
                match t {
                    Tag::Body(n) => Tag::Body(n + body_offset),
                    Tag::Exit(g, n) => Tag::Exit(g + exit_groups.len(), n + exit_offset),
                }
            }));
            if i < self.exit_groups.len() {
                exit_offset += self.exit_groups[i];
            }
            exit_groups.append(&mut other.exit_groups.clone());

            offsets[i] = body_offset;
            // Increase offsets by the number of elements we appended:
            body_offset += other.verts.len();
        }

        // All of the Exit face indices from 'self' need to be
        // modified to refer to Body vertices of something in
        // 'others':
        for i in 0..faces.len() {
            match faces[i] {
                Tag::Exit(g, n) => {
                    faces[i] = Tag::Body(n + offsets[g]);
                },
                _ => { },
            };
        }

        OpenMesh {
            verts: verts,
            faces: faces,
            exit_groups: exit_groups,
        }
    }
}
