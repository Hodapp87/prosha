use std::fs::OpenOptions;
use std::io;
use std::borrow::Borrow;

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

    fn to_meshfunc(&self) -> MeshFunc {
        MeshFunc {
            faces: self.faces.clone(),
            verts: self.verts.iter().map(|v| VertexUnion::Vertex(*v)).collect(),
        }
    }
}

#[derive(Clone, Debug)]
pub enum VertexUnion {
    /// A concrete vertex.
    Vertex(Vertex),
    /// A vertex argument - something like an argument to a function with
    /// the given positional index.
    ///
    /// The job of `MeshFunc.connect` is to bind these arguments to concrete
    /// vertices.
    Arg(usize),
}

/// A face-vertex mesh whose vertices can either be concrete values
/// (as in `Mesh`) or aliases (indices to other vertices of some
/// hypothetical mesh).  This can be turned to `mesh` only if those
/// aliases are resolved to concrete vertices with a call like
/// `connect()`.
#[derive(Clone, Debug)]
pub struct MeshFunc {
    //
    pub verts: Vec<VertexUnion>,
    /// Indices of triangles (taken as every 3 values).  Indices begin
    /// with `alias_verts`, and then continue into `verts` - that is,
    /// from `0..alias_verts.len()` they refer to `alias_verts`, and from
    /// `alias_verts.len()..(alias_verts.len()+verts.len())` they refer to
    /// `verts`.
    pub faces: Vec<usize>,
}

impl MeshFunc {

    pub fn to_mesh(&self) -> Mesh {
        Mesh {
            faces: self.faces.clone(),
            verts: self.verts.iter().map(|v| match *v {
                VertexUnion::Vertex(v) => v,
                VertexUnion::Arg(_) => panic!("Mesh still has vertex arguments!"),
            }).collect(),
        }
    }

    /// Returns a new `MeshFunc` whose concrete vertices have
    /// been transformed. Note that alias vertices are left untouched.
    pub fn transform(&self, xfm: &Transform) -> MeshFunc {
        let v = self.verts.iter().map(|v| {
            match v {
                VertexUnion::Vertex(v) => VertexUnion::Vertex(xfm.mtx * v),
                a@_ => a.clone(),
            }
        });

        MeshFunc {
            verts: v.collect(),
            faces: self.faces.clone(),
        }
    }

    /// Appends any number of meshes together.  Returns both a single
    /// mesh, and a vector which gives the offset by which each
    /// corresponding input mesh was shifted.  That is, for the i'th
    /// index in `meshes`, all of its triangle indices were shifted by
    /// the i'th offset in the resultant mesh.
    pub fn append<T, U>(meshes: T) -> (MeshFunc, Vec<usize>)
        where U: Borrow<MeshFunc>,
              T: IntoIterator<Item = U>
    {
        let mut offsets: Vec<usize> = vec![];
        let mut v: Vec<VertexUnion> = vec![];
        let mut f: Vec<usize> = vec![];
        for mesh_ in meshes {
            let mesh = mesh_.borrow();

            // Position in 'verts' at which we're appending
            // mesh.verts, which we need to know to shift indices:
            let offset = v.len();
            offsets.push(offset);

            // Copy all vertices:
            v.append(&mut mesh.verts.clone());
            // Append its faces, applying offset:
            f.extend(mesh.faces.iter().map(|n| n + offset));
        }

        (MeshFunc { verts: v, faces: f }, offsets)
    }

    /// Treat this mesh as a 'parent' mesh to connect with any number
    /// of 'child' meshes, all of them paired with their respective
    /// vertex argument values (i.e. `arg_vals` from `Child`).
    /// This returns a tuple of (new mesh, offsets), where 'offsets'
    /// gives the offset of where child meshes were shifted in the new
    /// mesh.
    ///
    /// That is, the vertices of 'children[i]' begin at vertex
    /// 'offset[i]' of the new mesh. This is needed in order to adjust
    /// references to vertices of a mesh in 'children' - such as
    /// 'arg_vals' of `rule::Child`.
    pub fn connect<T, U>(&self, children: T) -> (MeshFunc, Vec<usize>)
    where U: Borrow<MeshFunc>,
          T: IntoIterator<Item = (U, Vec<usize>)>
    //pub fn connect(&self, children: &Vec<(OpenMesh, Vec<usize>)>) -> (OpenMesh, Vec<usize>)
    {
        // TODO: Clean up this description a bit
        // TODO: Clean up Vec<usize> stuff

        // Copy body vertices & faces:
        let mut verts: Vec<VertexUnion> = self.verts.clone();
        let mut faces = self.faces.clone();

        let mut offsets: Vec<usize> = vec![];
        
        for (child_,mapping) in children {

            let child = child_.borrow();
            
            // offset corresponds to the position in 'verts' at
            // which we're appending everything in 'child.verts' -
            // thus, the offset we shift all indices in 'children' by.
            let offset = verts.len();

            // Concrete vertices in 'child.verts' need to be copied:
            verts.extend(child.verts.iter().filter_map(|v| {
                match v {
                    VertexUnion::Vertex(_) => Some(v.clone()),
                    VertexUnion::Arg(_) => None,
                }
            }));

            // All faces need copied, but if if the index was to
            // a concrete vertex, then it needs shifted by 'offset';
            // if an alias, it needs remapped.
            faces.extend(child.faces.iter().map(|n|
                match child.verts[*n] {
                    VertexUnion::Vertex(_) => n + offset,
                    VertexUnion::Arg(m) => mapping[m],
                }
            ));

            offsets.push(offset);
        }

        let m = MeshFunc {
            verts: verts,
            faces: faces,
        };
        (m, offsets)
    }
}
