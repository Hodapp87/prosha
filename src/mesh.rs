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

    pub fn to_meshfunc(&self) -> MeshFunc {
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

pub fn vert_args<T: IntoIterator<Item = usize>>(v: T) -> Vec<VertexUnion> {
    v.into_iter().map(|i| VertexUnion::Arg(i)).collect()
}

/// A face-vertex mesh, some of whose vertices may be 'vertex arguments'
/// rather than concrete vertices. The job of `connect()` is to resolve
/// these.
#[derive(Clone, Debug)]
pub struct MeshFunc {
    pub verts: Vec<VertexUnion>,
    /// Indices of triangles (taken as every 3 values).
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
                a @ _ => a.clone(),
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
              T: IntoIterator<Item=U>
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
    ///
    /// This returns a tuple of (new `MeshFunc`, new `arg_vals`), where
    /// `arg_vals[i]` is the new index of `self.verts[i]` in the
    /// returned `MeshFunc`.
    pub fn connect<T, U>(&self, children: T) -> (MeshFunc, Vec<Vec<usize>>)
        where U: Borrow<MeshFunc>,
              T: IntoIterator<Item=(U, Vec<usize>)>
    {
        // TODO: Clean up this description a bit
        // TODO: Clean up Vec<usize> stuff

        // Copy body vertices & faces:
        let mut verts: Vec<VertexUnion> = self.verts.clone();
        let mut faces = self.faces.clone();

        let mut remaps: Vec<Vec<usize>> = vec![];

        for (child_, arg_vals) in children {
            let child = child_.borrow();

            // 'offset' corresponds to the position in 'verts' at
            // which we begin appending everything in 'child.verts'.
            let offset = verts.len();

            // 'remap[i]' - if 'child.verts[i]' is a Vertex, not an Arg -
            // will contain the index in 'verts' that this vertex was
            // copied to.  (This is needed because below we copy only
            // the Vertex elements, not the Arg ones.)
            let mut remap = vec![0; child.verts.len()];
            let mut j = 0;

            // Like mentioned, copy just the Vertex in 'child.verts':
            verts.extend(child.verts.iter().enumerate().filter_map(|(i,v)| {
                match v {
                    VertexUnion::Vertex(_) => {
                        // TODO: A less-side-effectful way?
                        remap[i] = offset + j;
                        j += 1;
                        Some(v.clone())
                    },
                    VertexUnion::Arg(_) => None,
                }
            }));
            // So, note that:
            // 1st Vertex in 'child.verts' went to 'verts[offset + 0]'.
            // 2nd Vertex in 'child.verts' went to 'verts[offset + 1]'.
            // 3rd Vertex in 'child.verts' went to 'verts[offset + 2]'.
            // and so forth.
            // Since this skips Arg elements, we used 'remap' to
            // store the mapping: 'child.verts[i]' was copied to
            // 'verts[remap[i]].'

            // We then use 'remap' below to update the indices in
            // 'child.faces' to vertices in the returned mesh.
            //
            // Also, we didn't copy Arg elements from 'child.verts', but
            // we do use them below - if 'child.faces' makes a
            // reference to this Arg, we resolve it with the help of
            // 'arg_vals', passed in by the caller.
            faces.extend(child.faces.iter().map(|n| {
                let f = match child.verts[*n] {
                    VertexUnion::Vertex(_) => remap[*n],
                    VertexUnion::Arg(m) => arg_vals[m],
                };
                if f >= verts.len() {
                    panic!("face >= num_verts")
                }
                f
            }));

            // Since the caller might need this same remapping (there
            // could be an Arg that referred to these vertices), we
            // save this information to return it later:
            remaps.push(remap);
        }

        let m = MeshFunc {
            verts: verts,
            faces: faces,
        };
        (m, remaps)
    }
}
