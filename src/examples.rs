use nalgebra::*;
//pub mod examples;

use crate::openmesh::{OpenMesh, Tag, Mat4, Vertex, vertex};
use crate::rule::{Rule, RuleEval, Child};
use crate::prim;
use crate::util;

struct CurveHorn {
    seed: Vec<Vertex>,
    id_xform: Mat4,
    flip180: Mat4,
    incr: Mat4,
}

impl CurveHorn {

    fn init() -> (CurveHorn, Rule<CurveHorn>) {
        let y = &Vector3::y_axis();
        let c = CurveHorn {
            seed: vec![
                vertex(-0.5, -0.5, 0.0),
                vertex(-0.5,  0.5, 0.0),
                vertex( 0.5,  0.5, 0.0),
                vertex( 0.5, -0.5, 0.0),
            ],
            id_xform: nalgebra::geometry::Transform3::identity().to_homogeneous(),
            flip180: nalgebra::geometry::Rotation3::from_axis_angle(
                &nalgebra::Vector3::y_axis(),
                std::f32::consts::PI).to_homogeneous(),
            incr: geometry::Rotation3::from_axis_angle(y, 0.1).to_homogeneous() *
                Matrix4::new_scaling(0.95) *
                geometry::Translation3::new(0.0, 0.0, 0.2).to_homogeneous(),
        };
        (c, Rule { eval: Self::start })
    }
    
    fn start(&self) -> RuleEval<Self> {
        RuleEval {
            geom: OpenMesh {
                verts: self.seed.clone(),
                faces: vec![],
            },
            final_geom: prim::empty_mesh(),
            children: vec![
                Child {
                    rule: Rule { eval: Self::recur },
                    xf: self.id_xform,
                    vmap: vec![0,1,2,3],
                },
                Child {
                    rule: Rule { eval: Self::recur },
                    xf: self.flip180,
                    vmap: vec![3,2,1,0],
                },
            ],
        }
    }

    fn recur(&self) -> RuleEval<Self> {

        let verts = self.seed.clone();
        let next_verts: Vec<Vertex> = verts.iter().map(|v| self.incr * v).collect();
        
        let geom = OpenMesh {
            verts: next_verts.clone(),
            faces: vec![
                // The below is just connecting two groups of 4 vertices
                // each, straight across and then to the next.
                Tag::Body(1), Tag::Parent(0), Tag::Body(0),
                Tag::Parent(1), Tag::Parent(0), Tag::Body(1),
                Tag::Body(2), Tag::Parent(1), Tag::Body(1),
                Tag::Parent(2), Tag::Parent(1), Tag::Body(2),
                Tag::Body(3), Tag::Parent(2), Tag::Body(2),
                Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
                Tag::Body(0), Tag::Parent(3), Tag::Body(3),
                Tag::Parent(0), Tag::Parent(3), Tag::Body(0),
                // TODO: I should really generate these, not hard-code them.
            ],
        };

        // TODO: This could be made slightly nicer by taking it to a peak
        // instead of just flattening it in XY, but this is a pretty minor
        // change.
        let final_geom = OpenMesh {
            verts: vec![],
            faces: vec![
                Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
                Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
            ],
        };
        
        RuleEval{
            geom: geom,
            final_geom: final_geom,
            children: vec![
                Child {
                    rule: Rule { eval: Self::recur },
                    xf: self.incr,
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    }
}

struct CubeThing {
}

impl CubeThing {

    fn init() -> (CubeThing, Rule<CubeThing>) {
        (CubeThing {}, Rule { eval: Self::rec })
    }
    
    fn rec(&self) -> RuleEval<Self> {

        let mesh = prim::cube();

        // Quarter-turn in radians:
        let qtr = std::f32::consts::FRAC_PI_2;

        let y = &Vector3::y_axis();
        let z = &Vector3::z_axis();
        
        // Each element of this turns to a branch for the recursion:
        let turns: Vec<Mat4> = vec![
            geometry::Transform3::identity().to_homogeneous(),
            geometry::Rotation3::from_axis_angle(y, qtr).to_homogeneous(),
            geometry::Rotation3::from_axis_angle(y, qtr * 2.0).to_homogeneous(),
            geometry::Rotation3::from_axis_angle(y, qtr * 3.0).to_homogeneous(),
            geometry::Rotation3::from_axis_angle(z, qtr).to_homogeneous(),
            geometry::Rotation3::from_axis_angle(z, -qtr).to_homogeneous(),
        ];

        let gen_rulestep = |rot: &Mat4| -> Child<Self> {
            let m: Mat4 = rot *
                Matrix4::new_scaling(0.5) *
                geometry::Translation3::new(6.0, 0.0, 0.0).to_homogeneous();
            Child {
                rule: Rule { eval: Self::rec },
                xf: m,
                vmap: vec![],
            }
        };

        RuleEval {
            geom: mesh,
            final_geom: prim::empty_mesh(),
            children: turns.iter().map(gen_rulestep).collect(),
        }
    }
}

struct RamHorn {
}

impl RamHorn {

    fn init() -> (RamHorn, Rule<RamHorn>) {
        (RamHorn{}, Rule { eval: Self::start })
    }
    
    // Conversion from Python & automata_scratch
    fn start(&self) -> RuleEval<Self> {
        let opening_xform = |i| {
            let r = std::f32::consts::FRAC_PI_2 * i;
            ((geometry::Rotation3::from_axis_angle(
                &nalgebra::Vector3::z_axis(), r).to_homogeneous()) *
             geometry::Translation3::new(0.25, 0.25, 1.0).to_homogeneous() *
             Matrix4::new_scaling(0.5) *
             geometry::Translation3::new(0.0, 0.0, -1.0).to_homogeneous())
        };
        RuleEval {
            geom: OpenMesh {
                verts: vec![
                    // 'Top' vertices:
                    vertex(-0.5, -0.5, 1.0),  //  0 (above 9)
                    vertex(-0.5,  0.5, 1.0),  //  1 (above 10)
                    vertex( 0.5,  0.5, 1.0),  //  2 (above 11)
                    vertex( 0.5, -0.5, 1.0),  //  3 (above 12)
                    // Top edge midpoints:
                    vertex(-0.5,  0.0, 1.0),  //  4 (connects 0-1)
                    vertex( 0.0,  0.5, 1.0),  //  5 (connects 1-2)
                    vertex( 0.5,  0.0, 1.0),  //  6 (connects 2-3)
                    vertex( 0.0, -0.5, 1.0),  //  7 (connects 3-0)
                    // Top middle:
                    vertex( 0.0,  0.0, 1.0),  //  8
                    // 'Bottom' vertices:
                    vertex(-0.5, -0.5, 0.0),  //  9
                    vertex(-0.5,  0.5, 0.0),  // 10
                    vertex( 0.5,  0.5, 0.0),  // 11
                    vertex( 0.5, -0.5, 0.0),  // 12
                ],
                faces: vec![
                    // bottom face:
                    Tag::Body(9), Tag::Body(10), Tag::Body(11),
                    Tag::Body(9), Tag::Body(11), Tag::Body(12),
                    // two faces straddling edge from vertex 0:
                    Tag::Body(9), Tag::Body(0), Tag::Body(4),
                    Tag::Body(9), Tag::Body(7), Tag::Body(0),
                    // two faces straddling edge from vertex 1:
                    Tag::Body(10), Tag::Body(1), Tag::Body(5),
                    Tag::Body(10), Tag::Body(4), Tag::Body(1),
                    // two faces straddling edge from vertex 2:
                    Tag::Body(11), Tag::Body(2), Tag::Body(6),
                    Tag::Body(11), Tag::Body(5), Tag::Body(2),
                    // two faces straddling edge from vertex 3:
                    Tag::Body(12), Tag::Body(3), Tag::Body(7),
                    Tag::Body(12), Tag::Body(6), Tag::Body(3),
                    // four faces from edge (0,1), (1,2), (2,3), (3,0):
                    Tag::Body(9), Tag::Body(4), Tag::Body(10),
                    Tag::Body(10), Tag::Body(5), Tag::Body(11),
                    Tag::Body(11), Tag::Body(6), Tag::Body(12),
                    Tag::Body(12), Tag::Body(7), Tag::Body(9),
                ],
            },
            final_geom: prim::empty_mesh(),
            children: vec![
                Child {
                    rule: Rule { eval: Self::ram_horn },
                    xf: opening_xform(0.0),
                    vmap: vec![5,2,6,8],
                },
                Child {
                    rule: Rule { eval: Self::ram_horn },
                    xf: opening_xform(1.0),
                    vmap: vec![4,1,5,8],
                },
                Child {
                    rule: Rule { eval: Self::ram_horn },
                    xf: opening_xform(2.0),
                    vmap: vec![7,0,4,8],
                },
                Child {
                    rule: Rule { eval: Self::ram_horn },
                    xf: opening_xform(3.0),
                    vmap: vec![6,3,7,8],
                },
                // TODO: These vertex mappings appear to be right.
                // Explain *why* they are right.
            ],
        }
    }

    fn ram_horn(&self) -> RuleEval<Self> {
        let v = Unit::new_normalize(Vector3::new(-1.0, 0.0, 1.0));
        let incr: Mat4 = geometry::Translation3::new(0.0, 0.0, 0.8).to_homogeneous() *
            geometry::Rotation3::from_axis_angle(&v, 0.3).to_homogeneous() *
            Matrix4::new_scaling(0.9);
        let seed = vec![
            vertex(-0.5, -0.5, 1.0),
            vertex(-0.5,  0.5, 1.0),
            vertex( 0.5,  0.5, 1.0),
            vertex( 0.5, -0.5, 1.0),
        ];
        let next = seed.iter().map(|v| incr * v).collect();
        let geom = OpenMesh {
            verts: next,
            faces: vec![
                Tag::Body(1), Tag::Parent(0), Tag::Body(0),
                Tag::Parent(1), Tag::Parent(0), Tag::Body(1),
                Tag::Body(2), Tag::Parent(1), Tag::Body(1),
                Tag::Parent(2), Tag::Parent(1), Tag::Body(2),
                Tag::Body(3), Tag::Parent(2), Tag::Body(2),
                Tag::Parent(3), Tag::Parent(2), Tag::Body(3),
                Tag::Body(0), Tag::Parent(3), Tag::Body(3),
                Tag::Parent(0), Tag::Parent(3), Tag::Body(0),
            ],
        };
        let final_geom = OpenMesh {
            verts: vec![],
            faces: vec![
                Tag::Parent(0), Tag::Parent(2), Tag::Parent(1),
                Tag::Parent(0), Tag::Parent(3), Tag::Parent(2),
            ],
        };
        RuleEval {
            geom: geom,
            final_geom: final_geom,
            children: vec![
                Child {
                    rule: Rule { eval: Self::ram_horn },
                    xf: incr,
                    vmap: vec![0,1,2,3],
                },
            ],
        }
    }
}

struct Twist {
    seed: Vec<Vertex>,
    seed_sub: Vec<Vertex>,
    dx0: f32,
    dy: f32,
    ang: f32,
    count: usize,
    subdiv: usize,
}

impl Twist {

    pub fn init() -> (Twist, Rule<Twist>) {
        let subdiv = 2;
        let xf = geometry::Rotation3::from_axis_angle(&Vector3::x_axis(), -0.7).to_homogeneous();
        let seed = vec![
            vertex(-0.5,  0.0, -0.5),
            vertex( 0.5,  0.0, -0.5),
            vertex( 0.5,  0.0,  0.5),
            vertex(-0.5,  0.0,  0.5),
        ].iter().map(|v| xf * v).collect();
        let seed_sub = util::subdivide_cycle(&seed, subdiv);
        let t = Twist {
            dx0: 2.0,
            dy: 0.1,
            ang: 0.1,
            count: 4,
            seed: seed,
            seed_sub: seed_sub,
            subdiv: subdiv,
        };
        (t, Rule { eval: Self::start })
    }
    
    // Meant to be a copy of twist_from_gen from Python & automata_scratch
    pub fn start(&self) -> RuleEval<Twist> {

        let n = self.seed_sub.len();
        
        // Quarter-turn in radians:
        let qtr = std::f32::consts::FRAC_PI_2;
        let y = &Vector3::y_axis();
        let xform = |i| {
            (geometry::Rotation3::from_axis_angle(y, qtr * (i as f32)).to_homogeneous() *
             geometry::Translation3::new(self.dx0, 0.0, 0.0).to_homogeneous())
        };
        
        // First generate 'count' children, each one shifted/rotated
        // differently:
        let children: Vec<Child<Twist>> = (0..self.count).map(|i| {
            let xf = xform(i);
            Child {
                rule: Rule { eval: Self::recur },
                xf: xf,
                vmap: ((n+1)*i..(n+1)*(i+self.count)).collect(), // N.B.
                // note n+1, not n. the +1 is for the centroid below
            }
        }).collect();

        // Use byproducts of this to make 'count' copies of 'seed' with
        // this same transform:
        let meshes = children.iter().map(|child| {
            let mut vs = self.seed_sub.iter().map(|v| child.xf * v).collect();
            // and in the process, generate faces for these seeds:
            let (centroid, f) = util::connect_convex(&vs, false);
            vs.push(centroid);
            OpenMesh { verts: vs, faces: f }
        });
        
        RuleEval {
            geom: OpenMesh::append(meshes),
            final_geom: prim::empty_mesh(),
            children: children,
        }
    }

    pub fn recur(&self) -> RuleEval<Twist> {
        let y = &Vector3::y_axis();
        let incr = geometry::Translation3::new(-self.dx0, 0.0, 0.0).to_homogeneous() *
            geometry::Rotation3::from_axis_angle(y, self.ang).to_homogeneous() *
            geometry::Translation3::new(self.dx0, self.dy, 0.0).to_homogeneous();
        
        let seed_orig = self.seed.iter().map(|v| incr * v).collect();
        let seed_sub = util::subdivide_cycle(&seed_orig, self.subdiv);
        let n = seed_sub.len();

        let (vc, faces) = util::connect_convex(&seed_sub, true);
        
        RuleEval {
            geom: OpenMesh {
                verts: seed_sub,
                faces: util::parallel_zigzag_faces(n),
            },
            final_geom: OpenMesh { verts: vec![vc], faces },
            children: vec![
                Child {
                    rule: Rule { eval: Self::recur },
                    xf: incr,
                    vmap: (0..n).collect(),
                },
            ],
        }
    }
}

pub fn main() {

    {
        let vs = vec![
            vertex(-0.5,  0.0, -0.5),
            vertex( 0.5,  0.0, -0.5),
            vertex( 0.5,  0.0,  0.5),
            vertex(-0.5,  0.0,  0.5),
        ];
        let vs2 = util::subdivide_cycle(&vs, 2);
        println!("vs={:?}", vs);
        println!("vs2={:?}", vs2);
    }

    fn run_test<A>((a, r): (A, Rule<A>), iters: u32, name: &str) {
        println!("Running {}...", name);
        let (mesh, nodes) = r.to_mesh(&a, iters);
        println!("Evaluated {} rules", nodes);
        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        mesh.write_stl_file(&fname).unwrap();
    }

    fn run_test_iter<A>((a, r): (A, Rule<A>), iters: usize, name: &str) {
        println!("Running {}...", name);
        let (mesh, nodes) = r.to_mesh_iter(&a, iters);
        println!("Evaluated {} rules", nodes);
        let fname = format!("{}.stl", name);
        println!("Writing {}...", fname);
        mesh.write_stl_file(&fname).unwrap();
    }
    
    /*
    run_test(CubeThing::init(), Rule { eval: CubeThing::rec }, 3, "cube_thing");
    // this can't work on its own because the resultant OpenMesh still
    // has parent references:
    //run_test(Rule { eval: recur }, 100, "curve_horn_thing");
    run_test(CurveHorn::init(), Rule { eval: CurveHorn::start }, 100, "curve_horn2");
    run_test(RamHorn::init(), Rule { eval: RamHorn::start }, 200, "ram_horn");
    run_test(Twist::init(), Rule { eval: Twist::start }, 200, "twist");
    */

    run_test_iter(CubeThing::init(), 3, "cube_thing2");
    run_test_iter(CurveHorn::init(), 100, "curve_horn2_iter");
    run_test_iter(RamHorn::init(), 100, "ram_horn2");
    // TODO: If I increase the above from 100 to ~150, Blender reports
    // that the very tips are non-manifold.  I am wondering if this is
    // some sort of numerical precision issue.
    run_test_iter(Twist::init(), 100, "twist2");
}
