# This needs a title

This work was started as an attempt to make meshes in a more
"generative" style, described by recursive grammars and
replacement rules.  One goal was to make it easy to produce
manifold meshes by following certain rules, and do so in a
"correct-by-construction" manner rather than by having to
patch up or subdivide the meshes in post-processing.

These grammars by their nature worked in discrete steps,
but at one point I tried (unsuccessfully) to extend this
system to working in a more continuous and parametric
way.

I also ran into problems anytime I wanted to produce
meshes in a way that was more "refining" than "generative".
They're not completely distinct. However, the specific issue
I ran into is that the rules were explicitly designed around
'child' rules never being able to modify topology of geometry
from a 'parent' rule, besides being able to connect to its
vertices - and sometimes the "refining" part of things
required this in order to work right.

I've also disliked how much my model felt like it tied me
down to the "triangle mesh" representation. I haven't
found a good way to build up higher-level representations
to modularise and compose - but haven't given up yet on
this.  In some sense it is a conflict of goals because
the aim was correct-by-construction triangle meshes.

Also, I did this in order to learn the Rust language, and I
repeatedly kept bumping into the conclusion that Rust was
just not the right language for this. I was in need of things
like closures and first-class functions and I neglected to
consider how much those assume the presence of garbage
collection. Really, I wanted a Lisp, and then the presence of
a REPL would have been another bonus.

If I actually understood my goals enough to put better
constraints on my model, Rust probably would have been fine.
As it stands now, the lack of clarity in both my theory
and in my implementation is a far bigger issue than anything
related to Rust.

## Highest priority:

- Figure out the crash bug in `vec_indexed!` if I put a Vertex
  *after* an Arg.
- Just scrap `parametric_mesh` as much as possible and use existing
  tools (e.g. OpenSubdiv) because this DCEL method is just painful for
  what it is and I have some questions on how it can even work
  theoretically.
- Connect up the `parametric_mesh` stuff that remains, and worry about
  perfect meshes later.
- Get identical or near-identical meshes to `ramhorn_branch` from
  Python.  (Should just be a matter of tweaking parameters.)
- Look at performance.
  - Start at `to_mesh_iter()`. The cost of small appends/connects
    seems to be killing performance.
  - `connect()` is a big performance hot-spot: 85% of total time in
    one test, around 51% in `extend()`, 33% in `clone()`. It seems
    like I should be able to share geometry with the `Rc` (like noted
    above), defer copying until actually needed, and pre-allocate the
    vector to its size (which should be easy to compute).
- See `automata_scratch/examples.py` and implement some of the tougher
  examples.
  - `twisty_torus`, `spiral_nested_2`, & `spiral_nested_3` are all
    that remain.  To do them, I need to compose transformations (not
    in the matrix sense), but I also probably need to produce
    RuleEvals which always have `xf` of identity transformation since
    the Python code does not 'inherit' transforms unless I tell it to.

## Important but less critical:

- Docs on modules
- Compute global scale factor, and perhaps pass it to a rule (to
  eventually be used for, perhaps, adaptive subdivision).  Note that
  one can find the scale factors by taking the length of the first 3
  columns of the transform matrix (supposedly).
- swept-isocontour stuff from
  `/mnt/dev/graphics_misc/isosurfaces_2018_2019/spiral*.py`.  This
  will probably require that I figure out parametric curves
- Make an example that is more discrete-automata, less
  approximation-of-space-curve.

- Catch-alls:
  - Grep for all TODOs in code, really.
  - Look at everything in `README.md` in `automata_scratch`.

## If I'm bored:

- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
- Multithread!  This looks very task-parallel anywhere that I branch.
- Would being able to name a rule node (perhaps conditionally under
  some compile-time flag) help for debugging?
- Use an actual logging framework.
- Take a square.  Wrap it around to a torus. Now add a twist (about
  the axis that is normal to the square). This is simple, but it looks
  pretty cool.
- How can I take tangled things like the cinquefoil and produce more
  'iterative' versions that still weave around?

## Research Areas

- Can I use automatic differentiation in any way here to avoid the
  numerical annoyances?
- [Geometry and Algorithms for Computer Aided Design (Hartmann)](https://www2.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf)
- https://en.wikipedia.org/wiki/Surface_triangulation
- https://www.cs.cmu.edu/~quake/triangle.html

## Reflections & Quick Notes

- Generalizing to space curves moves this away from the "discrete
  automata" roots, but it still ends up needing the machinery I made
  for discrete automata.
- If you *pre* multiply a transformation: you are transforming the
  entire global space.  If you *post* multiply: you are transforming
  the current local space. 
