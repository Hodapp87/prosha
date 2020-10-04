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
way.  (See `parametric_mesh` and any DCEL code.)

I also ran into problems anytime I wanted to produce
meshes in a way that was more "refining" than "generative".
They're not completely distinct. However, the specific issue
I ran into is that the rules were explicitly designed around
'child' rules never being able to modify topology of geometry
from a 'parent' rule, besides being able to connect to its
vertices - and sometimes the "refining" part of things
required this in order to work right.

The problems with the parametric/continuous, and the
aforementioned "refining", were related. The issue is that
in order to get good meshes, I needed to be able to minimize
approximation error with the triangles and avoid triangles
with extreme angles, and there was seemingly no good way to
do this by incremental construction (like I was trying to
use elsewhere in my model) - and so its seems I just ended up
reinventing, badly, a lot of existing work with subdivision
and meshing.

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

On top of this, my implementation is pretty slow when it is
using a large number of rules each producing small geometry
(which is almost literally the only way it *can* be used
if you want to produce a fairly complex mesh). I did some
profiling some months ago that showed I was spending the
vast majority of my time in `extend()` and `clone()` for
`Vec` - and so I could probably see some huge performance
gains if I could simply pre-allocate vectors and share geometry
more. Also, I'm pretty sure this code does some very task-parallel
elements (e.g. anytime a rule branches), and multithreading should
be able to exploit this if I care.

If I actually understood my goals enough to put better
constraints on my model, Rust probably would have been fine.
As it stands now, the lack of clarity in both my theory
and in my implementation is a far bigger issue than anything
related to Rust.

## Highest priority:

- Fix `ramhorn_branch`.
- Once I've fixed that, see about a refactor that respects the
  same model, but involves much less ceremony and boilerplate.
- Look at performance.
  - Start at `to_mesh_iter()`. The cost of small appends/connects
    seems to be killing performance.
  - `connect()` is a big performance hot-spot: 85% of total time in
    one test, around 51% in `extend()`, 33% in `clone()`. It seems
    like I should be able to share geometry with the `Rc` (like noted
    above), defer copying until actually needed, and pre-allocate the
    vector to its size (which should be easy to compute).

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
  - Look at everything in `README.md` in `automata_scratch`,
    my old Python code from around 2019-09.

## If I'm bored:

- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
- Multithread!  This looks very task-parallel anywhere that I branch.
- Would being able to name a rule node (perhaps conditionally under
  some compile-time flag) help for debugging?
- Use an actual logging framework.
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
- Don't reinvent subdivision surfaces.