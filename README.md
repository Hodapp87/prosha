# This needs a title

## Highest priority:

- Continue to refine the 'barbs' example, which broke some new ground.
- Implement the continuous parametric transformations from 2020-05-07
  in my notes.  This will require some new abstractions.
- Try some non-deterministic examples.
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

- Elegance & succinctness:
  - Clean up `ramhorn_branch` because it's ugly.
  - What patterns can I factor out?  I do some things regularly, like:
    the clockwise boundaries, the zigzag connections.
  - Declarative macro to shorten this `Tag::Parent`, `Tag::Body`
    nonsense - and perhaps force to groups of 3?  Does this have any
    value, though, over just making helper functions like `p(...)` and
    `b(...)`?
  - I'm near certain a declarative macros can simplify some bigger
    things like my patterns with closures (e.g. the Y combinator like
    method for recursive calls).
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

- Fix links in tri_mesh docs that use relative paths & do a PR?
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

- When I have an iterated transform, that is basically transforming by
  M, MM=M^2, MMM=M^3, ..., and it seems to me that I should be able to
  compute its eigendecomposition and use this to compute fractional
  powers of the matrix.  Couldn't I then determine the continuous
  function I'm approximating by taking the `d/di (M^i)V` - i.e. the
  partial derivative of the result of transforming a vector `V` with
  `M^i`?  (See also:
  https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Functional_calculus
  and my 2020-04-20 paper notes.  My 2020-04-24 org notes have some
  things too - this relates to dynamical systems and eigenvalues.)
  Later note: I have a feeling I was dead wrong about a bunch of this.

## Reflections & Quick Notes

- My old Python version composed rules in the opposite order and I
  think this made things more complicated.  I didn't realize that I
  did it differently in this code, but it became much easier -
  particularly, more "inner" transformations are much easier to write
  because all that matters is that they work properly in the
  coordinate space they inherit.
- Generalizing to space curves moves this away from the "discrete
  automata" roots, but it still ends up needing the machinery I made
  for discrete automata.
- If you *pre* multiply a transformation: you are transforming the
  entire global space.  If you *post* multiply: you are transforming
  the current local space. 