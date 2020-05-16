# This needs a title

## Highest priority:

- Implement the continuous parametric transformations from 2020-05-07
  in my notes.  This will require some new abstractions.
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

## Reflections & Quick Notes

- Generalizing to space curves moves this away from the "discrete
  automata" roots, but it still ends up needing the machinery I made
  for discrete automata.
- If you *pre* multiply a transformation: you are transforming the
  entire global space.  If you *post* multiply: you are transforming
  the current local space. 