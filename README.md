# This needs a title

## Highest priority:

- Clean up `ramhorn_branch` because it's fugly.
- See `automata_scratch/examples.py` and implement some of the tougher
  examples.
  - `spiral_nested_2` & `spiral_nested_3` (how to compose
    efficiently?)
  - `twisty_torus`

## Important but less critical:

- Elegance & succinctness (my recent closure work may help with this):
  - Why must I repeat myself so much in these definitions?
  - What patterns can I factor out?  I do some things regularly, like:
    the clockwise boundaries, the zigzag connections
  - Procedural macro to shorten this `Tag::Parent`, `Tag::Body`
    nonsense - and perhaps force to groups of 3?
- Docs on modules
- Grep for all TODOs in code, really.
- Look at performance.  Can I save on copies of geometry by using
  `Rc<OpenMesh>` or the like?  In many cases I have nothing but copied
  geometry.  Can I pre-allocate vectors instead of
  extending/appending?
  - `connect()` is a big performance hot-spot: 85% of total time in
    one test, around 51% in `extend()`, 33% in `clone()`. It seems
    like I should be able to share geometry with the `Rc` (like noted
    above), defer copying until actually needed, and pre-allocate the
    vector to its size (which should be easy to compute).
  - The cost of small appends/connects seems to be killing
    performance.
- Look at everything in `README.md` in `automata_scratch`.
- Use an actual logging framework.
- Migrate tests to... well... actual tests.
- I am starting to see a pattern emerge in how I have to modularize
  things around closures.  What can a macro do for me here?
- swept-isocontour stuff from
  `/mnt/dev/graphics_misc/isosurfaces_2018_2019/spiral*.py`

## If I'm bored:

- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
- Multithread!  This looks very task-parallel anywhere that I branch.
- Would being able to name a rule node (perhaps conditionally under
  some compile-time flag) help for debugging?
