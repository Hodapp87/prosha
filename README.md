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
- Look at performance.
  - Start at `to_mesh_iter()`. The cost of small appends/connects
    seems to be killing performance.
  - `connect()` is a big performance hot-spot: 85% of total time in
    one test, around 51% in `extend()`, 33% in `clone()`. It seems
    like I should be able to share geometry with the `Rc` (like noted
    above), defer copying until actually needed, and pre-allocate the
    vector to its size (which should be easy to compute).
- Compute global scale factor, and perhaps pass it to a rule (to
  eventually be used for, perhaps, adaptive subdivision)
- swept-isocontour stuff from
  `/mnt/dev/graphics_misc/isosurfaces_2018_2019/spiral*.py`

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
