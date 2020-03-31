# This needs a title

## Highest priority:

- If my `closure_try2` branch seems to be working: start converting
  other things and cleaning everything up.  (`twist` is still ugly.
  Look at all my TODOs in it.)
- See `automata_scratch/examples.py` and implement some of the tougher
  examples.
  - `spiral_nested_2` & `spiral_nested_3` (how to compose
    efficiently?)
  - `twisty_torus`
  - `ram_horn_branch` - Can I pass depth via a closure?

## Important but less critical:

- Elegance & succinctness (my recent closure work may help with this):
  - Why must I repeat myself so much in these definitions?
  - The notation for transforms is really cumbersome.  Some syntactic
    sugar might go far.
  - What patterns can I factor out?  I do some things regularly, like:
    the clockwise boundaries, the zigzag connections, the iterating over
    a `Vec<Vertex>` to transform each element and make another vector.
- Docs on modules
- Grep for all TODOs in code, really.
- Look at performance.  Can I save on copies of geometry by using
  `Rc<OpenMesh>` or the like?  In many cases I have nothing but copied
  geometry.  Can I pre-allocate vectors instead of
  extending/appending?
- Look at everything in README.md in automata_scratch.

## If I'm bored:

- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
- Multithread!  This looks very task-parallel anywhere that I branch.
