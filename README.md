# This needs a title

## Highest priority:

- Continue fixing `to_mesh_iter`, which still doesn't yet handle
  branching because (for one thing) it never tracks depth properly in
  order to backtrack.
- Clean up `twist` - maybe make a struct or trait.
- Do transforms compose in the *reverse* of automata_scratch? This
  appears to be the case.

## Important but less critical:

- Why must I repeat myself so much in these definitions?
- The notation for transforms is really cumbersome.  Some syntactic
  sugar might go far.
- What patterns can I factor out?  I do some things regularly, like:
  the clockwise boundaries, the zigzag connections, the iterating over
  a `Vec<Vertex>` to transform each element and make another vector.
- Docs on modules
- Grep for all TODOs in code, really.
- Look at everything in README.md in automata_scratch.
- Implement some of the tougher examples from the above too, e.g. the
  triple nested spiral.  See `examples.py`.

## If I'm bored:

- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
