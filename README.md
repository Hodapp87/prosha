# This needs a title

## Highest priority:

- See `automata_scratch/examples.py` and implement some of the tougher
  examples.
  - `spiral_nested_2` & `spiral_nested_3` (how to compose
    efficiently?)
  - `twisty_torus`
  - `ram_horn_branch` - how do I pass depth in order to do this right?

## Important but less critical:

- Elegance & succinctness:
  - Why must I repeat myself so much in these definitions?
  - The notation for transforms is really cumbersome.  Some syntactic
    sugar might go far.
  - What patterns can I factor out?  I do some things regularly, like:
    the clockwise boundaries, the zigzag connections, the iterating over
    a `Vec<Vertex>` to transform each element and make another vector.
  - I have seen many of my bugs come from: all this arithmetic on
    indices.  I generate vertex maps more or less manually.
  - Helper method to transform `Vec<Vertex>` and such.  I do this
    often.
- Docs on modules
- Grep for all TODOs in code, really.
- Look at everything in README.md in automata_scratch.

## If I'm bored:

- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
- Multithread!  This looks very task-parallel anywhere that I branch.
