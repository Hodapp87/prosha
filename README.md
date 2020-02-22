# This needs a title

## Highest priority:

- See `ram_horn_*` and `curve_horn_*` TODOs: these are both solved by
  some sort of parent-vertex-mapping layer.  The only other way around
  this that I need is to require that rule functions exist in
  explictly separate forms.
- Consider trampolining `to_mesh`.  My call stack seems needlessly
  deep in spots.  Can I make tail-recursive?

## Important but less critical:

- Grep for all TODOs in code, really.
- Look at everything in README.md in automata_scratch.
- Implement some of the tougher examples from the above too, e.g. the
  triple nested spiral.  See `examples.py`.
- Actual Rust-style docs!

## If I'm bored:

- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
