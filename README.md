# This needs a title

## Highest priority:

- Consider trampolining `to_mesh`.  My call stack seems needlessly
  deep in spots.  Can I make tail-recursive?

## Important but less critical:

- Grep for all TODOs in code, really.
- Look at everything in README.md in automata_scratch.
- Implement some of the tougher examples from the above too, e.g. the
  triple nested spiral
- Lots of Rust-kosher refactoring (once I understand Rust better)
- Actual Rust-style docs!

## If I'm bored:

- See `curve_horn_start` comments; can I elegantly solve this issue of
  how to connect an exit vertex multiple places?
- Fix links in tri_mesh docs that use relative paths & do a PR?
- Look in https://www.nalgebra.org/quick_reference/# for "pour
  obtain".  Can I fix this somehow?  Looks like a French-ism that made
  its way in.
