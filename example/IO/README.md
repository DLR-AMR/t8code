# I/O Examples

This folder contains examples related to Input/Output operations in `t8code`. These examples demonstrate how to read and write `t8code`'s main data structures, such as the `t8_cmesh` and `t8_forest`, to and from disk.

## Subfolders

- [**cmesh/**](cmesh/README.md): Examples showing how to read a coarse mesh from a file and write it back out.
- [**forest/**](forest/README.md): Examples showing how to save an adapted forest to a file and load it back. This is useful for checkpointing and restarting simulations.

## Main Purpose

These examples are essential for understanding how to manage data persistence in `t8code`. They provide a starting point for integrating `t8code` into workflows that involve pre-existing mesh files or require saving simulation states.
