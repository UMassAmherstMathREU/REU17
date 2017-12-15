# REU17
PythonSage3Dptts

# 3D Partition Library for SageMath

This repo contains a package for doing "box counting" calculations,
with skew plane partitions and reverse plane partitions, for the
purpose of calculating Donaldson-Thomas and Pandharipende-Thomas
invariants.  The code in ptdt_package defines these objects, and can
calculate the 1-leg case on both sides, with insertions.  It also
implements the Hillman-Grassl correspondence, which is relevant to
enumerating these objects.  This code has documentation and tests, and
matches the style of the sage packages for similar combinatorial
objects.

The rest of the repo is a bit scattered at the moment.  This includes
code we've used for experiments.  Code to calculate the full 3-leg
case on both the PT and DT sides is a work in progress.

## Credits
Gabriel Clara

James Hagborg

Augustus Ijams

Edward McCormick
