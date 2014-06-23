This directory contains files with parameters of the q-SPC/FW water model and
an XYZ file with an equilibrated configuration.

The two subdirectories contain example build scripts for classical and PIMD
simulations of 64 water molecules in a bulk setup using the q-SPC/Fw water
potential. Simulation input directories are created by running the build script
with the name of the output directory as an argument:

./build.sh run

For the PIMD simulation, you also need to provide the number of beads as the
second argument:

./build.sh run 32

This requires that you have the PINY Python package available in your
PYTHONPATH. Refer to the main README.md for instructions.
