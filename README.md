PINY
====

Build
-----

Pick a suitable platform-specific subdirectory in the `compile` directory and use it either directly or make a copy to modify. The file `make_head.mk` holds all build settings. Run `make` and the resulting executable will be placed in `bin`.


Examples
---------

There are example simulations ready to run in the `examples` directory.


Python tools and interface
--------------------------

There are some Python tools that can help with preparing and post-processing a PINY simulation. They include command line utilities as well as a Python package. To use these, set your paths by sourcing the `env.sh` file in the `python` directory. For shells other than bash, set the required environment variables "by hand". For repeated use, it is probably a good idea to place these settings in your shell's RC file.

The tools require the numpy library.

Additionally, if you want to use the optional Python interface to PINY itself, you also need to have Cython. First, in one of the `compile` subdirectories, do `make libpiny` to link a shared library. Then, do `make` in the `python` directory to build the interface.

Acknowledgments:
----------------
Funding for the development of PINY_MD has come from a variety of sources, most recently, from grants from the National Science Foundation CHE-1301314 and CHE-1566085 
