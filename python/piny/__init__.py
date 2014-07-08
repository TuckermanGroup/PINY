"""Python tools and wrappers for the PINY molecular simulations package."""

from __future__ import absolute_import


# if the interface is not built, fail nicely
try:
    from .piny import PINYEmpirical
except ImportError:
    pass

# if ASE is available, add support for it
try:
    import ase
except ImportError:
    pass
else:
    del ase
    from . import ase

# tools are pure Python, import
from . import tools
