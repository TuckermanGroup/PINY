"""Python tools and wrappers for the PINY molecular simulations package."""

# if the interface is not built, fail nicely
try:
    from .piny import PINYEmpirical
except ImportError:
    pass

# tools are pure Python, import
from . import tools
