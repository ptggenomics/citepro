import sys

from . import recipe
from . import io
from . import pp
from . import nbgui



sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["recipe", "io", "pp", "nbgui"]})

__all__ = ["recipe", "io", "pp", "nbgui"]
