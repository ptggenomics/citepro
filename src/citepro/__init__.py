import sys

from . import recipe
from . import io
from . import pp



sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["recipe", "io", "pp"]})

__all__ = ["recipe", "io", "pp"]
