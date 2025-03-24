import sys
import logging
logger = logging.getLogger('citepro')
logger.setLevel(logging.INFO)

#try:
#    import cudf
#    use_gpu = True
#    logger.info("found cudf, try using cudf zero code change mode")
#    import cudf.pandas
#    cudf.pandas.install()
#except ModuleNotFoundError:
#    use_gpu = False
#    logger.info("rapids_singlecell not installed, fallback to CPU")

# rapids_single cell seems to have incompatible API calls with cuml. so this trial1 below failed.
#try:
#    from cuml.accel import install
#    install()
#except ModuleNotFoundError:
#    logger.info("cuml not installed, fallback to sk-learn")
    


from . import recipe
from . import io
from . import pp
from . import nbgui
from . import utils


sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["recipe", "io", "pp", "nbgui", "utils"]})

__all__ = ["recipe", "io", "pp", "nbgui", "utils"]
