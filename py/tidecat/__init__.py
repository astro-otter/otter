# get the version
from ._version import __version__


# explicitly set the package variable to ensure relative import work
__package__ = "tidecat"

# import important stuff
from .tde import TDE
from .tdecatalog import TDECatalog
