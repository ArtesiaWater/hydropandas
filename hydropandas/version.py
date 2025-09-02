from importlib import metadata
from sys import version as os_version

__version__ = "0.15.0"


def show_versions():
    """Method to print the versions of dependencies."""
    msg = (
        f"Python version      : {os_version}\n"
        f"Numpy version       : {metadata.version('numpy')}\n"
        f"Scipy version       : {metadata.version('scipy')}\n"
        f"Pandas version      : {metadata.version('pandas')}\n"
        f"Matplotlib version  : {metadata.version('matplotlib')}\n\n"
        f"Hydropandas version : {__version__}"
    )
    return print(msg)
