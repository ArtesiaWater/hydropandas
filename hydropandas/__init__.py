# flake8: noqa
import logging

from .data.fews_pid import get_fews_pid
from .extensions import geo as _geo
from .extensions import gwobs as _gwobs
from .extensions import plots as _plots
from .extensions import stats as _stats
from .obs_collection import (
    ObsCollection,
    read_bro,
    read_dino,
    read_fews,
    read_imod,
    read_knmi,
    read_menyanthes,
    read_modflow,
    read_waterinfo,
    read_wiski,
)
from .observation import (
    EvaporationObs,
    GroundwaterObs,
    WaterQualityObs,
    MeteoObs,
    ModelObs,
    Obs,
    PrecipitationObs,
    WaterlvlObs,
)
from .version import __version__

logging.getLogger("hydropandas").addHandler(logging.NullHandler())
