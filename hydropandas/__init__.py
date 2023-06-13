# flake8: noqa
import logging

from .extensions import geo as _geo
from .extensions import gwobs as _gwobs
from .extensions import plots as _plots
from .extensions import stats as _stats
from .io.fews import get_fews_pid
from .obs_collection import (
    ObsCollection,
    read_bro,
    read_bronhouderportaal_bro,
    read_dino,
    read_excel,
    read_fews,
    read_imod,
    read_knmi,
    read_menyanthes,
    read_modflow,
    read_pickle,
    read_waterinfo,
    read_wiski,
)
from .observation import (
    EvaporationObs,
    GroundwaterObs,
    MeteoObs,
    ModelObs,
    Obs,
    PrecipitationObs,
    WaterlvlObs,
    WaterQualityObs,
)
from .version import __version__

logging.getLogger("hydropandas").addHandler(logging.NullHandler())
