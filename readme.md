<img src="/docs/_static/Artesia_logo.jpg" alt="Artesia" width="200" align="right">

[![PyPi](https://img.shields.io/pypi/v/hydropandas.svg)](https://pypi.python.org/pypi/hydropandas)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/hydropandas)](https://pypi.python.org/pypi/hydropandas)
[<img src="https://github.com/codespaces/badge.svg" height="20">](https://codespaces.new/ArtesiaWater/hydropandas?quickstart=1)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![hydropandas](https://github.com/ArtesiaWater/hydropandas/workflows/hydropandas/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions?query=workflow%3Ahydropandas)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Coverage)
[![Documentation Status](https://readthedocs.org/projects/hydropandas/badge/?version=latest)](https://hydropandas.readthedocs.io/en/latest/?badge=latest)

# HydroPandas

Hydropandas is a Python package for reading, analyzing and writing
(hydrological) timeseries data.

## Introduction

The HydroPandas package allows users to store a timeseries and metadata in a
single object. This object inherits from a pandas DataFrame, with all its
wonderful features, and is extended with custom methods and attributes related
to hydrological timeseries.

The HydroPandas package also provides convenient read functions for Dutch hydrological data from:

- [BRO](https://www.broloket.nl)
- [DINO](https://www.dinoloket.nl)
- FEWS PI-XML
- [KNMI](https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script)
- [Lizard](https://vitens.lizard.net/)
- MODFLOW groundwater models
- IMOD groundwater models
- [Pastastore](https://github.com/pastas/pastastore)
- [Waterinfo](https://waterinfo.rws.nl/)
- WISKI csv files

## Install

Install the module with pip:

`pip install hydropandas`

HydroPandas requires `pandas`, `scipy`, `matplotlib`, `tqdm`, `requests` and `colorama`.

For some functionality additional packages are required:

- `geopandas`: for dealing with shapefiles
- `pastastore`: for reading or storing data from PastaStore
- `bokeh`, `branca`, `folium`: for interactive maps
- `flopy`: for reading data from MODFLOW models
- `xarray`: for loading data from REGIS

For installing in development mode, clone the repository and install by
typing `pip install -e .` from the module root directory.
For installing all the optional packages use `pip install -e .[full]`.

## Get in touch

- Questions on HydroPandas ("How can I?") can be asked and answered on [Github Discussions](https://github.com/ArtesiaWater/hydropandas/discussions).
- Bugs, feature requests and other improvements can be posted as [Github Issues](https://github.com/ArtesiaWater/hydropandas/issues).
- Find out how to contribute to HydroPandas at our [Contribution page](https://hydropandas.readthedocs.io/en/stable/contribute.html).

## Examples

Importing a groundwater time series from the BRO using the BRO-id and the tube number:

```python
import hydropandas as hpd
gw_bro = hpd.GroundwaterObs.from_bro("GMW000000041261", 1)
```

Or import all groundwater time series from the BRO within a certain extent:

```python
oc = hpd.read_bro(extent=(117850, 118180, 439550, 439900))
```

## The Obs class

The Obs class holds the measurements and metadata for one timeseries. There are
currently 5 specific Obs classes for different types of measurements:

- GroundwaterObs: for groundwater measurements
- WaterQualityObs: for groundwater quality measurements
- WaterlvlObs: for surface water level measurements
- ModelObs: for "observations" from a MODFLOW model
- MeteoObs: for meteorological observations
- PrecipitationObs: for precipitation observations, subclass of MeteoObs
- EvaporationObs: for evaporation observations, subclass of MeteoObs

Each of these Obs classes is essentially a pandas DataFrame with additional
methods and attributes related to the type of measurement that it holds.
Each Obs object also contain specific methods to read data from specific sources.

## The ObsCollection class

The ObsCollection class, as the name implies, represents a collection of Obs
classes, e.g. 10 timeseries of the groundwater level in a certain area. The
ObsCollection is also a pandas DataFrame in which each timeseries is stored
in a different row. Each row contains metadata (e.g. latitude and longitude
of the observation point) and the Obs object (DataFrame) that holds the
measurements. It is recommended to let an ObsCollection contain only one Obs
type, e.g. to create an ObsCollection for 10 GroundwaterObs, and a separate
ObsCollection for 5 PrecipitationObs.

Like the Obs class, the ObsCollection class contains a bunch of methods for
reading data from different sources. See the next section for supported data
sources.

## Authors

- Onno Ebbens, Artesia
- Ruben Caljé, Artesia
- Davíd Brakenhoff, Artesia
- Martin Vonk, Artesia
