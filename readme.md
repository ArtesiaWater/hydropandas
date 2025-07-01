<img src="/docs/_static/Artesia_logo.jpg" alt="Artesia" width="200" align="right">

[![PyPi](https://img.shields.io/pypi/v/hydropandas.svg)](https://pypi.python.org/pypi/hydropandas)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/hydropandas)](https://pypi.python.org/pypi/hydropandas)
[<img src="https://github.com/codespaces/badge.svg" height="20">](https://codespaces.new/ArtesiaWater/hydropandas?quickstart=1)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![hydropandas](https://github.com/ArtesiaWater/hydropandas/actions/workflows/on_pr_master.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/on_pr_master.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Coverage)
[![Documentation Status](https://readthedocs.org/projects/hydropandas/badge/?version=latest)](https://hydropandas.readthedocs.io/en/latest/?badge=latest)

# HydroPandas

Hydropandas is a Python package for reading, analyzing and writing
(hydrological) timeseries.

## Reading

The HydroPandas package provides convenient read functions from various sources.
The table below lists all API-accessible sources. Click a link in the first column
for the documentation. The "API available" column indicates current availability
(updated weekly).


| source          | observations                       | API available | location             |
|-----------------|------------------------------------|---------------|----------------------|
| [BRO](https://hydropandas.readthedocs.io/en/stable/examples/01_groundwater_observations.html) | Groundwater                  | [![BRO](https://github.com/ArtesiaWater/hydropandas/actions/workflows/bro.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/bro.yml) | Netherlands          |
| [KNMI](https://hydropandas.readthedocs.io/en/stable/examples/02_knmi_observations.html) | Meteorological                 | [![KNMI](https://github.com/ArtesiaWater/hydropandas/actions/workflows/knmi.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/knmi.yml) | Netherlands          |
| [Lizard](https://hydropandas.readthedocs.io/en/stable/examples/06_lizard.html) | Groundwater                  | [![Lizard](https://github.com/ArtesiaWater/hydropandas/actions/workflows/lizard.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/lizard.yml) | Netherlands (Vitens) |
| [Waterconnect](https://hydropandas.readthedocs.io/en/stable/examples/09_waterconnect.html) | Groundwater                  | [![Waterconnect](https://github.com/ArtesiaWater/hydropandas/actions/workflows/waterconnect.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/waterconnect.yml) | South Australia      |
| [Waterinfo](https://hydropandas.readthedocs.io/en/stable/examples/08_waterinfo.html) | Surface water quantity and quality | [![Waterinfo](https://github.com/ArtesiaWater/hydropandas/actions/workflows/waterinfo.yml/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions/workflows/waterinfo.yml) | Netherlands          |
---

Some sources also provide files readable by HydroPandas.


| source          | observations                       | file format          | location             |
|-----------------|------------------------------------|----------------------|----------------------|
| [BRO](https://hydropandas.readthedocs.io/en/stable/examples/01_groundwater_observations.html) | Groundwater                  | xml          | Netherlands          |
| [DINO](https://hydropandas.readthedocs.io/en/stable/examples/01_groundwater_observations.html) | Groundwater / surface water                  | csv          | Netherlands          |
| [FEWS](https://hydropandas.readthedocs.io/en/stable/examples/07_fews.html) | Groundwater / surface water                  | xml          | Netherlands          |
| [KNMI](https://hydropandas.readthedocs.io/en/stable/examples/02_knmi_observations.html) | Meteorological                 | txt          | Netherlands          |
| [Pastastore](https://hydropandas.readthedocs.io/en/stable/examples/03_hydropandas_and_pastas.html) | Time series models                  | NA      | NA      |
| [Waterinfo](https://hydropandas.readthedocs.io/en/stable/examples/08_waterinfo.html) | Surface water quantity and quality | csv / zip          | Netherlands          |
| Wiski (no docs available)                | Groundwater | csv          | Netherlands          |
---
## Install

Install the module with pip:

`pip install hydropandas`

For some functionality additional packages are required. Install all optional packages:

`pip install hydropandas[full]`

For installing in development mode, clone the repository and install by
typing `pip install -e .[full]` from the module root directory.

## Documentation

-   Documentation is provided on the dedicated website
    [hydropandas.readthedocs.io](https://hydropandas.readthedocs.io/en/stable/)
-   Examples are available in the [examples directory on the documentation website](https://hydropandas.readthedocs.io/en/stable/examples.html)
-   View and edit the example notebooks of hydropandas in
    [GitHub Codespaces](https://codespaces.new/hydropandas/hydropandas?quickstart=1)

## Get in touch

- Questions on HydroPandas ("How can I?") can be asked and answered on [Github Discussions](https://github.com/ArtesiaWater/hydropandas/discussions).
- Bugs, feature requests and other improvements can be posted as [Github Issues](https://github.com/ArtesiaWater/hydropandas/issues).
- Find out how to contribute to HydroPandas at our [Contribution page](https://hydropandas.readthedocs.io/en/stable/contribute.html).


## Structure

The HydroPandas package allows users to store a timeseries and metadata in a
single object (Obs class). Or store a collection of timeseries with metadata
in a single object (ObsCollection class). Both inheret from a pandas DataFrame
and are extended with custom methods and attributes related to hydrological timeseries.

### The Obs class

The Obs class holds the measurements and metadata for one timeseries. There are
currently 7 specific Obs classes for different types of measurements:

- GroundwaterObs: for groundwater measurements
- WaterQualityObs: for groundwater quality measurements
- WaterlvlObs: for surface water level measurements
- ModelObs: for "observations" from a MODFLOW model
- MeteoObs: for meteorological observations
- PrecipitationObs: for precipitation observations, subclass of MeteoObs
- EvaporationObs: for evaporation observations, subclass of MeteoObs

Each of these Obs classes is essentially a pandas DataFrame with additional
methods and attributes related to the type of measurement that it holds.
Each Obs object also contains specific methods to read data from specific sources.

### The ObsCollection class

The ObsCollection class hold the data for a collection of Obs classes, e.g. 
10 timeseries of the groundwater level in a certain area. The
ObsCollection is essentialy a pandas DataFrame in which each timeseries is stored
in a different row. Each row contains metadata (e.g. latitude and longitude
of the observation point) and the Obs object that holds the
measurements. It's recommended to use one ObsCollection per observation type — for 
example, group 10 GroundwaterObs in one collection and 5 PrecipitationObs in another.

More information on dealing with Obs and ObsCollection objects in [the documentation](https://hydropandas.readthedocs.io/en/stable/examples/00_hydropandas_objects.html)

## Authors

- Onno Ebbens, Artesia
- Ruben Caljé, Artesia
- Davíd Brakenhoff, Artesia
- Martin Vonk, Artesia
