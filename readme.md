<img src="/docs/_static/Artesia_logo.jpg" alt="Artesia" width="200" align="right">

[![PyPi](https://img.shields.io/pypi/v/hydropandas.svg)](https://pypi.python.org/pypi/hydropandas)
[![PyPi Supported Python Versions](https://img.shields.io/pypi/pyversions/hydropandas)](https://pypi.org/project/spei/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ArtesiaWater/hydropandas/master)

[![hydropandas](https://github.com/ArtesiaWater/hydropandas/workflows/hydropandas/badge.svg)](https://github.com/ArtesiaWater/hydropandas/actions?query=workflow%3Ahydropandas)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/c1b99f474bdc49b0a47e00e4e9f66c2f)](https://www.codacy.com/gh/ArtesiaWater/hydropandas/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ArtesiaWater/hydropandas&utm_campaign=Badge_Coverage)
[![Documentation Status](https://readthedocs.org/projects/hydropandas/badge/?version=latest)](https://hydropandas.readthedocs.io/en/latest/?badge=latest)

[![Format: isort](https://img.shields.io/badge/imports-isort-ef8336)](https://pycqa.github.io/isort/index.html)
[![Format: Black](https://img.shields.io/badge/code_style-black-black)](https://github.com/psf/black)
[![Linter: flake8](https://img.shields.io/badge/linter-flake8-yellowgreen)](https://flake8.pycqa.org/)
[![Linter: ruff](https://img.shields.io/badge/linter-ruff-red)](https://github.com/charliermarsh/ruff)

# hydropandas

The HydroPandas package allows users to store observation data at multiple locations in a single object (ObsCollection).
An ObsCollection is a pandas DataFrame extended with custom methods and attributes related to hydrological timeseries.
The HydroPandas package also provides convenient read functions for Dutch hydrological data from:
-   [BRO](https://www.broloket.nl)
-   [DINO](https://www.dinoloket.nl)
-   FEWS PI-XML
-   [KNMI](https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script)
-   MODFLOW groundwater models
-   IMOD groundwater models
-   [Pastastore](https://github.com/pastas/pastastore)
-   [waterinfo](https://waterinfo.rws.nl/)
-   WISKI csv files


## Installation

Install the module with pip:

`pip install hydropandas`

Hydropandas requires `numpy`, `scipy`, `matplotlib`, `pandas`, `requests` and `zeep`.

For some functionality additional packages are required:

-   `geopandas`: for dealing with shapefiles
-   `pastastore`: for reading or storing data from PastaStore
-   `bokeh`, `branca`, `folium`: for interactive maps
-   `flopy`: for reading data from MODFLOW models
-   `xarray`: for loading data from REGIS

For installing in development mode, clone the repository and install by
typing `pip install -e .` from the module root directory.
For installing all the optional packages use `pip install -e .[full]`

If you have trouble installing HydroPandas, refer to the
[Dependencies section](#dependencies) below.

## Example usage

Importing a single DINO csv file:

```python
import hydropandas as hpd
fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
gw = hpd.GroundwaterObs.from_dino(fname=fname)
```

Or for a zipfile:

```python
import hydropandas as hpd
dinozip = './tests/data/2019-Dino-test/dino.zip'
dino_gw = hpd.read_dino(dirname=dinozip,
			subdir='Grondwaterstanden_Put',
                        suffix='1.csv',
                        ObsClass=hpd.GroundwaterObs,
                        keep_all_obs=False)
```

## The Obs class

The Obs class holds the measurements and metadata for one timeseries. There are
currently 5 specific Obs classes for different types of measurements:

-   GroundwaterObs: for groundwater measurements
-   WaterQualityObs: for groundwater quality measurements
-   WaterlvlObs: for surface water level measurements
-   ModelObs: for "observations" from a MODFLOW model
-   MeteoObs: for meteorological observations
-   PrecipitationObs: for precipitation observations, subclass of MeteoObs
-   EvaporationObs: for evaporation observations, subclass of MeteoObs

Each of these Obs classes is essentially a pandas DataFrame with additional
methods and attributes related to the type of measurement that it holds.
The classes also contain specific methods to read data from specific sources.

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

## Dependencies

Hydropandas (indirectly) uses some packages that cannot be installed
automatically with `pip` on Windows. These packages are:

-   GDAL
-   Fiona
-   Shapely

If you do not have these packages already it is recommended to first try
installing them with `conda install <pkg>`. Otherwise, read the instructions
below how to install them manually.

Download the packages from [Christoph Gohlke's website](https://www.lfd.uci.edu/~gohlke/pythonlibs).
Use CTRL+F to find the download link on the page. Be sure to download the
correct version of the package. The Python version should match your Python
version. Also the architecture should match (i.e. 64bits vs 32bits).
For example:

-   GDAL-3.1.4-cp38-cp38-win_amd64.whl

This is the GDAL version for Python 3.8 (as can be seen from the cp38 in the
name), for 64-bits Python (as derived from the amd64 in the name).

Once you have downloaded the correct files, navigate to the directory in which
you saved your downloads. Now type the following commands (the order is
important):

1.  `pip install GDAL-3.1.4-cp38-cp38-win_amd64.whl`
2.  `pip install Fiona-1.8.17-cp38-cp38-win_amd64.whl`
3.  `pip install Shapely-1.7.1-cp38-cp38-win_amd64.whl`

After you've done this you can install hydropandas using
`pip install hydropandas`.

## Authors

-   Onno Ebbens, Artesia
-   Ruben Caljé, Artesia
-   Davíd Brakenhoff, Artesia
-   Martin Vonk, Artesia
