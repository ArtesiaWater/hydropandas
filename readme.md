<img src="/doc/_static/Artesia_logo.jpg" alt="Artesia" width="200" align="right">

[![Build Status](https://travis-ci.com/ArtesiaWater/hydropandas.svg?branch=master)](https://travis-ci.com/ArtesiaWater/hydropandas)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3cfd43a484ac4591a9ce1c117d062bab)](https://www.codacy.com/manual/ArtesiaWater/hydropandas?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ArtesiaWater/hydropandas&amp;utm_campaign=Badge_Grade)

# hydropandas
The hydropandas module is a Python package for reading timeseries data into DataFrames. The basic idea behind the package is to allow users to manipulate data using all of the wonderful features included in pandas, but to also allow the addition of custom methods and attributes related to the timeseries. The hydropandas module extends pandas.DataFrame with extra functionality and stores metadata related to the type of measurements.


## The Obs class
The Obs class holds the measurements and metadata for one timeseries. There are currently 5 specific Obs classes for different types of measurements:
- GroundwaterObs: for groundwater measurements
- GroundwaterQualityObs: for groundwater quality measurements
- WaterlvlObs: for surface water level measurements
- ModelObs: for hydropandas from a MODFLOW model
- KnmiObs: for (daily) KNMI hydropandas

Each of these Obs classes is essentially a pandas DataFrame with additional methods and attributes related to the type of measurement that it holds. The classes also contain specific methods to read data from specific sources.

## The ObsCollection class
The ObsCollection class, as the name implies, represents a collection of Obs classes, e.g. 10 timeseries of the groundwater level in a certain area. The ObsCollection is also a pandas DataFrame in which each timeseries is stored in a different row. Each row contains metadata (e.g. latitude and longitude of the observation point) and the Obs object (DataFrame) that holds the measurements. It is recommended to let an ObsCollection contain only one Obs type, e.g. to create an ObsCollection for 10 GroundwaterObs, and a separate ObsCollection for 5 KnmiObs.

Like the Obs class, the ObsCollection class contains a bunch of methods for reading data from different sources. See the next section for supported data sources.

## Supported data sources
Currently supported datasources that can be read:
- FEWS PI-XML
- [DINO](www.dinoloket.nl) csv
- WISKI csv
- Artesia Fieldlogger for [Android](https://play.google.com/store/apps/details?id=nl.artesia.fieldlogger&hl=en) and [iOS](https://apps.apple.com/nl/app/fieldlogger/id924565721)
- [Pastas](https://github.com/pastas/pastas) projects (deprecated)
- [Pastastore](https://github.com/pastas/pastastore), for managing Pastas timeseries and models
- [PyStore](https://github.com/ranaroussi/pystore), a fast datastore for pandas timeseries
- [Arctic](https://github.com/man-group/arctic), a timeseries / dataframe database that sits atop MongoDB
- [KNMI](https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script) data
- MODFLOW groundwater models
- IMOD groundwater models

ObsCollection can be exported to:
- Artesia Fieldlogger
- Shapefile
- Pastas projects (deprecated)
- Pastastore
- Arctic
- Pystore

## Example usage
Importing a single DINO csv file:
```python
import observation as obs
fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
gw = obs.GroundwaterObs.from_dino(fname=fname, verbose=True)
```

Or for a zipfile:
```python
import observation as obs
dinozip = './tests/data/2019-Dino-test/dino.zip'
dino_gw = obs.ObsCollection.from_dino(dirname=dinozip,
                                      subdir='Grondwaterstanden_Put',
                                      suffix='1.csv',
                                      ObsClass=obs.GroundwaterObs,
                                      keep_all_obs=False,
                                      verbose=False)
```

## Authors
 - Onno Ebbens, Artesia
 - Ruben Caljé, Artesia
 - Davíd Brakenhoff, Artesia
