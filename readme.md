# observations

The observations module is a Python package for reading timeseries data into DataFrames. The basic idea behind the package is to allow users to manipulate data using all of the wonderful features included in pandas, but to also allow the addition of custom methods and attributes related to the timeseries. The observations module extends pandas.DataFrame with extra functions and metadata related to the type of measurements that are loaded.

## The Obs class
The Obs class holds the measurements and metadata for one timeseries. There are currently 4 specific Obs classes for different types of measurements:
- GroundwaterObs: for groundwater measurements
- GroundwaterQualityObs: for groundwater quality measurements
- WaterlvlObs: for surface water level measurements
- ModelObs: for observations from a MODFLOW model

Each of these Obs classes is essentially a pandas DataFrame with additional methods and attributes related to the type of measurement that it holds. The classes also contain specific methods to read data from specific sources.

## The ObsCollection class
The ObsCollection class, as the name implies, represents a collection of Obs classes, e.g. 10 timeseries of the groundwater level in a certain area. The ObsCollection is also a pandas DataFrame in which each timeseries is stored in a different row. Each row contains metadata (e.g. latitude and longitude of the observation point) and the Obs class (DataFrame) that holds the measurements.

Like the Obs class, the ObsCollection class contains a bunch of methods for reading data from different sources. See the next section for supported data sources.

## Supported data sources
Currently supported datasources that can be read:
- FEWS PI-XML
- [DINO](www.dinoloket.nl) csv
- WISKI csv
- Artesia Fieldlogger for [Android](https://play.google.com/store/apps/details?id=nl.artesia.fieldlogger&hl=en) and [iOS](https://apps.apple.com/nl/app/fieldlogger/id924565721)
- [Pastas](https://github.com/pastas/pastas) projects
- MODFLOW groundwater models
- IMOD

ObsCollection can be exported to:
- Artesia Fieldlogger
- Shapefile
- Pastas projects

## Example usage
Importing a single DINO csv file:
```python
from observations import observation as obs
fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
gw = obs.GroundwaterObs.from_dino_file(fname=fname, verbose=True)
```

Or for a zipfile:
```python
from observations import obs_collection as oc
dinozip = './tests/data/2019-Dino-test/Dino.zip'
dino_gw = oc.ObsCollection.from_dino_dir(dirname=dinozip,
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
