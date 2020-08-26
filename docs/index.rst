Welcome to hydropandas's documentation!
=======================================
The hydropandas module is a Python package for reading timeseries data into 
DataFrames. The basic idea behind the package is to allow users to manipulate 
data using all of the wonderful features included in pandas, but to also allow 
the addition of custom methods and attributes related to the timeseries. The 
hydropandas module extends pandas.DataFrame with extra functionality and 
stores metadata related to the type of measurements.

Supported data sources
----------------------

Currently supported datasources that can be read:

- `FEWS PI-XML <https://publicwiki.deltares.nl/display/FEWSDOC/The+Delft-Fews+Published+Interface>`_
- `DINOLoket CSV <www.dinoloket.nl>`_
- WISKI CSV
- Artesia Fieldlogger for 
  `Android <https://play.google.com/store/apps/details?id=nl.artesia.fieldlogger&hl=en>`_ and 
  `iOS <https://apps.apple.com/nl/app/fieldlogger/id924565721>`_
- `Pastas projects <https://github.com/pastas/pastas>`_ (deprecated, see pastastore)
- `Pastastore <https://github.com/pastas/pastastore>`_, for managing Pastas timeseries and models
- `PyStore <https://github.com/ranaroussi/pystore>`_, a fast datastore for pandas timeseries
- `Arctic <https://github.com/man-group/arctic>`_, a timeseries / dataframe database that sits atop MongoDB
- `KNMI data <https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script>`_
- MODFLOW groundwater models
- IMOD groundwater models

ObsCollection can be exported to:

- Artesia Fieldlogger
- Shapefile
- Pastas projects (deprecated)
- Pastastore
- Arctic
- Pystore

See the table of contents to get started with `hydropandas`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Getting started <getting_started>
   Example usage <examples>
   User guide <user_guide>
   hydropandas API-docs <source/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
