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
- `Pastastore <https://github.com/pastas/pastastore>`_, for managing Pastas timeseries and models
- `KNMI data <https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script>`_
- MODFLOW groundwater models
- IMOD groundwater models

ObsCollection can be exported to:

- Shapefile
- Pastastore

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
