Welcome to hydropandas's documentation!
=======================================
Hydropandas is a Python package for reading, analyzing and writing
(hydrological) timeseries data. Users can store a timeseries and metadata in a
single object. This object inherits from a pandas DataFrame, with all its
wonderful features, and is extended with custom methods and attributes related
to hydrological timeseries.

Supported data sources
----------------------

Currently supported datasources that can be read:

- `BRO <https://www.broloket.nl>`_
- `DINOLoket CSV <www.dinoloket.nl>`_
- `FEWS PI-XML <https://publicwiki.deltares.nl/display/FEWSDOC/The+Delft-Fews+Published+Interface>`_
- IMOD groundwater models
- `KNMI data <https://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script>`_
- `Lizard <https://vitens.lizard.net>`_
- MODFLOW groundwater models
- `Pastastore <https://github.com/pastas/pastastore>`_, for managing Pastas timeseries and models
- `Waterinfo <https://waterinfo.rws.nl>`_
- WISKI csv files

ObsCollection can be exported to:

- Excel (with one tab for each time series and a single tab with all metadata)
- Shapefile
- Pastastore

See the table of contents to get started with `hydropandas`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Getting started <getting_started>
   Examples <examples/index>
   User guide <user_guide>
   Hydropandas API-docs <source/modules>
   Contribute <contribute>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
