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

- `BRO <examples/01_groundwater_observations.html>`_
- `DINOLoket CSV <examples/01_groundwater_observations.html>`_
- `FEWS PI-XML <examples/07_fews.html>`_
- IMOD groundwater models
- `KNMI data <examples/02_knmi_observations.html>`_
- `Matroos data <examples/11_matroos.html>`_
- `Lizard <examples/06_lizard.html>`_
- MODFLOW groundwater models
- `Pastastore <examples/03_hydropandas_and_pastas.html>`_, for managing Pastas timeseries and models
- `Water connect <examples/09_waterconnect.html>`_
- `Waterinfo <examples/08_waterinfo.html>`_
- WISKI csv files

ObsCollection can be exported to:
- Excel (with one tab for each time series and a single tab with all metadata)
- Pastastore
- Pickle
- Shapefile/geopackage/geojson

See the table of contents to get started with `hydropandas`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Getting started <getting_started>
   Examples gallery <examples/index>
   User guide <user_guide>
   Hydropandas API-docs <source/modules>
   Contribute <contribute>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
