Welcome to hydropandas's documentation!
=======================================

Hydropandas is a Python package for reading, analyzing and writing
(hydrological) timeseries data. Users can store a timeseries and metadata in a
single object. This object inherits from a pandas DataFrame, with all its
wonderful features, and is extended with custom methods and attributes related
to hydrological timeseries.

Supported data sources
----------------------

**Dutch Sources**

* **BRO** (Basisregistratie Ondergrond) - `Tutorial <examples/01_groundwater_observations.html>`_
* **DINO Loket** CSV files - `Tutorial <examples/01_groundwater_observations.html>`_
* **KNMI** weather data - `Tutorial <examples/02_knmi_observations.html>`_  
* **Waterinfo** (RWS) - `Tutorial <examples/08_waterinfo.html>`_
* **Lizard** platform - `Tutorial <examples/06_lizard.html>`_
* **Matroos** system - `Tutorial <examples/11_matroos.html>`_

**International Sources**

* **Water Connect** (Australia) - `Tutorial <examples/09_waterconnect.html>`_

**File Formats**

* **FEWS PI-XML** - `Tutorial <examples/07_fews.html>`_
* **WISKI** CSV exports
* **Excel/CSV** files
* **MODFLOW/IMOD** model outputs

**Integration Platforms**

* **Pastastore** - `Tutorial <examples/03_hydropandas_and_pastas.html>`_

Export Capabilities
-------------------

**Export Formats**

* **Excel** (multi-sheet with metadata)
* **GeoPackage/Shapefile** (with spatial data)
* **JSON/CSV** (for data exchange)
* **Pickle** (for Python workflows)

**Integration**

* **geopandas** (full DataFrame compatibility)
* **Pastas** (time series modeling)
* **Pastastore** (bulk time series modeling)

See the table of contents to get started with hydropandas.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Getting started <getting_started>
   Examples gallery <examples/index>
   User guide <user_guide>
   FAQ <faq>
   Hydropandas API-docs <source/modules>
   Contribute <contribute>

.. toctree::
   :maxdepth: 1 
   :caption: Quick Links:
   
   Installation Guide <getting_started>
   Quick Start Tutorial <examples/00_hydropandas_objects>
   Interactive Examples <examples/index>
   FAQ & Troubleshooting <faq>
   API Reference <source/modules>
   Contributing <contribute>
   Report Issues <https://github.com/ArtesiaWater/hydropandas/issues>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
