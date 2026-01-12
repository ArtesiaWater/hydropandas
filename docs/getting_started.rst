===============
Getting Started
===============

Welcome to hydropandas! This guide will help you get up and running quickly with hydrological time series analysis in Python.

What is hydropandas?
--------------------

Hydropandas is a Python package designed specifically for working with hydrological time series data. It extends pandas DataFrame functionality with specialized methods for:

* Reading data from multiple hydrological data sources
* Managing metadata alongside time series data
* Visualizing temporal and spatial patterns
* Exporting data to various formats

Key Features
------------

* **Multiple Data Sources**: BRO, KNMI, Lizard, Waterinfo, FEWS, and more
* **Rich Metadata**: Store location, units, and source information with your data
* **Spatial Analysis**: Built-in GIS capabilities with geopandas integration
* **Interactive Visualizations**: Maps and plots for data exploration
* **Export Options**: Excel, csv, GeoPackage, Shapefile, JSON, and more

Installation
------------

**Basic Installation**

Install the package::

    pip install hydropandas

**Full Installation (Recommended)**

For complete functionality, install with all optional dependencies::

    pip install hydropandas[full]

**Development Installation**

If you want to contribute or use the latest features::

    git clone https://github.com/ArtesiaWater/hydropandas.git
    cd hydropandas
    pip install -e .[full]

Quick Start
-----------

Import hydropandas and start exploring::

    import hydropandas as hpd
    
    # Enable informative logging
    hpd.util.get_color_logger("INFO")

Common Use Cases
----------------

**1. Single Time Series Analysis**

Read a groundwater observation from DINO and visualize it::

   import hydropandas as hpd

   # Read single groundwater observation from a Dino csv file
   fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
   gw = hpd.GroundwaterObs.from_dino(path=fname)
   
   # Quick plot
   ax = gw['stand_m_tov_nap'].plot(title=gw.name)
   ax.set_ylabel(gw.unit)

**2. Multiple Observations**

Work with multiple observation points simultaneously::

   import hydropandas as hpd
   import contextily as ctx
   
   # Read collection of observations
   dinozip = './tests/data/2019-Dino-test/dino.zip'
   obs_collection = hpd.read_dino(
       dirname=dinozip,
       subdir='Grondwaterstanden_Put',
       suffix='1.csv',
       ObsClass=hpd.GroundwaterObs,
       keep_all_obs=False
   )
   
   # Create spatial map
   ax = obs_collection.to_gdf().plot('location',figsize=(10, 8), legend=True)
   ctx.add_basemap(ax=ax, crs=28992)
   ax.set_title('Observation Locations')

**3. API Data Access**

Fetch real-time data from web services::

   # Get KNMI weather data
   knmi_obs = hpd.read_knmi(
       stns=[260], # De Bilt
       meteo_vars=['RH'],  # Precipitation
       starts='2020-01-01',
       ends='2020-12-31'
   )
   
   # Get groundwater data from BRO
   bro_obs = hpd.read_bro(
       extent=[110000, 120000, 418000, 420000]  # Extent in RD coordinates
       tmin='2020-01-01',
       tmax='2020-12-31'
   )

Next Steps
----------

Ready to dive deeper? Check out:

* :doc:`Examples gallery <examples/index>` - Comprehensive tutorials and use cases
* :doc:`User guide <user_guide>` - Detailed explanation of core concepts  
* :doc:`API documentation <source/modules>` - Complete reference for all functions
* :doc:`Contributing guide <contribute>` - Help improve hydropandas

Need Help?
----------

* **Found a bug?** Report it on `GitHub Issues <https://github.com/ArtesiaWater/hydropandas/issues>`_
* **Have a question?** Start a `GitHub Discussion <https://github.com/ArtesiaWater/hydropandas/discussions>`_  
* **Want to contribute?** See our :doc:`contribute` guide

For more examples please see the :doc:`Examples gallery <examples/index>`.