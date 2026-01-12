
.. _UserGuide:

==========
User Guide  
==========

This guide covers the core concepts and features of hydropandas. This guide 
will help you make the most of hydropandas.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

Hydropandas is built around two main classes that work together to provide a powerful 
framework for hydrological data analysis:

* **Obs**: Individual time series with metadata
* **ObsCollection**: Collections of multiple time series

Both classes inherit from pandas DataFrame, giving you access to all pandas functionality 
while adding specialized methods for hydrological data.

The ``Obs`` Class
-----------------

The Obs class represents a single time series of hydrological measurements at a specific 
location. It combines your time series data with rich metadata in a single, easy-to-use object.

**Available Obs Types**

Hydropandas provides specialized Obs classes for different measurement types:

- **GroundwaterObs**: Groundwater level and quality measurements
- **WaterlvlObs**: Surface water level measurements  
- **WaterQualityObs**: Water quality parameters (temperature, pH, conductivity, etc.)
- **ModelObs**: Model output from MODFLOW and other groundwater models
- **MeteoObs**: Meteorological observations (general weather data)
- **PrecipitationObs**: Precipitation measurements (subclass of MeteoObs)
- **EvaporationObs**: Evaporation measurements (subclass of MeteoObs)

**Key Features**

Each Obs object includes:

- **Time series data**: Your measurements as a pandas DataFrame
- **Metadata**: Location coordinates, measurement units, data source, etc.
- **Specialized methods**: Plotting, analysis, and export functions

**Example: Working with a GroundwaterObs**

.. code-block:: python

    import hydropandas as hpd
    
    # Create from various sources
    fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
    gw = hpd.GroundwaterObs.from_dino(path=fname)
    gw = hpd.GroundwaterObs.from_bro("GMW000000041261", tube_nr=1, tmin='2020-01-01')
    
    # Access data and metadata
    print(f"Location: {gw.x:.0f}, {gw.y:.0f}")
    print(f"Unit: {gw.unit}")
    print(f"Data range: {gw.index[0]} to {gw.index[-1]}")
    
    # Analyze the data
    gw.plot()  # Quick visualization
    print('statistics:\n',gw.describe())  # Statistical summary
    
    # Export results
    gw.to_csv('output.csv')

The ``ObsCollection`` Class
---------------------------

The ObsCollection class manages collections of Obs objects, enabling analysis across 
multiple measurement points simultaneously. It's particularly powerful for:

* **Spatial analysis**: Compare measurements across different locations
* **Bulk operations**: Apply operations to multiple time series at once
* **Data management**: Organize and export large datasets efficiently

**Best Practices**

* **Consistent Obs Types**: Keep one observation type per collection (e.g., all GroundwaterObs)
* **Logical Grouping**: Group related measurements (same region, project, or time period)  
* **Metadata Management**: Ensure consistent metadata across observations

**Example: Working with ObsCollection**

.. code-block:: python

    import hydropandas as hpd
    from shapely.geometry import Point
    
    # Read multiple observations
    dinozip = './tests/data/2019-Dino-test/dino.zip'
    obs_coll = hpd.read_dino(
        dirname=dinozip,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        ObsClass=hpd.GroundwaterObs,
        keep_all_obs=False
    )
    
    # Spatial operations
    gdf = obs_coll.to_gdf()  # Convert to GeoDataFrame
    print('\ndistance to point:\n',gdf.distance(Point(213000, 473500)))
    
    # Bulk analysis
    print('\nmean value:\n', obs_coll.obs.apply(lambda x: x['stand_m_tov_nap'].mean()))

    # Export to excel file
    obs_coll.to_excel('all_observations.xlsx')

Data Sources and Import Methods
------------------------------

Hydropandas supports numerous data sources with specialized import methods:

**Dutch Data Sources**
* **BRO (Basisregistratie Ondergrond)**: ``read_bro()``
* **DINO Loket**: ``read_dino()``  
* **KNMI**: ``read_knmi()``
* **Waterinfo (RWS)**: ``read_waterinfo()``
* **Lizard**: ``read_lizard()``

**International Sources**  
* **Water Connect (Australia)**: ``read_waterconnect()``

**File Formats**
* **CSV files**: ``from_csv()``
* **Excel files**: ``read_excel()``
* **JSON files**: ``from_json()``
* **FEWS PI-XML**: ``read_fews()``
* **WISKI exports**: ``read_wiski()``

**Models**
* **MODFLOW**: ``read_modflow()``
* **IMOD**: ``read_imod()``

**Example: Reading Different Sources**

.. code-block:: python

    extent = [110000, 120000, 418000, 420000]

    # Government monitoring network
    bro_obs = hpd.read_bro(extent=extent)
    
    # Surface water
    rws_data = hpd.read_matroos(extent=extent)

    # Weather data  
    knmi_data = hpd.read_knmi(locations=rws_data, meteo_vars=['RH'])


Advanced Features
-----------------

**Spatial Analysis**

Hydropandas integrates seamlessly with GeoPandas for spatial analysis:

.. code-block:: python

    # Convert to GeoDataFrame for spatial operations
    gdf = obs_collection.to_gdf()
    
    # Spatial filtering
    locations = gdf.clip_by_rect(110000, 400000, 115000, 420000)
    print('\nlocations in rectangle:\n',locations.index)

**Temporal Analysis**

Built-in methods for temporal analysis:

.. code-block:: python

    # Resample to different frequencies
    monthly = obs.resample('ME').mean()
    annual = obs.resample('YE').mean()

**Visualization**

Multiple plotting options for data exploration:

.. code-block:: python

    # Single observation plots
    obs.plot()
    obs.plot(kind='hist')
    obs.plot(kind='box')
    
    # Collection plots
    obs_coll.plots.interactive_map()

Export and Integration
----------------------

**Export Formats**

.. code-block:: python

    # Spreadsheet formats
    obs_coll.to_excel('data.xlsx')  # Multi-sheet Excel
    obs.to_csv('single_series.csv')
    
    # GIS formats
    obs_coll.to_gdf().to_file('data.gpkg')  # GeoPackage
    obs_coll.to_gdf().to_file('data.shp')   # Shapefile
    
    # Archive formats
    obs_coll.to_pickle('data.pkl')  # Python pickle
    obs_coll.to_json('data.json')   # JSON format

**Integration with Other Tools**

.. code-block:: python

    # Pastas (time series modeling)
    import pastas as ps
    ml = ps.Model(obs)
    
    # xarray for multi-dimensional analysis
    xr_data = obs_coll.to_xarray()
    
    # Direct pandas access
    df = obs[['value']]  # Standard pandas operations

Configuration and Settings
--------------------------

**Global Settings**

.. code-block:: python

    import hydropandas as hpd
    
    # Configure logging level
    hpd.util.get_color_logger('DEBUG')


**Getting Help**

* Check the :doc:`examples gallery <examples/index>` for similar use cases
* Search existing `GitHub Issues <https://github.com/ArtesiaWater/hydropandas/issues>`_
* Ask questions in `GitHub Discussions <https://github.com/ArtesiaWater/hydropandas/discussions>`_
* Contribute improvements via :doc:`contributing guide <contribute>`