===============================
Frequently Asked Questions (FAQ)
===============================

This page answers common questions about using hydropandas effectively.

Installation & Setup
--------------------

**Q: Should I install hydropandas or hydropandas[full]?**

A: For most users, we recommend ``pip install hydropandas[full]`` to get all optional dependencies. 
This ensures all data sources and visualization features work properly.

**Q: I get import errors after installation. What's wrong?**

A: This usually means you're missing optional dependencies. Try:

.. code-block:: bash

    pip install hydropandas[full]

If problems persist, check our `installation troubleshooting guide <getting_started.html#installation>`_.

Data Sources & APIs
-------------------

**Q: Which data sources work without API keys?**

A: These sources work immediately:

* BRO (Dutch groundwater registry)
* KNMI (Dutch weather service) 
* Waterinfo (Dutch surface water)
* CSV/Excel files

**Q: How do I get API access for Lizard?**

A: Contact your Lizard organization administrator for API credentials. Then:

.. code-block:: python

    auth = ('username', 'password')  # or API key
    obs = hpd.read_lizard(extent=extent, auth=auth)

**Q: API calls are very slow. How can I speed them up?**

A: Try these optimization strategies:

.. code-block:: python

    # Use smaller time ranges
    obs = hpd.read_bro(extent=extent, tmin='2020-01-01', tmax='2020-12-31')

Data Analysis
-------------

**Q: Can I resample time series to different frequencies?**

A: Yes, use pandas resampling methods:

.. code-block:: python

    # Monthly averages
    monthly = obs.resample('ME').mean()
    
    # Daily maximum values
    daily_max = obs.resample('D').max()
    
    # Annual statistics
    annual_stats = obs.resample('YE').agg(['mean', 'min', 'max', 'std'])


Visualization
-------------

**Q: How do I create interactive maps?**

A: Use the built-in plotting methods:

.. code-block:: python

    # Basic interactive map
    obs_collection.plots.interactive_map()
    
    # Customized map with popup info
    obs_collection.plots.interactive_map(
        popup_width=400,
        tiles='OpenStreetMap'
    )
    
    # Static map with basemap
    import contextily as ctx
    ax = obs_collection.to_gdf().plot()
    ctx.add_basemap(ax=ax, crs=28992)

**Q: Can I customize the plotting style?**

A: Yes, hydropandas integrates with matplotlib and supports customization:

.. code-block:: python

    # Customize plot appearance
    obs.plot(figsize=(12, 6), color='blue', linewidth=2)

Export & Integration  
--------------------

**Q: How do I export data to Excel with proper formatting?**

A: Use the enhanced Excel export:

.. code-block:: python

    # Basic export
    obs_collection.to_excel('data.xlsx')

**Q: Can I integrate hydropandas with other Python packages?**

A: Absolutely! Hydropandas works seamlessly with:

.. code-block:: python

    # Pastas (time series modeling)
    import pastas as ps
    ml = ps.Model(obs['value'])
    
    # Scikit-learn (machine learning)
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(obs[['value']])
    
    # Xarray (multi-dimensional arrays)
    xr_data = obs_collection.to_xarray()

Error Handling
--------------

**Q: I get "No data found" errors. What should I check?**

A: Common causes and solutions:

1. **Check extent coordinates**: Ensure they're in the correct CRS
2. **Verify time range**: Some sources have limited historical data
3. **Check internet connection**: API calls require stable internet
4. **Validate parameters**: Ensure location names and codes are correct

.. code-block:: python

    # Debug mode for detailed error messages
    hpd.util.get_color_logger('DEBUG')
    
    # Check data availability first
    locations = hpd.read_bro(extent=extent, only_metadata=True)
    print(f"Found {len(locations)} locations")

Still Need Help?
----------------

If your question isn't answered here:

1. **Check the examples**: Browse our `examples gallery <examples/index.html>`_
2. **Search existing issues**: Look through `GitHub Issues <https://github.com/ArtesiaWater/hydropandas/issues>`_
3. **Ask the community**: Start a `GitHub Discussion <https://github.com/ArtesiaWater/hydropandas/discussions>`_
4. **Report bugs**: Create a new `GitHub Issue <https://github.com/ArtesiaWater/hydropandas/issues/new>`_

When asking for help, please include:

* Your hydropandas version: ``hpd.show_versions()``
* Complete error messages
* Minimal code example that reproduces the issue
* Your operating system and Python version