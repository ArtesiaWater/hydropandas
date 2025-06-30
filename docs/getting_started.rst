===============
Getting Started
===============

On this page you will find all the information to get started with `hydropandas`.

Getting Python
--------------
To install `hydropandas`, a working version of Python 3.9 or higher has to be
installed on your computer. We recommend using the
`Anaconda Distribution <https://www.continuum.io/downloads>`_
of Python.

Installing hydropandas
----------------------

Install the module by typing:: 

    pip install hydropandas

For installing in development mode, clone the repository and install by
typing the following from the module root directory::

    pip install -e .

Using `hydropandas`
-------------------

Start Python and import the module::

    import hydropandas as hpd


Read a single CSV-file downloaded from DINOLoket and plot the measurements::

   import hydropandas as hpd
   fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
   gw = hpd.GroundwaterObs.from_dino(path=fname)
   gw['stand_m_tov_nap'].plot()

Or read a zipfile and plot the location of the measurements on a map::

   import hydropandas as hpd
   import contextily as ctx
   dinozip = './tests/data/2019-Dino-test/dino.zip'
   dino_gw = hpd.ObsCollection.from_dino(dirname=dinozip,
                                           subdir='Grondwaterstanden_Put',
                                           suffix='1.csv',
                                           ObsClass=hpd.GroundwaterObs,
                                           keep_all_obs=False)
   ax = dino_gw.to_gdf().plot()
   ctx.add_basemap(ax=ax, crs=28992)

For more examples please see the `Examples gallery <examples/index>`_.