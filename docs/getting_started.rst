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

Dependencies
------------
This module has several optional dependencies that have to be installed. 
These include:

- pastastore (create pastas models of an ObsCollection)
- folium and bokeh (make an interactive map of an ObsCollection)
- xarray and netCDF4 (get regis layers for groundwater observations)
- flopy (interaction with modflow data)

See the :ref:`examples` section for some quick examples on how to get started.

