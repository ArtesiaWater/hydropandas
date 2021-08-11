===============
Getting Started
===============

On this page you will find all the information to get started with `hydropandas`.

Getting Python
--------------
To install `hydropandas`, a working version of Python 3.7 or higher has to be
installed on your computer. We recommend using the
`Anaconda Distribution <https://www.continuum.io/downloads>`_
of Python.

Installing hydropandas
----------------------

Install the module by typing:: 

    pip install hydropandas
	
Please note that some of the dependencies cannot be installed automatically 
on Windows. If you do not have these packages already you can install them 
manually using the following instructions:

Download these packages from `Christoph Gohlke's website <https://www.lfd.uci.edu/~gohlke/pythonlibs>`_ :

- GDAL
- Fiona
- Shapely
- Python-snappy
- Fastparquet

Use CTRL+F to find the download link on the page. Be sure to download the 
correct version of the package. The Python version should match your Python 
version. Also the architecture should match (i.e. 64bits vs 32bits). For example:

- GDAL-3.1.4-cp38-cp38-win_amd64.whl

This is the GDAL version for Python `3.8` (as can be seen from the `cp38` in the name), 
for `64-bits` Python (as derived from the `amd64` in the name).

Once you have downloaded the correct files, open Anaconda Prompt, and navigate to 
the directory in which you saved your downloads. Now use the following commands 
to install the packages (the order is important)::

	pip install GDAL-3.1.4-cp38-cp38-win_amd64.whl
	pip install Fiona-1.8.17-cp38-cp38-win_amd64.whl
	pip install Shapely-1.7.1-cp38-cp38-win_amd64.whl
	pip install python_snappy-0.5.4-cp38-cp38-win_amd64.whl
	pip install fastparquet-0.4.1-cp38-cp38-win_amd64.whl

After you've done this you can install hydropandas using::

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

