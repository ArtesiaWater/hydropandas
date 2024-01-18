========
Examples
========
This page provides some short example code snippets and links to 
Jupyter Notebooks with more detailed examples.

Example snippets
----------------

Importing a groundwater time series from the BRO using the BRO-id and the tube number::

   import hydropandas as hpd
   gw_bro = hpd.GroundwaterObs.from_bro("GMW000000041261", 1)


Or import all groundwater time series from the BRO within a certain extent::

   oc = hpd.read_bro(extent=(117850, 118180, 439550, 439900))

Example notebooks
-----------------

The links below link to Jupyter Notebooks with explanation and examples of the
usage of the `hydropandas` module:

.. toctree::
   :maxdepth: 2
   :glob:

   examples/*
