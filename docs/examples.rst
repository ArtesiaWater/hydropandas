========
Examples
========
This page provides some short example code snippets and links to 
Jupyter Notebooks with more detailed examples.

Example snippets
----------------

Importing a single CSV-file downloaded from DINOLoket::

   import observation as obs
   fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
   gw = obs.GroundwaterObs.from_dino(fname=fname, verbose=True)


Or for a zipfile::

   import observation as obs
   dinozip = './tests/data/2019-Dino-test/dino.zip'
   dino_gw = obs.ObsCollection.from_dino(dirname=dinozip,
                                         subdir='Grondwaterstanden_Put',
                                         suffix='1.csv',
                                         ObsClass=obs.GroundwaterObs,
                                         keep_all_obs=False,
                                         verbose=False)

Example notebooks
-----------------

The links below link to Jupyter Notebooks with explanation and examples of the
usage of the `hydropandas` module:

.. toctree::
   :maxdepth: 2
   :glob:

   examples/*
