========
Examples
========
This page provides some short example code snippets and links to 
Jupyter Notebooks with more detailed examples.

Example snippets
----------------

Reading a single CSV-file downloaded from DINOLoket::

   import hydropandas as hpd
   fname = './tests/data/2019-Dino-test/Grondwaterstanden_Put/B33F0080001_1.csv'
   gw = hpd.GroundwaterObs.from_dino(path=fname)


Or for a zipfile::

   import hydropandas as hpd
   dinozip = './tests/data/2019-Dino-test/dino.zip'
   dino_gw = hpd.ObsCollection.from_dino(dirname=dinozip,
                                         subdir='Grondwaterstanden_Put',
                                         suffix='1.csv',
                                         ObsClass=hpd.GroundwaterObs,
                                         keep_all_obs=False)

Example gallery
---------------

.. nbgallery::
    :name: examples

    ../examples/00_hydropandas_objects.ipynb
    ../examples/01_groundwater_observations.ipynb
    ../examples/02_knmi_observations.ipynb
    ../examples/03_hydropandas_and_pastas.ipynb
    ../examples/04_merging_observations.ipynb
    ../examples/05_bronhouderportaal_bro.ipynb
    ../examples/06_lizard.ipynb
    ../examples/07_fews.ipynb
    ../examples/08_waterinfo.ipynb
    ../examples/09_water_connect.ipynb

