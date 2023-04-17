
.. _UserGuide:

==========
User guide
==========

The `Obs` class
---------------

The Obs class holds the measurements and metadata for one timeseries. There 
are currently 5 specific Obs classes for different types of measurements:

- GroundwaterObs: for groundwater measurements
- WaterQualityObs: for groundwater quality measurements
- WaterlvlObs: for surface water level measurements
- ModelObs: for "observations" from a MODFLOW model
- MeteoObs: for meteorological observations
- PrecipitationObs: for precipitation observations, subclass of MeteoObs
- EvaporationObs: for evaporation observations, subclass of MeteoObs

Each of these Obs classes is essentially a pandas DataFrame with additional 
methods and attributes related to the type of measurement that it holds. The 
classes also contain specific methods to read data from specific sources.

The `ObsCollection` class
-------------------------

The ObsCollection class, as the name implies, represents a collection of Obs 
classes, e.g. 10 timeseries of the groundwater level in a certain area. The 
ObsCollection is also a pandas DataFrame in which each timeseries is stored 
in a different row. Each row contains metadata (e.g. latitude and longitude 
of the observation point) and the Obs object that holds the 
measurements. It is recommended to let an ObsCollection contain only one Obs 
type, e.g. to create an ObsCollection for 10 GroundwaterObs, and a separate 
ObsCollection for 5 KnmiObs.

Like the Obs class, the ObsCollection class contains a bunch of methods for 
reading data from different sources.