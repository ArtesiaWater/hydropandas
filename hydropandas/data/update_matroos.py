import json

import requests
from pyproj import Proj, Transformer

from hydropandas.io import matroos
from hydropandas.util import EPSG_28992

# %% parse locations, sources and units
URL = "https://noos.matroos.rws.nl/direct/get_series.php?"


def request_api_sources(location, unit):
    """get a list with all available sources at a location for a specific unit

    Parameters
    ----------
    location : str
        location
    unit : str
        unit

    Returns
    -------
    list of str
        sources
    """

    params = {"loc": location, "unit": unit}

    r = requests.get(URL, params=params)
    r.raise_for_status()

    return parse_sources(r.text)


def request_api_units(location):
    """get a list with all available units at a location.

    Parameters
    ----------
    location : str
        location

    Returns
    -------
    list of str
        units
    """
    params = {"loc": location}

    r = requests.get(URL, params=params)
    r.raise_for_status()

    return parse_units(r.text)


def request_api_locations():
    """get a list with all available locations.

    Returns
    -------
    list of str
        locations
    """
    r = requests.get(URL)
    r.raise_for_status()

    return parse_locations(r.text)


def parse_locations(result_str):
    """parse locations from api result

    Parameters
    ----------
    result_str : str
        api output

    Returns
    -------
    list of str
        locations
    """
    lsu = result_str.split("Locations:")[1]
    llocation, _su = lsu.split("Sources:")
    return llocation.split("\n")[2:-2]


def parse_units(result_str):
    """parse units from api result

    Parameters
    ----------
    result_str : str
        api output

    Returns
    -------
    list of str
        units
    """
    _, su = result_str.split("Sources:")
    _s, u = su.split("Units:")

    return u.split("\n")[2:-2]


def parse_sources(result_str):
    """parse sources from api result

    Parameters
    ----------
    result_str : str
        api output

    Returns
    -------
    list of str
        sources
    """
    _, su = result_str.split("Sources:")

    return su.split("\n")[2:-2]


# get all locations
locations = request_api_locations()

# prepare conversion from lat/lon to RD (28992)
proj_from = Proj("EPSG:4326")
proj_to = Proj(EPSG_28992)
transformer = Transformer.from_proj(proj_from, proj_to)

# get possible combinations of location, unit and source (takes approx 30 minutes)
params_dic = {}
for location in locations:
    udic = {}
    units = request_api_units(location)
    for unit in units:
        sources = request_api_sources(location, unit)
        if len(sources) > 0:
            udic[unit] = sources
        else:
            print(f"no sources for {location=} and {unit=}")
    params_dic[location] = {"units": udic}

    # add xy coordinates
    unit = list(params_dic[location]["units"].keys())[0]
    source = params_dic[location]["units"][unit][0]
    _, meta = matroos.get_matroos_obs(location, source, unit, only_metadata=True)
    xy = transformer.transform(meta["lat"], meta["lon"])
    params_dic[location]["coords"] = xy

# write to json
with open("matroos_params.json", "w") as f:
    json.dump(params_dic, f, indent=4)
