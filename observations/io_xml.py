import numpy as np
import pandas as pd
import xml.etree.ElementTree as etree


def read_xml(fname, ObsClass, to_mnap=False, remove_nan=False, verbose=False):
    """read a FEWS XML-file with measurements, return list of ObsClass objects

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    to_mnap : boolean, optional
        if True a column with 'Stand_m_tov_NAP' is added to the dataframe
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects

    """
    tree = etree.parse(fname)
    root = tree.getroot()
    obs_list = []
    for i in range(len(root)):
        if root[i].tag.endswith('series'):
            series = {}
            date = []
            time = []
            events = []
            for j in range(len(root[i])):
                if root[i][j].tag.endswith('header'):
                    for k in range(len(root[i][j])):
                        prop = root[i][j][k].tag.split('}')[-1]
                        val = root[i][j][k].text
                        if prop == 'x' or prop == 'y' or prop == 'lat' or prop == 'lon':
                            val = float(val)
                        series[prop] = val
                        if verbose:
                            if prop == 'locationId':
                                print('read {}'.format(val))
                elif root[i][j].tag.endswith('event'):
                    date.append(root[i][j].attrib['date'])
                    # flag.append(root[i][j].attrib['flag'])
                    time.append(root[i][j].attrib['time'])
                    # value.append(root[i][j].attrib['value'])
                    events.append({**root[i][j].attrib})
            # combine events in a dataframe
            index = pd.to_datetime(
                [d + ' ' + t for d, t in zip(date, time)])
            ts = pd.DataFrame(events, index=index, dtype=float)
            if remove_nan:
                ts.dropna(subset=['value'], inplace=True)
            if to_mnap:
                ts['Stand_m_tov_NAP'] = ts['value']

            if "x" in series.keys():
                x = series["x"]
            else:
                x = np.nan
            if "y" in series.keys():
                y = series["x"]
            else:
                y = np.nan

            obs_list.append(ObsClass(ts, name=series['locationId'],
                                     x=x, y=y, meta=series))
    return obs_list


def read_xml_alternative(fname, ObsClass, to_mnap=False, remove_nan=False, verbose=False):
    """read a FEWS XML-file with measurements, return list of ObsClass objects

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    to_mnap : boolean, optional
        if True a column with 'Stand_m_tov_NAP' is added to the dataframe
    remove_nan : boolean, optional
        remove nan values from measurements, flag information about the
        nan values is also lost
    verbose : boolean, optional
        print additional information to the screen (default is False).

    Returns
    -------
    list of ObsClass objects
        list of timeseries stored in ObsClass objects

    """

    from lxml import etree
    context = etree.iterparse(fname, events=('start',),
                              tag=("{http://www.wldelft.nl/fews/PI}header",
                                   "{http://www.wldelft.nl/fews/PI}event"))

    header_list = []
    events = []

    for event, element in context:
        header = {}
        if element.tag.endswith("header") and event == "start":
            for h_attr in element:
                tag = h_attr.tag.replace("{{{0}}}".format(
                    "http://www.wldelft.nl/fews/PI"), "")
                if h_attr.text is not None:
                    header[tag] = h_attr.text
                elif len(h_attr.attrib) != 0:
                    header[tag] = {**h_attr.attrib}
                else:
                    header[tag] = None
                if verbose:
                    if tag.startswith("locationId"):
                        print("read {}".format(header[tag]))
            label = header["locationId"] + "|" + header["parameterId"]
            header_list.append(header)
        elif element.tag.endswith("event") and event == "start":
            events.append({**element.attrib, "loc_param": label})

        # Free memory.
        element.clear()
        while element.getprevious() is not None:
            del element.getparent()[0]

    # combine events in a dataframe
    ts = pd.DataFrame(events)
    ts["datetime"] = pd.to_datetime(ts["date"] + " " + ts["time"])
    ts.set_index("datetime", inplace=True)
    ts.drop(["date", "time"], axis=1, inplace=True)
    ts["value"] = pd.to_numeric(ts["value"], errors="coerce")
    ts["flag"] = pd.to_numeric(ts["flag"], errors="coerce")

    if remove_nan:
        ts.dropna(subset=['value'], inplace=True)
    if to_mnap:
        ts['Stand_m_tov_NAP'] = ts['value']

    obs_list = []
    for i, (name, _) in enumerate(ts.groupby(by="loc_param")):
        if "x" in header_list[i].keys():
            x = header_list[i]["x"]
        else:
            x = np.nan
        if "y" in header_list[i].keys():
            y = header_list[i]["x"]
        else:
            y = np.nan

        obs_list.append(ObsClass(ts, name=header_list[i]['locationId'],
                                 x=x, y=y, meta=header_list[i]))

    return obs_list
