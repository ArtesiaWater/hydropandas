import numpy as np
import pandas as pd
import xml.etree.ElementTree as etree
from lxml.etree import iterparse, parse


def read_xml(fname, ObsClass, translate_dic={'locationId': 'locatie'}, 
             to_mnap=False, remove_nan=False, verbose=False):
    """read a FEWS XML-file with measurements, return list of ObsClass objects

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
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
                    date.append(root[i][j].attrib.pop('date'))
                    time.append(root[i][j].attrib.pop('time'))
                    events.append({**root[i][j].attrib})
            # combine events in a dataframe
            index = pd.to_datetime(
                [d + ' ' + t for d, t in zip(date, time)],
                errors="coerce")
            ts = pd.DataFrame(events, index=index, dtype=float)
            if remove_nan:
                ts.dropna(subset=['value'], inplace=True)
            if to_mnap:
                ts['stand_m_tov_nap'] = ts['value']

            if "x" in series.keys():
                x = series["x"]
            else:
                x = np.nan
            if "y" in series.keys():
                y = series["y"]
            else:
                y = np.nan
                
            for key, item in translate_dic.items():
                series[item] = series.pop(key)

            obs_list.append(ObsClass(ts, name=series['locatie'],
                                     locatie=series['locatie'],
                                     x=x, y=y, meta=series))
    return obs_list


def iterparse_pi_xml(fname, ObsClass, translate_dic={'locationId': 'locatie'}, 
                     locationIds=None, return_events=True,
                     keep_flags=(0, 1), return_df=False,
                     tags=('series', 'header', 'event'),
                     skip_errors=True, verbose=False):
    """read a FEWS XML-file with measurements, memory efficient

    Parameters
    ----------
    fname : str
        full path to file
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    locations : tuple or list of str, optional
        list of locationId's to read from XML file, others are skipped.
        If None (default) all locations are read.
    return_events : bool, optional
        return all event-information in a DataFrame per location, instead of
        just a Series (defaults to False). Overrules keep_flags kwarg.
    keep_flags : list of ints, optional
        keep the values with these flags (defaults to 0 and 1). Only used
        when return_events is False.
    tags : list of strings, optional
        Select the tags to be parsed. Defaults to series, header and event
    return_df : bool, optional
        return a DataFame with the data, instead of two lists (defaults to
        False)
    skip_errors: bool, optional
        if True, continue after error, else raise error

    Returns
    -------

    if return_df == True:
        df : pandas DataFrame
            a DataFrame containing the metadata and the series
    else:
        header_list : list of dictionaries
            list of metadata
        series_list : list of pandas Series
            list of timeseries

    """

    tags = ['{{http://www.wldelft.nl/fews/PI}}{}'.format(tag) for tag in tags]

    context = iterparse(fname, tag=tags)

    header_list = []
    series_list = []

    keep_flags = [str(flag) for flag in keep_flags]

    try:
        for _, element in context:
            if element.tag.endswith("header"):
                header = {}
                for h_attr in element:
                    tag = h_attr.tag.replace("{{{0}}}".format(
                        "http://www.wldelft.nl/fews/PI"), "")
                    # if specific locations are provided only read those
                    if locationIds is not None and tag.startswith("locationId"):
                        loc = h_attr.text
                        if loc not in locationIds:
                            continue
                    if h_attr.text is not None:
                        header[tag] = h_attr.text
                    elif len(h_attr.attrib) != 0:
                        header[tag] = {**h_attr.attrib}
                    else:
                        header[tag] = None
                    if verbose and tag.startswith("locationId"):
                        print("reading {}".format(header[tag]))
                events = []
            elif element.tag.endswith("event"):
                # if specific locations are provided only read those
                if locationIds is not None:
                    if loc not in locationIds:
                        continue
                events.append({**element.attrib})
            elif element.tag.endswith('series'):
                # if specific locations are provided only read those
                if locationIds is not None:
                    if loc not in locationIds:
                        continue
                if len(events) == 0:
                    if return_events:
                        s = pd.DataFrame()
                    else:
                        s = pd.Series()
                else:
                    df = pd.DataFrame(events)
                    df.index = pd.to_datetime(
                        [d + " " + t for d, t in zip(df['date'], df['time'])],
                        errors="coerce")
                    df.drop(columns=["date", "time"], inplace=True)
                    if return_events:
                        df['value'] = pd.to_numeric(df['value'], errors="coerce")
                        df['flag'] = pd.to_numeric(df['flag'])
                        s = df
                    else:
                        mask = df['flag'].isin(keep_flags)
                        s = pd.to_numeric(df.loc[mask, 'value'], errors="coerce")

                for key, item in translate_dic.items():
                    header[item] = header.pop(key)

                o = ObsClass(s, name=header['locatie'], 
                             locatie=header['locatie'],
                             meta=header)
                header_list.append(header)
                series_list.append(o)

            # Free memory.
            element.clear()
            # while element.getprevious() is not None:
            #    del element.getparent()[0]

    except Exception as e:
        if skip_errors:
            print("ERROR! Skipped {}".format(fname))
        else:
            raise(e)

    if return_df:
        for h, s in zip(header_list, series_list):
            h['series'] = s
        return pd.DataFrame(header_list)
    else:
        return header_list, series_list


def write_pi_xml(obs_coll, fname, timezone=1.0, version="1.24"):
    """
    Write PiTimeSeries object to PI-XML file.

    Parameters
    ----------
    fname: path
        path to XML file to be written

    """

    assert fname.endswith(
        ".xml"), "Output file should have '.xml' extension!"

    # first line of XML file
    line0 = '<?xml version="1.0" encoding="UTF-8"?>\n'

    # some definitions for timeseries XML file
    NS = r"http://www.wldelft.nl/fews/PI"
    FS = r"http://www.wldelft.nl/fews/fs"
    XSI = r"http://www.w3.org/2001/XMLSchema-instance"
    schemaLocation = (
        r"http://fews.wldelft.nl/schemas/version1.0"
        r"/Pi-schemas/pi_timeseries.xsd"
    )
    timeseriesline = ('<TimeSeries xmlns="{NS}" xmlns:xsi="{XSI}" '
                      'xsi:schemaLocation="{NS} {schema}" version="{version}" '
                      'xmlns:fs="{FS}">\n')

    # line templates
    paramline = "<{tag}>{param}</{tag}>\n"

    # write file
    with open(fname, "w") as f:
        f.write(line0)
        f.write(
            timeseriesline.format(
                NS=NS, FS=FS, XSI=XSI, schema=schemaLocation, version=version
            )
        )
        tzline = "\t" + \
            paramline.format(tag="timeZone", param=timezone)
        f.write(tzline)

        for o in obs_coll.obs:
            # start series
            start = "\t" + "<series>\n"
            f.write(start)
            # write header
            hlines = []
            hstart = 2 * "\t" + "<header>\n"
            hlines.append(hstart)
            for htag, hval in o.meta.items():
                if htag.endswith("Date"):
                    try:
                        hdate = hval.strftime("%Y-%m-%d")
                        htime = hval.strftime("%H:%M:%S")
                    except AttributeError as e:
                        if htag.startswith("start"):
                            hdate = o.index[0].strftime("%Y-%m-%d")
                            htime = o.index[0].strftime("%H:%M:%S")
                        elif htag.startswith("end"):
                            hdate = o.index[-1].strftime("%Y-%m-%d")
                            htime = o.index[-1].strftime("%H:%M:%S")
                        else:
                            raise(e)
                    hline = '<{tag} date="{date}" time="{time}"/>\n'.format(
                        tag=htag, date=hdate, time=htime
                    )
                elif htag.endswith("timeStep"):
                    hline = '<{tag} unit="{unit}"/>\n'.format(
                        tag=htag, unit=hval)
                else:
                    hline = paramline.format(tag=htag, param=hval)
                hlines.append(3 * "\t" + hline)
            hlines.append(2 * "\t" + "</header>\n")
            f.writelines(hlines)

            # write timeseries
            dates = o.reset_index()["index"].apply(
                lambda s: pd.datetime.strftime(s, "%Y-%m-%d")
            )
            times = o.reset_index()["index"].apply(
                lambda s: pd.datetime.strftime(s, "%H:%M:%S")
            )
            # set date and time attributes
            events = (
                2 * "\t"
                + '<event date="'
                + dates.values
                + '" time="'
                + times.values
            )
            # loop through columns and add to event
            for icol in o.columns:
                val = o[icol].astype(str)
                events += (
                    '" {}="'.format(icol)
                    + val.values
                )
            # close event
            events += '"/>\n'
            # write to file
            f.writelines(events)
            # end series
            f.write("\t" + "</series>\n")
        # end Timeseries
        f.write("</TimeSeries>\n")
