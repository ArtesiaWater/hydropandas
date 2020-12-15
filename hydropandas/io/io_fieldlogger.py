# -*- coding: utf-8 -*-
"""Created on Wed Sep 12 12:15:42 2018.

@author: Artesia
"""
import warnings

import pandas as pd

from .. import observation as obs


def _fieldlogger_groundwater_settings(df, verbose=False):
    """Get default settings for groundwater observations.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe with observation point data

    Returns
    -------
    heading
    name
    subname
    inputfield
    properties
    group
    """

    if verbose:
        print('use default settings for groundwater observations')

    heading = ('use_standard_time;NO\nNAME;INPUTTYPE;HINT'
               '\nStand;number;GWS t.o.v. bkpb (cm)'
               '\nOpmerking;text;maak opmerking...\n')
    name = df.locatie
    subname = df.index
    inputfield = 'Stand|Opmerking'
    properties = df.obs.apply(
        lambda o: "bovenkant peilbuis|{:+.2f} mNAP|"
        "onderkant peilbuis|{:+.2f} mNAP|"
        "Maximaal gemeten|{:+.2f} mNAP|"
        "Gemiddeldeld gemeten|{:+.2f} mNAP|"
        "Minimaal gemeten|{:+.2f} mNAP".format(
            o.meta['bovenkant_filter'],
            o.meta['onderkant_filter'],
            o['stand_m_tov_nap'].max(),
            o['stand_m_tov_nap'].mean(),
            o['stand_m_tov_nap'].min()))
    group = df.name

    return heading, name, subname, inputfield, properties, group


def _fieldlogger_waterlevel_settings(df, verbose=False):
    """Get default  settings for waterlevel observations.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe with observation point data

    Returns
    -------
    name
    subname
    inputfield
    properties
    group
    """

    if verbose:
        print('use default settings for waterlevel observations')

    heading = ('use_standard_time;NO\nNAME;INPUTTYPE;HINT\n'
               'Stand;number;oppw t.o.v. kruinhoogte stuw (cm)\n'
               'Opmerking;text;maak opmerking...\n')
    name = df.index
    subname = df.index
    inputfield = 'Stand|Opmerking'
    properties = df.obs.apply(
        lambda o: "Maximaal gemeten|{:+.2f} mNAP|"
        "Gemiddeldeld gemeten|{:+.2f} mNAP|"
        "Minimaal gemeten|{:+.2f} mNAP".format(
            o['stand_m_tov_nap'].max(),
            o['stand_m_tov_nap'].mean(),
            o['stand_m_tov_nap'].min()))
    group = df.name

    return heading, name, subname, inputfield, properties, group


def fieldlogger_csv_to_obs_list(fname, ObsClass=obs.GroundwaterObs):
    """Read a fieldlogger file into a list of observation objects.

    Parameters
    ----------
    fname : str
        name of the fieldlogger location csv file
    ObsClass : observation class
        type of fieldlogger observations

    Returns
    -------
    obs_list : list
        list of observation objects
    """
    # read header
    fieldlogger_meta = {}
    with open(fname, 'r') as f:
        lcounter = 0
        line = f.readline()
        fieldlogger_meta[line.split(';')[0]] = line.strip().split(';')[1]
        while len(line.split(';')) < 4:
            line = f.readline()
            if line == 'NAME;INPUTTYPE;HINT\n':
                columns = line.strip().split(';')
                inputfields_df = pd.DataFrame(columns=columns)
                line = f.readline()
                i = 0
                while line != 'GROUP;COLOR\n':
                    inputfields_df.loc[i] = line.strip().split(';')
                    line = f.readline()
                    i += 1
                fieldlogger_meta['inputfields'] = inputfields_df
                if line == 'GROUP;COLOR\n':
                    columns = line.strip().split(';')
                    groups_df = pd.DataFrame(columns=columns)
                    line = f.readline()
                    i = 0
                    while len(line.split(';')) < 3:
                        groups_df.loc[i] = line.strip().split(';')
                        line = f.readline()
                        i += 1
                    fieldlogger_meta['groups'] = groups_df
            lcounter += 1

        # read body
        df = pd.read_csv(f, sep=';', names=line.strip('\n').split(';'),
                         header=None, index_col='SUBNAME')

    # create obs objects
    obs_list = []
    for o in df.iterrows():
        name = o[0]
        meta = o[1].to_dict()

        # find coordinates
        if 'XCOOR' in meta and 'YCOOR' in meta:
            x = meta['XCOOR']
            y = meta['YCOOR']
        elif 'LAT' in meta and 'LON' in meta:
            warnings.warn('X and Y coordinates are set to lat/lon values!')
            x = meta['LON']
            y = meta['LAT']
        else:
            raise ValueError('could not find xy coordinates in .csv file')

        if 'PROPERTIES' in meta:
            properties = meta.pop('PROPERTIES').split('|')
            meta.update(dict(zip(properties[::2], properties[1::2])))

        obs_list.append(ObsClass(pd.DataFrame(), name=name,
                                 x=x, y=y,
                                 meta=meta))

    return obs_list, fieldlogger_meta


def df_to_fieldlogger_csv(df, fname, otype=None,
                          use_default_otype_settings=True,
                          heading=None,
                          name=None, subname=None, inputfield=None,
                          properties=None, group=None, group_color='blue',
                          verbose=False):
    """Write a csv file that can be read by the fieldlogger app.

    Notes
    -----
    1. option to change input fields for measurements is not yet added
    2. multiple groups are not yet supported

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe with observation point data
    fname : str
        name of the fieldlogger csv file
    otype : str, optional
        observation type, 'groundwater' and 'surfacewater' are supported
    use_default_otype_settings : boolean, optional
        use the default settings for this measurement type
    heading : str, optional
        heading of the csv, if None use default heading
    name : series or str, optional
        1st column in fieldlogger csv, if None use index of ObsCollection
    subname : series or str, optional
        2nd column in fieldlogger csv, if None use name + _sub
    inputfield : series or str, optional
        3th column in fieldlogger csv, if None use 'Stand|Opmerking'
    properties : series or str, optional
        4th column in fieldlogger csv, if None use empty string
    group : series or str, optional
        5th column in fieldlogger csv, if None do not add group to csv
    group_color : series or str, optional
        color of the group
    """

    # get default settings
    if otype == obs.GroundwaterObs and use_default_otype_settings:
        heading, name, subname, inputfield, properties, group = _fieldlogger_groundwater_settings(
            df, verbose)
    elif otype == obs.WaterlvlObs and use_default_otype_settings:
        heading, name, subname, inputfield, properties, group = _fieldlogger_waterlevel_settings(
            df, verbose)

    # csv heading
    if heading is None:
        heading = ('use_standard_time;NO\nNAME;INPUTTYPE;HINT\n'
                   'Stand;number;GWS t.o.v. bkpb (cm)\n'
                   'Opmerking;text;maak opmerking...\n')

    columns = ['SUBNAME', 'XCOOR', 'YCOOR', 'INPUTFIELD', 'PROPERTIES']

    # add group to heading
    if type(group) == str:
        heading += 'GROUP;COLOR\n{};{}\n'.format(group, group_color)
        columns.append('GROUP')
        if verbose:
            print('add data as group {}'.format(group))
    elif group is not None:
        raise NotImplementedError('multiple groups are not yet supported')
    else:
        if verbose:
            print('no group is specified')

    # csv body
    fieldlogger_df = df[['x', 'y']].copy()
    if name is not None:
        fieldlogger_df.index = name
    fieldlogger_df.index.name = 'NAME'

    if subname is None:
        fieldlogger_df['SUBNAME'] = df.index + '_sub'
    else:
        fieldlogger_df['SUBNAME'] = subname

    if inputfield is None:
        fieldlogger_df['INPUTFIELD'] = 'Stand|Opmerking'
    else:
        fieldlogger_df['INPUTFIELD'] = inputfield

    if properties is None:
        fieldlogger_df['PROPERTIES'] = ''
    else:
        fieldlogger_df['PROPERTIES'] = properties.values

    if group is not None:
        fieldlogger_df['GROUP'] = group

    fieldlogger_df.rename(columns={'x': 'XCOOR', 'y': 'YCOOR'}, inplace=True)
    fieldlogger_df.index.name = 'NAME'

    # write file
    with open(fname, 'w') as fo:
        fo.write(heading)
        fieldlogger_df.to_csv(fo, header=True, columns=columns, sep=';')

    return fieldlogger_df
