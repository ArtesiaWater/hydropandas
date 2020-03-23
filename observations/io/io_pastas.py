# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 12:15:42 2018

@author: Artesia
"""

import pandas as pd


def get_model_checks(results_df, evp_threshold=70,
                     verdampingsfactor_min=0.5, verdampingsfactor_max=2.0,
                     timeseries_memory_factor=5):

    check_df = pd.DataFrame(index=results_df.index)
    check_df['EVP>{}%'.format(evp_threshold)
             ] = results_df['EVP [%]'] > evp_threshold
    check_df['{} < Evaporation-factor < {}'.format(verdampingsfactor_min, verdampingsfactor_max)] = (
        results_df['Evaporation-factor [-]'] > verdampingsfactor_min) & (results_df['Evaporation-factor [-]'] < verdampingsfactor_max)
    check_df['Timeseries length > {} x Memory'.format(
        timeseries_memory_factor)] = results_df['Length calibration period [days]'] > results_df['Memory [days]'] * timeseries_memory_factor
    check_df['Passed all checks'] = check_df[['EVP>{}%'.format(evp_threshold), '{} < Evaporation-factor < {}'.format(
        verdampingsfactor_min, verdampingsfactor_max), 'Timeseries length > {} x Memory'.format(timeseries_memory_factor)]].all(axis=1)

    return check_df


def get_pr_results(pr):
    """

    Notes
    -----
    This function assumes that the project contains 1 model for every oseries




    """
    results_df = pd.DataFrame(index=pr.oseries.index)
    results_df['number of observations used in calibration'] = get_number_of_observations_in_calibration(
        pr)
    results_df = get_meteostation_evap(pr, results_df)
    results_df = get_meteostation_precip(pr, results_df)
    results_df['EVP [%]'] = get_evp(pr)
    results_df['Evaporation-factor [-]'], results_df['std Evaporation-factor [-]'] = get_rch_fct(
        pr)
    results_df['Constant [m NAP]'], results_df['Constant [m]'] = get_drainage_basis(
        pr)
    results_df['Memory [days]'] = get_memory(pr)
    results_df['Length calibration period [days]'] = get_calibration_period(pr)

    return results_df.round(2)


def get_number_of_observations_in_calibration(pr):

    return [ml.oseries_calib.index.size for ml in pr.models.values()]


def get_meteostation_evap(pr, results_df):

    evap_s = pr.get_nearest_stresses(kind="evap")
    results_df.loc[evap_s.index, 'meteostation verdamping'] = evap_s.values

    return results_df


def get_meteostation_precip(pr, results_df):

    evap_s = pr.get_nearest_stresses(kind="prec")
    results_df.loc[evap_s.index, 'meteostation neerslag'] = evap_s.values

    return results_df


def get_evp(pr, colname='EVP [%]'):

    return pr.get_statistics(['evp'])


def get_rch_fct(pr, add_stderr=True):

    parameters_df = pr.get_parameters(['recharge_f'])
    if add_stderr:
        rch_f_std_s = pr.get_parameters(['recharge_f'], param_value='stderr')
        return [-parameters_df, rch_f_std_s]
    else:
        return -parameters_df


def get_drainage_basis(pr, add_stderr=True):

    parameters_df = pr.get_parameters(['constant_d'])
    if add_stderr:
        db_std_s = pr.get_parameters(['constant_d'], param_value='stderr')
        return [parameters_df, db_std_s]
    else:
        return parameters_df


def get_memory(pr, stresstype='recharge',
               cutoff=0.9):

    memory = []
    for m in pr.models.values():
        rfunc = m.stressmodels[stresstype].rfunc
        memory.append(rfunc.get_tmax(
            m.get_parameters(stresstype), cutoff=cutoff))

    return memory


def get_calibration_period(pr):

    calibration_period = []
    for m in pr.models.values():
        calibration_period.append((m.get_tmax() - m.get_tmin()).days)

    return calibration_period


def read_project(pr, ObsClass, rename_dic={}):
    """Read pastas.Project into ObsCollection.

    Parameters
    ----------
    pr : pastas.Project
        Project to read
    ObsClass : Obs
        ObsClass to read data as, usually GroundwaterObs
    rename_dic : dict, optional
        rename columns in Project oseries dictionary, by default empty dict

    Returns
    -------
    list : list of Obs
        list of Obs containing oseries data

    """
    obs_list = []
    for index, row in pr.oseries.iterrows():
        metadata = row.to_dict()
        for key in rename_dic.keys():
            if key in metadata.keys():
                metadata[rename_dic[key]] = metadata.pop(key)

        s = pd.DataFrame(metadata.pop('series').series_original)
        s.rename(columns={index: 'stand_m_tov_nap'}, inplace=True)

        keys_o = ['name', 'x', 'y', 'locatie', 'filternr',
                  'metadata_available', 'maaiveld', 'meetpunt',
                  'bovenkant_filter', 'onderkant_filter']
        meta_o = {k: metadata[k] for k in keys_o if k in metadata}

        o = ObsClass(s, meta=metadata, **meta_o)
        obs_list.append(o)
    return obs_list
