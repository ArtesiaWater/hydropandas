import datetime as dt

import numpy as np
import pandas as pd

from . import accessor


@accessor.register_obscollection_accessor("stats")
class StatsAccessor:
    def __init__(self, oc_obj):
        self._obj = oc_obj

    @property
    def n_observations(self):
        return self._obj.obs.apply(lambda o: o.shape[0])

    @property
    def dates_first_obs(self):
        return self._obj.obs.apply(lambda o: o.index[0])

    @property
    def dates_last_obs(self):
        return self._obj.obs.apply(lambda o: o.index[-1])

    @property
    def obs_periods(self):
        return self._obj.dates_last_obs - self._obj.dates_first_obs

    def obs_per_year(self, col="stand_m_tov_nap"):
        pblist = {o.name: o.stats.obs_per_year(col=col) for o in self._obj.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def consecutive_obs_years(self, min_obs=12, col="stand_m_tov_nap"):
        """get the number of consecutive years with more than a minimum of
        observations.

        Parameters
        ----------
        min_obs : int or str, optional
            if min_obs is an integer it is the minimum number of observations
            per year. If min_obs is a string it is the column name of the
            obs_collection with minimum number of observation per year per
            observation.
        col : str, optional
            the column of the obs dataframe to get measurements from.

        Returns
        -------
        df : pd.DataFrame
            dataframe with the observations as column, the years as rows,
            and the values are the number of consecutive years.
        """

        if isinstance(min_obs, str):
            pblist = {o.name: o.stats.consecutive_obs_years(
                min_obs=self._obj.loc[o.name, min_obs], col=col) for o in self._obj.obs}
        else:
            pblist = {o.name: o.stats.consecutive_obs_years(
                min_obs=min_obs, col=col) for o in self._obj.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def mean_in_period(self, tmin=None, tmax=None, col="stand_m_tov_nap"):
        """get the mean value of one column (col) in all observations within a
        period defined by tmin and tmax. If both tmin and tmax are None the
        whole period in which there are observations is used.

        Parameters
        ----------
        tmin : datetime, optional
            start of averaging period. The default is None.
        tmax : datetime, optional
            end of averaging period. The default is None.
        col : str, optional
            name of the column in the obs object. The default is "stand_m_tov_nap".

        Returns
        -------
        pd.Series
            mean values for each observation.
        """
        if tmin is None:
            tmin = self._obj.stats.dates_first_obs.min()
        if tmax is None:
            tmax = self._obj.stats.dates_last_obs.max()

        return self._obj.obs.apply(lambda o: o.loc[tmin:tmax, col].mean())

    def get_no_of_observations(self, column_name='stand_m_tov_nap',
                               after_date=None, before_date=None):
        """get number of non-nan values of a column in the observation df.

        Parameters
        ----------
        column_name : str, optional
            column name of the  observation data you want to count
        after_date : dt.datetime, optional
            get the number of observations after this date
        before_date : dt.datetime, optional
            get the number of observations before this date

        Returns
        -------
        pandas series with the number of observations for each row in the obs
        collection.
        """

        def get_num_obs(o, column_name=column_name, after_date=after_date,
                        before_date=before_date):
            try:
                return o[after_date:before_date][column_name].dropna().count()
            except KeyError:
                return 0

        return self._obj.obs.apply(get_num_obs)

    def get_seasonal_stat(self, column_name='stand_m_tov_nap', stat='mean',
                          winter_months=[1, 2, 3, 4, 11, 12],
                          summer_months=[5, 6, 7, 8, 9, 10]):
        """get statistics per season.

        Parameters
        ----------
        column_name : str, optional
            column name of the  observation data you want stats for
        stat : str, optional
            type of statistics, all statisics from df.describe() are available
        winter_months : list of int, optional
            month number of winter months
        summer_months : list of int, optional
            month number of summer months


        Returns
        -------
        DataFrame with stats for summer and winter
        """

        df_list = []
        for o in self._obj.obs.values:
            df_list.append(o.stats.get_seasonal_stat(column_name, stat,
                                                     winter_months, summer_months))

        return pd.concat(df_list)

    def get_first_last_obs_date(self):
        """get the date of the first and the last measurement.

        Returns
        -------
        DataFrame with 2 columns with the dates of the first and the last
        measurement
        """

        date_first_measurement = [o.index.min() for o in self._obj.obs.values]
        date_last_measurement = [o.index.max() for o in self._obj.obs.values]

        first_last_obs = pd.DataFrame(index=self._obj.index,
                                      data={'date_first_measurement': date_first_measurement,
                                            'date_last_measurement': date_last_measurement})

        return first_last_obs

    def get_min(self, column_name='stand_m_tov_nap',
                after_date=None, before_date=None):
        """get the minimum value of every obs object.

        Parameters
        ----------
        column_name: str, optional
            column name for which min is calculated
        after_date : dt.datetime, optional
            get the min value after this date
        before_date : dt.datetime, optional
            get the min value before this date

        Returns
        -------
        pandas series with the minimum of each observation in the obs
        collection.
        """

        def get_min_obs(o, column_name=column_name, after_date=after_date,
                        before_date=before_date):
            try:
                return o[after_date:before_date][column_name].dropna().min()
            except KeyError:
                return np.nan

        return self._obj.obs.apply(get_min_obs)

    def get_max(self, column_name='stand_m_tov_nap',
                after_date=None, before_date=None):
        """get the maximum value of every obs object.

        Parameters
        ----------
        column_name: str, optional
            column name for which max is calculated
        after_date : dt.datetime, optional
            get the max value after this date
        before_date : dt.datetime, optional
            get the max value before this date

        Returns
        -------
        pandas series with the maximum of each observation in the obs
        collection.
        """

        def get_max_obs(o, column_name=column_name, after_date=after_date,
                        before_date=before_date):
            try:
                return o[after_date:before_date][column_name].dropna().max()
            except KeyError:
                return np.nan

        return self._obj.obs.apply(get_max_obs)


@accessor.register_obs_accessor("stats")
class StatsAccessorObs:
    def __init__(self, obs):
        self._obj = obs

    def get_seasonal_stat(self, column_name='stand_m_tov_nap', stat='mean',
                          winter_months=[1, 2, 3, 4, 11, 12],
                          summer_months=[5, 6, 7, 8, 9, 10]):
        """get statistics per season.

        Parameters
        ----------
        column_name : str, optional
            column name of the  observation data you want stats for
        stat : str, optional
            type of statistics, all statisics from df.describe() are available
        winter_months : list of int, optional
            month number of winter months
        summer_months : list of int, optional
            month number of summer months


        Returns
        -------
        winter_stats, summer_stats
            two lists with the statistics for the summer and the winter.
        """

        if self._obj.empty:
            df = pd.DataFrame(index=[self._obj.name], data={'winter_{}'.format(stat): [np.nan],
                                                            'summer_{}'.format(stat): [np.nan]})
        else:
            winter_stat = self._obj.loc[self._obj.index.month.isin(
                winter_months)].describe().loc[stat, column_name]
            summer_stat = self._obj.loc[self._obj.index.month.isin(
                summer_months)].describe().loc[stat, column_name]
            df = pd.DataFrame(index=[self._obj.name], data={'winter_{}'.format(stat): [winter_stat],
                                                            'summer_{}'.format(stat): [summer_stat]})

        return df

    def obs_per_year(self, col):
        if self._obj.empty:
            return pd.Series(dtype=float)
        else:
            return self._obj.groupby(self._obj.index.year).count()[col]

    def consecutive_obs_years(self, col, min_obs=12):

        obs_per_year = self._obj.stats.obs_per_year(col=col)

        # Add missing years
        if obs_per_year.empty:
            # no obs, set series to current year with 0 obs
            obs_per_year_all = pd.Series(
                index=[dt.datetime.now().year], data=0)
        else:
            obs_per_year_all = pd.Series(dtype=float,
                                         index=range(obs_per_year.index[0],
                                                     obs_per_year.index[-1] + 1))
            obs_per_year_all.loc[obs_per_year.index] = obs_per_year

        mask_obs_per_year = obs_per_year_all >= min_obs
        mask_obs_per_year.loc[obs_per_year_all.isna()] = np.nan
        mask_obs_per_year.loc[mask_obs_per_year == 0] = np.nan
        cumsum = mask_obs_per_year.cumsum().fillna(method="pad")
        reset = -cumsum.loc[mask_obs_per_year.isnull()].diff().fillna(cumsum)
        result = mask_obs_per_year.where(
            mask_obs_per_year.notnull(), reset).cumsum()

        return result
