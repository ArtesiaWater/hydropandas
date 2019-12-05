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

    def obs_per_year(self, col="Stand_m_tov_NAP"):
        pblist = {o.name: o.obs_per_year(col=col) for o in self._obj.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def consecutive_obs_years(self, min_obs=12, col="Stand_m_tov_NAP"):
        """ get the number of consecutive years with more than a minimum of
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
            pblist = {o.name: o.consecutive_obs_years(
                min_obs=self._obj.loc[o.name, min_obs], col=col) for o in self._obj.obs}
        else:
            pblist = {o.name: o.consecutive_obs_years(
                min_obs=min_obs, col=col) for o in self._obj.obs}
        df = pd.DataFrame.from_dict(pblist)
        return df

    def mean_in_period(self, tmin=None, tmax=None, col="Stand_m_tov_NAP"):
        if tmin is None:
            tmin = self._obj.dates_first_obs.min()
        if tmax is None:
            tmax = self._obj.dates_last_obs.max()

        return self._obj.obs.apply(lambda o: o.loc[tmin:tmax, col].mean())

    def get_no_of_observations(self, column_name='Stand_m_tov_NAP',
                               after_date=None, before_date=None):
        """get number of non-nan values of a column in the observation df

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
        pandas series with the number of observaitons for each row in the obs
        collection.
        """
        no_obs = []
        for o in self._obj.obs.values:
            try:
                no_obs.append(o[after_date:before_date]
                              [column_name].dropna().count())
            except KeyError:
                no_obs.append(0)

        self._obj['no_obs'] = no_obs

        return self._obj['no_obs']

    def get_seasonal_stat(self, column_name='Stand_m_tov_NAP', stat='mean',
                          winter_months=[1, 2, 3, 4, 11, 12],
                          summer_months=[5, 6, 7, 8, 9, 10]):
        """get statistics per season

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
            df_list.append(o.get_seasonal_stat(column_name, stat,
                                               winter_months, summer_months))

        return pd.concat(df_list)

    def get_first_last_obs_date(self):
        """adds two columns to the ObsCollection with the date of the first
        and the last measurement


        """

        self._obj['date_first_measurement'] = [o.index.min()
                                               for o in self._obj.obs.values]
        self._obj['date_last_measurement'] = [o.index.max()
                                              for o in self._obj.obs.values]

    def get_min_max(self, obs_column='Stand_m_tov_NAP'):
        """adds two columns to the Obscollection with the minimum and the
        maximum of a column (defined by obs_column)

        returns the absolute minimum and maximum of all observations

        Parameters
        ----------
        obs_column: str, optional
            column name for which min and max is calculated

        Returns
        -------
        min, max: float
            the minimum and maximum of the obs_column in all observations

        """

        self._obj['max'] = [o[obs_column].max() for o in self._obj.obs.values]
        self._obj['min'] = [o[obs_column].min() for o in self._obj.obs.values]

        return(self._obj['min'].min(), self._obj['max'].max())
