import numpy as np
from scipy import interpolate

from . import accessor, util


@accessor.register_obscollection_accessor("geo")
class GeoAccessor:
    def __init__(self, oc_obj):
        self._obj = oc_obj

    def get_bounding_box(self, xcol='x', ycol='y', buffer=0):
        """returns the bounding box of all observations

        Parameters
        ----------
        xcol : str, optional
            column name with x values
        ycol : str, optional
            column name with y values
        buffer : int or float, optional
            add a buffer around the bouding box from the observations

        Returns
        -------
        xmin, ymin, xmax, ymax
            coordinates of bouding box

        """

        xmin = self._obj._obj[xcol].min() - buffer
        xmax = self._obj._obj[xcol].max() + buffer
        ymin = self._obj._obj[ycol].min() - buffer
        ymax = self._obj._obj[ycol].max() + buffer

        return (xmin, ymin, xmax, ymax)

    def get_extent(self, xcol='x', ycol='y', buffer=0):
        """returns the extent of all observations

        Parameters
        ----------
        xcol : str, optional
            column name with x values
        ycol : str, optional
            column name with y values
        buffer : int or float, optional
            add a buffer around the bouding box from the observations

        Returns
        -------
        xmin, xmax, ymin, ymax
            coordinates of bouding box

        """

        xmin = self._obj[xcol].min() - buffer
        xmax = self._obj[xcol].max() + buffer
        ymin = self._obj[ycol].min() - buffer
        ymax = self._obj[ycol].max() + buffer

        return (xmin, xmax, ymin, ymax)

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326',
                    add_to_meta=True, add_to_df=True):
        """get lattitude and longitude from x and y attributes

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon
        add_to_meta : boolean, optional
            if True the new coordinates are added to the meta dictionary
        add_to_df : boolean, optional
            if True lon and lat columns are added to the ObsCollection


        """

        for o in self._obj.obs.values:
            o.geo.get_lat_lon(in_epsg, out_epsg, add_to_meta)

        if add_to_df:
            self._obj.add_meta_to_df('lat')
            self._obj.add_meta_to_df('lon')

    def get_nearest_point(self, obs_collection2=None, gdf2=None,
                          xcol_obs1='x', ycol_obs1='y',
                          xcol_obs2='x', ycol_obs2='y', verbose=False):
        """get nearest point of another obs collection for each
        point in the current obs collection.

        Parameters
        ----------
        obs_collection2 : ObsCollection, optional
            collection of observations of which the nearest point is found
        gdf2 : GeoDataFrame, optional
            dataframe to look for nearest observation point
        xcol_obs1 : str, optional
            x column in self._obj used to get geometry
        ycol_obs1 : str, optional
            y column in self._obj used to get geometry
        xcol_obs2 : str, optional
            x column in obs_collection2 used to get geometry
        ycol_obs2 : str, optional
            y column in self._obj used to get geometry
        verbose : boolean, optional
            Print additional information to the screen (default is False).

        Returns
        -------
        pandas.DataFrame
            with columns 'nearest point' and 'distance nearest point'
        """

        from shapely.ops import nearest_points

        gdf1 = self._obj.to_gdf(xcol=xcol_obs1, ycol=ycol_obs1)
        if obs_collection2 is not None:
            gdf2 = obs_collection2.to_gdf(xcol=xcol_obs2, ycol=ycol_obs2)
        elif gdf2 is None:
            raise ValueError('obs_collecction2 or gdf2 should be defined')

        pts_gdf2 = gdf2.geometry.unary_union

        def nearest_point(point_gdf1, pts=pts_gdf2):
            # find the nearest point and return the corresponding Place value
            nearest_point_gdf2 = nearest_points(point_gdf1, pts_gdf2)[1]
            nearest = gdf2[gdf2.geometry == nearest_point_gdf2]
            return nearest.index.values[0]

        def distance_nearest_point(point_gdf1, pts=pts_gdf2):
            # find the nearest point and return the corresponding Place value
            nearest_point_gdf2 = nearest_points(point_gdf1, pts_gdf2)[1]
            distance = point_gdf1.distance(nearest_point_gdf2)
            return distance

        gdf1['nearest point'] = gdf1.apply(
            lambda row: nearest_point(row.geometry), axis=1)
        gdf1['distance nearest point'] = gdf1.apply(
            lambda row: distance_nearest_point(row.geometry), axis=1)

        return gdf1[['nearest point', 'distance nearest point']]

    def get_distance_to_point(self, point, xcol='x', ycol='y'):
        """get distance of every observation to a point.

        Parameters
        ----------
        point : shapely.geometry.point.Point
            point geometry
        xcol : str, optional
            x column in self._obj used to get x coordinates
        ycol : str, optional
            y column in self._obj used to get y coordinates

        Returns
        -------
        pd.Series
            distance to the point for every observation in self._obj

        """

        gdf = self._obj[[xcol, ycol]].to_gdf(xcol=xcol, ycol=ycol)

        return gdf.distance(point)

    def within_extent(self, extent, inplace=False):
        """Slice ObsCollection by extent

        Parameters
        ----------
        extent : tuple
            format (xmin, xmax, ymin, ymax), default dis.sr.get_extent() format

        Returns
        -------
        new_oc : obs_collection.ObsCollection
            ObsCollection with observations within extent
        """

        new_oc = self._obj[(self._obj.x > extent[0]) & (self._obj.x < extent[1])
                           & (self._obj.y > extent[2]) & (self._obj.y < extent[3])]
        if inplace:
            self._obj._update_inplace(new_oc)
        else:
            return new_oc

    def within_polygon(self, gdf=None, shapefile=None, inplace=False,
                       **kwargs):
        """Slice ObsCollection by checking if points are within a shapefile

        Parameters
        ----------
        gdf : GeoDataFrame, optional
            geodataframe containing a single polygon
        shapefile : str, optional
            Not yet implemented
        inplace : bool, default False
            Modify the ObsCollection in place (do not create a new object).
        **kwargs :
            kwargs will be passed to the self._obj.to_gdf() method

        Returns
        -------
        new_oc : obs_collection.ObsCollection
            ObsCollection with observations within polygon
        """

        if gdf is not None:
            if gdf.shape[0] != 1:
                raise NotImplementedError(
                    'cannot handle zero or multiple polygons')
            gdf_oc = self._obj.to_gdf(**kwargs)
            new_oc = self._obj.loc[gdf_oc.within(gdf.geometry.values[0])]
        elif shapefile is not None:
            raise NotImplementedError('shapefiles are not yet implemented')
        else:
            raise ValueError('shapefile or gdf must be specified')

        if inplace:
            self._obj._update_inplace(new_oc)
        else:
            return new_oc

    def get_surface_level(self, xcol='x', ycol='y', buffer=10.,
                          add_to_oc=False, if_exists='error', **kwargs):
        """get maaiveld at the observation points in the observation collection

        Parameters
        ----------
        xcol : str, optional
            column name with x coordinates, by default 'x'
        ycol : str, optional
            column name with y coordinates, by default 'y'
        buffer: int or float, optional
            buffer used to get surrounding ahn values
        add_to_oc : bool, optional
            if True the maaiveld is added to the observation collection, default is False
        if_exists : str, optional
            what to do if an observation point already has a maaiveld, options:
            'error', 'replace' or 'keep', by default 'error'
        **kwargs:
            are passed to art.get_ahn_within_extent() function

        Returns
        -------
        np.array or None
            list of ahn values at observation locations or None (if add_to_oc = True)
        """

        if self._obj.shape[0] < 2:
            raise NotImplementedError(
                'this method will probabaly not work on collections with 0 or 1 observation points')

        # attempt art_tools import
        art = util._import_art_tools()

        # get x and y values from oc_col
        xp = self._obj[xcol].values.astype(float)
        yp = self._obj[ycol].values.astype(float)

        extent = self._obj.geo.get_extent(buffer=buffer)

        ahn = art.get_ahn_within_extent(extent, **kwargs)

        z = art.rasters.get_values(ahn)

        # use griddata (is slow, buit can handle NaNs)
        xc, yc = art.rasters.get_xy_mid(ahn)
        xc, yc = np.meshgrid(xc, yc)
        mask = ~np.isnan(z)
        points = np.column_stack((xc[mask], yc[mask]))
        zp = interpolate.griddata(points, z[mask], np.column_stack((xp, yp)))

        if add_to_oc:
            if if_exists == 'error' and 'maaiveld' in self._obj.columns:
                raise KeyError(
                    "maaiveld already in columns set if_exists to 'keep' or 'replace' to overwrite")
            elif if_exists == 'replace':
                self._obj['maaiveld'] = np.nan
            elif 'maaiveld' not in self._obj.columns:
                self._obj['maaiveld'] = np.nan
            obs_new_maaiveld = self._obj['maaiveld'].isna()
            self._obj.loc[obs_new_maaiveld, 'maaiveld'] = zp[obs_new_maaiveld]
            for o in self._obj.loc[obs_new_maaiveld, 'obs'].values:
                o.maaiveld = self._obj.loc[o.name, 'maaiveld']
                o.meta['maaiveld'] = o.maaiveld

            return None

        else:
            return zp


@accessor.register_obs_accessor("geo")
class GeoAccessorObs:
    def __init__(self, obs):
        self._obj = obs

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326',
                    add_to_meta=True):
        """get lattitude and longitude from x and y attributes

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon
        add_to_meta : boolean, optional
            if True the new coordinates are added to the meta dictionary

        Returns
        -------
        lon, lat : longitude and lattitude of x, y coordinates

        """

        from pyproj import Proj, transform

        inProj = Proj(init=in_epsg)
        outProj = Proj(init=out_epsg)

        if np.isnan(self._obj.x) or np.isnan(self._obj.y):
            lon, lat = np.nan, np.nan
        else:
            lon, lat = transform(inProj, outProj, self._obj.x, self._obj.y)

        if add_to_meta:
            self._obj.meta['lon'] = lon
            self._obj.meta['lat'] = lat

        return lat, lon

    def get_surface_level(self, buffer=10, **kwargs):
        """returns maaiveld at observation point

        Parameters
        ----------
        buffer: int or float, optional
            buffer used to get surrounding ahn values
        **kwargs:
            are passed to art.get_ahn_within_extent() function

        Returns
        -------
        zp: float
            ahn value at location
        """

        # attempt art_tools import
        art = util._import_art_tools()

        extent = [self._obj.x - buffer, self._obj.x + buffer,
                  self._obj.y - buffer, self._obj.y + buffer]
        ds = art.get_ahn_within_extent(extent, **kwargs)
        z = art.rasters.get_values(ds)
        xc, yc = art.rasters.get_xy_mid(ds)
        xc, yc = np.meshgrid(xc, yc)
        mask = ~np.isnan(z)
        points = np.column_stack((xc[mask], yc[mask]))

        zp = float(interpolate.griddata(
            points, z[mask], ((self._obj.x, self._obj.y))))

        return zp
