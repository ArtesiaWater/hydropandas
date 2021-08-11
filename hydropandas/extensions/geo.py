import numpy as np
import pandas as pd

from . import accessor


@accessor.register_obscollection_accessor("geo")
class GeoAccessor:
    def __init__(self, oc_obj):
        self._obj = oc_obj

    def get_bounding_box(self, xcol='x', ycol='y', buffer=0):
        """returns the bounding box of all observations.

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

        xmin = self._obj[xcol].min() - buffer
        xmax = self._obj[xcol].max() + buffer
        ymin = self._obj[ycol].min() - buffer
        ymax = self._obj[ycol].max() + buffer

        return (xmin, ymin, xmax, ymax)

    def get_extent(self, xcol='x', ycol='y', buffer=0):
        """returns the extent of all observations.

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

    def set_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326',
                    add_to_meta=True):
        """create columns with lat and lon values of the observation points.

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon
        add_to_meta : bool, optional
            if True the lat and lon values are added to the observation meta
            dictionary. The default is True.

        Returns
        -------
        None.
        """

        df_lat_lon = self._obj.geo.get_lat_lon(in_epsg, out_epsg)
        for iname in df_lat_lon.index:
            self._obj._set_metadata_value(iname, 'lat',
                                          df_lat_lon.loc[iname,
                                                         'lat'], add_to_meta)
            self._obj._set_metadata_value(iname, 'lon',
                                          df_lat_lon.loc[iname,
                                                         'lon'], add_to_meta)

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326'):
        """get lattitude and longitude from x and y attributes.

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon

        Returns
        -------
        pandas.DataFrame
            with columns 'lat' and 'lon'
        """

        df_lat_lon = pd.DataFrame(
            index=self._obj.index, columns=['lat', 'lon'])
        for iname in self._obj.index:
            o = self._obj.loc[iname, 'obs']
            df_lat_lon.loc[iname, ['lat', 'lon']
                           ] = o.geo.get_lat_lon(in_epsg, out_epsg)

        return df_lat_lon

    def get_nearest_point(self, obs_collection2=None, gdf2=None,
                          xcol_obs1='x', ycol_obs1='y',
                          xcol_obs2='x', ycol_obs2='y'):
        """get nearest point of another obs collection for each point in the
        current obs collection.

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
            return nearest.index[0]

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

    def _get_nearest_geometry(self, gdf=None,
                              xcol_obs='x', ycol_obs='y',
                              geometry_type='polygon',
                              multiple_geometries='error'):
        """get nearest geometry for each point in the obs collection. Function
        works for line and polygon geometries.

        Parameters
        ----------
        gdf : GeoDataFrame
            dataframe with polygon features
        xcol_obs : str, optional
            x column in self._obj used to get geometry
        ycol_obs : str, optional
            y column in self._obj used to get geometry
        geometry_type : str, optional
            geometry type, can be 'polygon' or 'line'.
        multiple_geometries : str, optional
            keyword on how to deal with multiple geometries being nearest.
            Options are:
                'error' -> raise a ValueError
                'keep_all' -> return the indices of multiple geometries as a
                string seperated by a comma
                'keep_first' -> return the index of the first geometry

        Returns
        -------
        pandas.DataFrame
            with columns 'nearest geometry' and 'distance nearest geometry'
        """

        gdf_obs = self._obj.to_gdf(xcol=xcol_obs, ycol=ycol_obs)
        gdf = gdf.copy()
        for i, point in gdf_obs.geometry.items():
            distances = [point.distance(pol) for pol in gdf.geometry.values]
            if (np.array(distances) == np.min(distances)).sum() > 1:
                if multiple_geometries == 'error':
                    raise ValueError(f'multiple {geometry_type}s are nearest')
                elif multiple_geometries == 'keep_all':
                    ids = []
                    for i_min in np.where(np.array(distances) == np.min(distances))[0]:
                        ids.append(gdf.index[i_min])
                    gdf_obs.loc[i, f'nearest {geometry_type}'] = ', '.join(ids)
                    gdf_obs.loc[i, f'distance nearest {geometry_type}'] = np.min(distances)
                elif multiple_geometries == 'keep_first':
                    gdf_obs.loc[i, f'nearest {geometry_type}'] = gdf.iloc[np.argmin(
                    distances)].name
                    gdf_obs.loc[i, f'distance nearest {geometry_type}'] = np.min(distances)
                else:
                    raise ValueError(f'invalid value for multiple_geometries -> {multiple_geometries}')
            else:
                gdf_obs.loc[i, f'nearest {geometry_type}'] = gdf.iloc[np.argmin(
                    distances)].name
                gdf_obs.loc[i, f'distance nearest {geometry_type}'] = np.min(distances)

        return gdf_obs[[f'nearest {geometry_type}', f'distance nearest {geometry_type}']]

    def get_nearest_line(self, gdf=None,
                         xcol_obs='x', ycol_obs='y',
                         multiple_lines='error'):
        """get nearest line for each point in the obs collection. Function
        calls the nearest_polygon function.


        Parameters
        ----------
        gdf : GeoDataFrame
            dataframe with line features
        xcol_obs : str, optional
            x column in self._obj used to get geometry
        ycol_obs : str, optional
            y column in self._obj used to get geometry
        multiple_lines : str, optional
            keyword on how to deal with multiple lines being nearest.
            Options are:
                'error' -> raise a ValueError
                'keep_all' -> return the indices of multiple lines as a
                string seperated by a comma
                'keep_first' -> return the index of the first line

        Returns
        -------
        pandas.DataFrame
            with columns 'nearest polygon' and 'distance nearest polygon'
        """
        return self._get_nearest_geometry(gdf=gdf,
                                               xcol_obs=xcol_obs,
                                               ycol_obs=ycol_obs,
                                               multiple_geometries=multiple_lines,
                                               geometry_type='line')

    def get_nearest_polygon(self, gdf=None,
                            xcol_obs='x', ycol_obs='y',
                            multiple_polygons='error'):
        """get nearest polygon for each point in the obs collection. Function
        also works for lines instead of polygons


        Parameters
        ----------
        gdf : GeoDataFrame
            dataframe with polygon features
        xcol_obs : str, optional
            x column in self._obj used to get geometry
        ycol_obs : str, optional
            y column in self._obj used to get geometry
        multiple_polygons : str, optional
            keyword on how to deal with multiple polygons being nearest.
            Options are:
                'error' -> raise a ValueError
                'keep_all' -> return the indices of multiple polygons as a
                string seperated by a comma
                'keep_first' -> return the index of the first polygon

        Returns
        -------
        pandas.DataFrame
            with columns 'nearest polygon' and 'distance nearest polygon'
        """
        return self._get_nearest_geometry(gdf=gdf,
                                          xcol_obs=xcol_obs,
                                          ycol_obs=ycol_obs,
                                          multiple_geometries=multiple_polygons,
                                          geometry_type='polygon')

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
        """Slice ObsCollection by extent.

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
        """Slice ObsCollection by checking if points are within a shapefile.

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


@accessor.register_obs_accessor("geo")
class GeoAccessorObs:
    def __init__(self, obs):
        self._obj = obs

    def get_lat_lon(self, in_epsg='epsg:28992', out_epsg='epsg:4326'):
        """get lattitude and longitude from x and y attributes.

        Parameters
        ----------
        in_epsg : str, optional
            epsg code of current x and y attributes, default (RD new)
        out_epsg : str, optional
            epsg code of desired output, default lat/lon

        Returns
        -------
        lon, lat : longitude and lattitude of x, y coordinates
        """

        from pyproj import Proj, transform

        inProj = Proj(in_epsg)
        outProj = Proj(out_epsg)

        if np.isnan(self._obj.x) or np.isnan(self._obj.y):
            lat, lon = np.nan, np.nan
        else:
            lat, lon = transform(inProj, outProj, self._obj.x, self._obj.y)

        return lat, lon
