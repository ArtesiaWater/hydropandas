# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 12:15:42 2018

@author: Artesia
"""

import os
import time
import zipfile

import flopy
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
from flopy.utils import Util2d, Util3d, reference
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shapely.geometry import Point, Polygon
from shapely.strtree import STRtree


def colorbar_inside(mappable=None, ax=None, width=0.2, height="90%", loc=5, **kw):
    if ax is None:
        ax = plt.gca()
    cax = inset_axes(ax, width=width, height=height, loc=loc)
    cb = plt.colorbar(mappable, cax=cax, ax=ax, **kw)
    if loc == 1 or loc == 4 or loc == 5:
        cax.yaxis.tick_left()
        cax.yaxis.set_label_position("left")
    return cb


def title_inside(title, ax=None, x=0.5, y=0.98, **kwargs):
    if ax is None:
        ax = plt.gca()
    return ax.text(x, y, title,
                   horizontalalignment='center', verticalalignment='top',
                   transform=ax.transAxes, **kwargs)


def geodataframe2grid(mf, shp_in):
    # cut a geodataframeby a modflow-grid
    geom_col_name = shp_in._geometry_column_name

    # make a polygon for each of the grid-cells
    # this may result in lines on the edge of the polygon counting double
    polys = []
    for row in range(mf.sr.nrow):
        for col in range(mf.sr.ncol):
            vert = mf.sr.get_vertices(row, col)
            pol = Polygon(vert)
            pol.row = row
            pol.col = col
            polys.append(pol)

    s = STRtree(polys)

    shp_list = []
    # cut the lines with the grid
    for index, row in shp_in.iterrows():
        g = row[geom_col_name]
        result = s.query(g)
        for r in result:
            i = g.intersection(r)
            add_shapes_to_list(i, g.geometryType(),
                               geom_col_name, row, shp_list, r)
    shp_out = gpd.GeoDataFrame(pd.DataFrame(shp_list), geometry=geom_col_name)
    return shp_out


def add_shapes_to_list(i, geometryType, geom_col_name, row, shp_list, r):
    if geometryType == 'LineString':
        if not i.is_empty:
            it = i.geometryType()
            if it == 'GeometryCollection':
                for im in i.geoms:
                    add_shapes_to_list(
                        im, geometryType, geom_col_name, row, shp_list, r)
            elif it == 'MultiLineString':
                for im in i.geoms:
                    add_shapes_to_list(
                        im, geometryType, geom_col_name, row, shp_list, r)
            elif it == 'LineString':
                rown = row.copy()
                rown[geom_col_name] = i
                rown['row'] = r.row
                rown['col'] = r.col
                shp_list.append(rown)
            elif it == 'Point':
                # endpoint of the linestring is on the cell-edge
                pass
            elif it == 'MultiPoint':
                # mutiple endpoints of the linestring are on the cell-edge
                pass
            else:
                raise NotImplementedError(
                    'geometryType ' + it + ' not yet supprted in geodataframe2grid')
    elif geometryType == 'Polygon' or geometryType == 'MultiPolygon':
        it = i.geometryType()
        if it == 'GeometryCollection':
            for im in i.geoms:
                add_shapes_to_list(
                    im, geometryType, geom_col_name, row, shp_list, r)
        elif it == 'MultiPolygon':
            for im in i.geoms:
                add_shapes_to_list(
                    im, geometryType, geom_col_name, row, shp_list, r)
        elif it == 'Polygon':
            rown = row.copy()
            rown[geom_col_name] = i
            rown['row'] = r.row
            rown['col'] = r.col
            shp_list.append(rown)
        elif it == 'Point':
            # endpoint of the polygon is on the cell-edge
            pass
        elif it == 'LineString' or it == 'MultiLineString':
            # one of the edges of the polygon is on a cell-egde
            pass
        else:
            raise NotImplementedError(
                'geometryType ' + it + ' not yet supprted in geodataframe2grid')
    else:
        raise NotImplementedError(
            'geometryType ' + geometryType + ' not yet supprted in geodataframe2grid')


def unzip_file(src, dst, force=False, preserve_datetime=False):
    """Unzip file

    Parameters
    ----------
    src : str
        source zip file
    dst : str
        destination directory
    force : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        use date of the zipfile for the destination file

    Returns
    -------
    int
        1 of True

    """
    if os.path.exists(dst):
        if not force:
            print(
                "File not unzipped. Destination already exists. Use 'force=True' to unzip.")
            return
    if preserve_datetime:
        zipf = zipfile.ZipFile(src, 'r')
        for f in zipf.infolist():
            zipf.extract(f, path=dst)
            date_time = time.mktime(f.date_time + (0, 0, -1))
            os.utime(os.path.join(dst, f.filename), (date_time, date_time))
        zipf.close()
    else:
        zipf = zipfile.ZipFile(src, 'r')
        zipf.extractall(dst)
        zipf.close()
    return 1


def df2gdf(df, xcol='x', ycol='y'):
    """Convert DataFrame to a point GeoDataFrame

    Parameters
    ----------
    df : pd.DataFrame
        dataframe
    xcol : str
        column name with x values
    ycol : str
        column name with y values

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    gdf = gpd.GeoDataFrame(df.copy(), geometry=[Point(
        (s[xcol], s[ycol])) for i, s in df.iterrows()])
    return gdf


def script_newer_than_output(fscript, foutput):
    if not os.path.exists(foutput):
        return True

    tm_script = os.path.getmtime(fscript)

    if isinstance(foutput, str):
        tm_output = os.path.getmtime(foutput)
    else:
        tm_output = np.zeros(len(foutput))
        for i in range(len(foutput)):
            tm_output[i] = os.path.getmtime(foutput[i])

    return np.all(tm_script > tm_output)


def unzip_changed_files(zipname, pathname, check_time=True, check_size=False,
                        debug=False):
    # Extract each file in a zip-file only when the properties are different
    # With the default arguments this method only checks the modification time
    with zipfile.ZipFile(zipname) as zf:
        infolist = zf.infolist()
        for info in infolist:
            fname = os.path.join(pathname, info.filename)
            extract = False
            if os.path.exists(fname):
                if check_time:
                    tz = time.mktime(info.date_time + (0, 0, -1))
                    tf = os.path.getmtime(fname)
                    if tz != tf:
                        extract = True
                if check_size:
                    sz = info.file_size
                    sf = os.path.getsize(fname)
                    if sz != sf:
                        extract = True
            else:
                extract = True
            if extract:
                if debug:
                    print('extracting {}'.format(info.filename))
                zf.extract(info.filename, pathname)
                # set the correct modification time
                # (which is the time of extraction by default)
                tz = time.mktime(info.date_time + (0, 0, -1))
                os.utime(os.path.join(pathname, info.filename), (tz, tz))


def interp_weights(xy, uv, d=2):
    """Calculate interpolation weights

    Parameters
    ----------
    xy : np.array
        array containing x-coordinates in first column and y-coordinates
        in second column
    uv : np.array
        array containing coordinates at which interpolation weights should
        be calculated, x-data in first column and y-data in second column
    d : int, optional
        dimension of data? (the default is 2, which works for 2D data)

    Returns
    -------
    vertices: np.array
        array containing interpolation vertices

    weights: np.array
        array containing interpolation weights per point

    Reference
    ---------
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids

    """

    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def interpolate(values, vtx, wts):
    """interpolate values at locations defined by vertices and points,
       as calculated by interp_weights function.

    Parameters
    ----------
    values : np.array
        array containing values to interpolate
    vtx : np.array
        array containing interpolation vertices, see interp_weights()
    wts : np.array
        array containing interpolation weights, see interp_weights()

    Returns
    -------
    arr: np.array
        array containing interpolated values at locations as given by
        vtx and wts

    Reference
    ---------
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids

    """

    return np.einsum('nj,nj->n', np.take(values, vtx), wts)
