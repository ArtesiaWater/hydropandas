# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 12:15:42 2018

@author: Artesia
"""

import os
import time
import zipfile
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import tempfile


def _import_art_tools():
    try:
        import art_tools as art
        return art
    except ModuleNotFoundError as e:
        print(
            'This function is not available please contact Artesia for more information')
        raise(e)


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


def get_files(file_or_dir, ext, unpackdir=None, force_unpack=False,
              preserve_datetime=False):
    """internal method to get list of files with specific
    extension from dirname.

    Parameters
    ----------
    file_or_dir : str
        file or path to data
    ext : str
        extension of filenames to store in list
    force_unpack : bool, optional
        force unzip, by default False
    preserve_datetime : bool, optional
        preserve datetime of unzipped files, by default False
        (useful for checking whether data has changed)

    """
    # check if unpackdir is same as file_or_dir, if same, this can cause
    # problems when the unpackdir still contains zips that will be unpacked
    # again.
    if unpackdir is not None:
        if os.path.normcase(unpackdir) == os.path.normcase(file_or_dir):
            raise ValueError("Please specify a different folder to unpack"
                             " files!")

    # unzip dir
    if file_or_dir.endswith('.zip'):
        zipf = file_or_dir
        if unpackdir is None:
            file_or_dir = tempfile.TemporaryDirectory().name
        else:
            file_or_dir = unpackdir
        unzip_file(zipf, file_or_dir, force=force_unpack,
                   preserve_datetime=preserve_datetime)

    # file_or_dir is directory
    if os.path.isdir(file_or_dir):
        # check for zips in dir
        zip_fnames = [i for i in os.listdir(
            file_or_dir) if i.endswith(".zip")]
        if len(zip_fnames) > 0:
            # unzip zips
            if unpackdir is None:
                dirname = tempfile.TemporaryDirectory().name
            else:
                dirname = unpackdir
            for zipf in zip_fnames:
                unzip_file(os.path.join(file_or_dir, zipf), dirname,
                           force=force_unpack,
                           preserve_datetime=preserve_datetime)
        else:
            dirname = file_or_dir
        # get all files with extension ext
        unzip_fnames = [i for i in os.listdir(
            dirname) if i.endswith(ext)]
    elif os.path.isfile(file_or_dir):
        # file_or_dir is actually an xml
        unzip_fnames = [os.path.basename(file_or_dir)]  # get file name
        dirname = os.path.dirname(file_or_dir)  # get directory path
    else:
        raise NotImplementedError("Cannot parse 'file_or_dir'!")

    return dirname, unzip_fnames


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
