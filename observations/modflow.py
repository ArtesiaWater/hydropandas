import numpy as np
import scipy.interpolate as intp


def findrowcolumn(xpts, ypts, ml):
    """
    find row and column number for points in an irregular grid. Adapted from Flopy.

    Parameters
    ----------
    xpts : 1d array
        x-coordinates of points
    ypts : 1d array
        y-coordinates of points
    ml : flopy.modflow.mf.Modflow
        modflow model

    Returns
    -------
    r, c: arrays with dtype: int
        array containing row and column indices of the points

    """
    xedge, yedge = ml.modelgrid.xyedges

    if not isinstance(xpts, np.ndarray):
        xpts = np.asarray(xpts)
    if not isinstance(ypts, np.ndarray):
        ypts = np.asarray(ypts)
    if not isinstance(xedge, np.ndarray):
        xedge = np.asarray(xedge)
    if not isinstance(yedge, np.ndarray):
        yedge = np.asarray(yedge)

    fx = intp.interp1d(xedge, range(len(xedge)), bounds_error=False)
    fy = intp.interp1d(yedge, range(len(yedge)), bounds_error=False)

    c = np.asarray(np.floor(fx(xpts)), dtype=int)
    r = np.asarray(np.floor(fy(ypts)), dtype=int)

    return r, c


def get_model_layer(x, y, z, ml, dis, zgr=None):
    """get index of model layer based on elevation
    at a specific location.

    Parameters
    ----------
    x : np.array
        x-coordinate of locations
    y : np.array
        y-coordinate of locations
    z : np.array
        elevations corresponding to x, y locations
    ml : flopy.modflow.mf.Modflow
        modflow model
    zgr : np.3darray, optional
        3D array containing layer elevations for if dis
        does not contain this information. (the default is
        None, which uses dis to determine layer elevations)

    Returns
    -------
    ilay: np.array
        array containing 0-indexed model layer for each point

    TODO:
        - avoid loop if model layer elevations are equal everywhere
        -
    """

    r, c = findrowcolumn(x, y, ml)

    if zgr is None:
        zgr = np.concatenate(
            [dis.top.array[np.newaxis], dis.botm.array], axis=0)

    ilay = np.nan * np.ones(x.shape, dtype=np.int)

    for i in range(x.shape[0]):
        ir = r[i]
        ic = c[i]

        if np.isnan(z[i]):
            continue
        elif (ir >= 0) & (ic >= 0):
            zvec = zgr[:, ir, ic]
            ilay[i] = int(np.interp(z[i], zvec[::-1], np.arange(len(zvec))[::-1],
                                    left=-999, right=999))

    return ilay


def get_pb_modellayer(x, y, ftop, fbot, ml, zgr=None, verbose=False):
    """get index of model layer based on filterdepth of
    piezometer.

    Parameters
    ----------
    x : np.array
        x-coordinate of piezometers
    y : np.array
        y-coordinate of piezometers
    ftop : np.array
        top elevation of filter
    fbot : np.array
        bottom elevation of filter
    ml : flopy.modflow.mf.Modflow
        modflow model
    zgr : np.3darray, optional
        array containing model layer elevation
        information (the default is None, which
        gets this information from the dis object)
    verbose : bool, optional
        print information about how layer is
        selected when piezometer is screened in
        multiple layers (the default is False,
        which prints no info)

    Returns
    -------
    ilay: np.array
        array containing index of model layer per piezometer

    TODO:
        - some things are done double in this function and in get_model_layer, necessary?
        - speed up if model layer elevation is constant everywhere?

    """
    dis = ml.get_package('DIS')
    
    ilay_ftop = get_model_layer(x, y, ftop, ml, dis, zgr=zgr)
    ilay_fbot = get_model_layer(x, y, fbot, ml, dis, zgr=zgr)
    
    

    if zgr is None:
        zgr = np.concatenate(
            [dis.top.array[np.newaxis], dis.botm.array], axis=0)

    r, c = findrowcolumn(x, y, ml)

    ilay = np.nan * np.ones(x.shape, dtype=np.int)

    for i in range(x.shape[0]):
        ir = r[i]
        ic = c[i]

        if (ir >= 0) & (ic >= 0):
            zvec = zgr[:, ir, ic]

            fti = ilay_ftop[i]
            fbi = ilay_fbot[i]

            if np.isnan(fti) and np.isnan(fbi):
                continue
            elif np.isnan(fti):
                ilay[i] = fbi
            elif np.isnan(fbi):
                ilay[i] = fti
            # Warn if piezometer is on boundary. Piezometer in layer in which majority of screen is located
            elif fbi != fti:
                fbi = int(fbi)
                fti = int(fti)
                if verbose:
                    print("filter on layer boundary:")
                    print("-layers: {0}, {1}".format(fbi, fti))
                    print("-layer elev in between: {0:.2f}".format(zvec[fbi]))
                    print("-fb, ft: {0:.2f}, {1:.2f}".format(fbot[i], ftop[i]))
                    print("-length in layer: {0:.2f}, {1:.2f}".format(
                          zvec[fbi] - fbot[i], ftop[i] - zvec[fbi]))
                if fti - fbi > 2:
                    ilay[i] = np.nan  # set to unknown
                    if verbose:
                        print("Piezometer filter spans {} layers. "
                              "Not sure which layer to select".format(fti - fbi))
                elif fti - fbi == 2:  # if filter spans 3 layers
                    ilay[i] = fbi + 1  # use middle layer
                elif np.abs(fbot[i] - zvec[fbi]) > np.abs(ftop[i] - zvec[fbi]):
                    ilay[i] = fbi
                    if verbose:
                        print("selected layer:", fbi)
                elif np.abs(fbot[i] - zvec[fbi]) < np.abs(ftop[i] - zvec[fbi]):
                    ilay[i] = fti
                    if verbose:
                        print("selected layer:", fti)
            else:
                ilay[i] = fbi

    return ilay
