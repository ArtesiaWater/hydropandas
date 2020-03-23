from tqdm import tqdm
import arctic
from timeit import default_timer


def read_arctic(connstr, libname, ObsClass, verbose=False):
    """Read all timeseries from Arctic database.

    Parameters
    ----------
    connstr : str
        connection string to MongoDB instance
    libname : str
        name of library to read data from
    ObsClass : Obs
        type of Obs to load data as
    verbose : bool, optional
        progressbar and timing info, by default False

    Returns
    -------
    obs_list : list of Obs
        list of Obs

    """
    arc = arctic.Arctic(connstr)
    lib = arc.get_library(libname)

    start = default_timer()
    obs_list = []
    rows_read = 0
    for sym in (tqdm(lib.list_symbols())
                if verbose else lib.list_symbols()):
        item = lib.read(sym)
        o = ObsClass(item.data, name=sym, meta=item.metadata)
        obs_list.append(o)
        if verbose:
            rows_read += len(item.data.index)
    end = default_timer()
    if verbose:
        print("Symbols: {0:.0f}  "
              "Rows: {1:,.0f}  "
              "Time: {2:.2f}s  "
              "Rows/s: {3:,.1f}".format(len(lib.list_symbols()),
                                        rows_read,
                                        (end - start),
                                        rows_read / (end - start)))

    return obs_list
