import logging
import os
import re
import tempfile
from io import FileIO, TextIOWrapper
from pathlib import Path
from typing import Union
from zipfile import ZipFile

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _read_dino_groundwater_header(f):
    line = f.readline()
    header = {}
    while line.strip() not in ["", ",,,,,,,,,,,,"]:
        if "," in line:  # for csv from dinoloket
            propval = line.replace('"', "").split(",")
            prop = propval[0]
            prop = prop.replace(":", "").strip()
            val = propval[1]
            if propval[2] != "":
                val = val + " " + propval[2].replace(":", "") + " " + propval[3]
            header[prop] = val
        else:  # for artdiver dino-csv
            propval = line.split(":")
            prop = propval[0]
            val = ":".join(propval[1:])
            header[prop] = val
        line = f.readline()

    return line, header


def _read_empty(f, line):
    while line.strip() in ["", ",,,,,,,,,,,,"] and line != "":
        line = f.readline()
    return line


def _read_dino_groundwater_referencelvl(f, line):
    ref = {}
    while line.strip() not in ["", ",,,,,,,,,,,,"]:
        propval = line.replace('"', "").strip().split(",")
        prop = propval[0]
        prop = prop.replace(":", "")
        if len(propval) > 1:
            val = propval[1]
            ref[prop] = val
        line = f.readline()
    return line, ref


def _read_dino_groundwater_metadata(f, line):
    if "Peildatum" in line:
        return line, {"metadata_available": False}, {}

    _translate_dic_float = {
        "x-coordinaat": "x",
        "y-coordinaat": "y",
        "filternummer": "tube_nr",
    }
    _translate_dic_div_100 = {
        "meetpunt (cm t.o.v. nap)": "tube_top",
        "bovenkant filter (cm t.o.v. nap)": "screen_top",
        "onderkant filter (cm t.o.v. nap)": "screen_bottom",
        "maaiveld (cm t.o.v. nap)": "ground_level",
    }
    metalist = list()
    line = line.strip()
    properties = line.replace('"', "").split(",")

    line = f.readline()

    while line.strip() not in ["", ",,,,,,,,,,,,"]:
        meta = {}
        values = line.replace('"', "").strip().split(",")
        for i, val in enumerate(values):
            meta[properties[i].lower()] = val
        metalist.append(meta)
        line = f.readline()

    meta_ts = {}
    if metalist:
        # add time dependent metadata to meta_ts
        for i, meta in enumerate(metalist):
            meta_tsi = {}
            start_date = pd.to_datetime(meta.pop("startdatum"), dayfirst=True)
            # end_date = pd.to_datetime(meta.pop('einddatum'), dayfirst=True)
            meta_tsi["monitoring_well"] = meta["locatie"]
            for key, item in _translate_dic_float.items():
                if meta[key] == "":
                    meta_tsi[item] = np.nan
                else:
                    meta_tsi[item] = float(meta[key])

            for key, item in _translate_dic_div_100.items():
                if meta[key] == "":
                    meta_tsi[item] = np.nan
                else:
                    meta_tsi[item] = float(meta[key]) / 100.0
            if i == 0:
                for key, item in meta_tsi.items():
                    meta_ts[key] = pd.Series(name=key, dtype=type(item))

            for key, item in meta_tsi.items():
                meta_ts[key].loc[start_date] = item

        # remove series with non time variant metadata from meta_ts
        ts_keys = (
            ["monitoring_well"]
            + list(_translate_dic_float.values())
            + list(_translate_dic_div_100.values())
        )
        for key in ts_keys:
            unique_values = meta_ts[key].unique()
            if len(unique_values) == 1:
                meta_ts.pop(key)

        obs_att = meta_tsi.copy()
        obs_att["name"] = f'{obs_att["monitoring_well"]}-{int(obs_att["tube_nr"]):03d}'
        obs_att["metadata_available"] = True
    else:
        # no metadata
        obs_att = {}
        obs_att["monitoring_well"] = ""
        obs_att["tube_nr"] = np.nan
        obs_att["name"] = "unknown"
        obs_att["x"] = np.nan
        obs_att["y"] = np.nan
        obs_att["tube_top"] = np.nan
        obs_att["ground_level"] = np.nan
        obs_att["screen_top"] = np.nan
        obs_att["screen_bottom"] = np.nan
        obs_att["metadata_available"] = False

    return line, obs_att, meta_ts


def _read_dino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(",")
    while "" in titel:
        titel.remove("")

    if line != "":
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(
                f,
                header=None,
                names=titel,
                parse_dates=["peildatum"],
                index_col="peildatum",
                dayfirst=True,
                usecols=usecols,
            )
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_dino_groundwater_quality_txt(f: Union[str, Path, FileIO]):
    """Read dino groundwater quality (grondwatersamenstelling) from a dinoloket
    txt file.

    Notes
    -----
    this function has not been tested thoroughly

    Parameters
    ----------
    filepath_or_buffer : str
        path to txt file

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """

    if isinstance(f, str):
        fname = f
    else:
        fname = f.name.split(os.sep)[-1]
    if isinstance(f, (str, Path)):
        if isinstance(f, str):
            f = Path(f)
        fname = str(f.stem)
        f = f.open("r")

    logger.info("reading -> {}".format(fname))

    # LOCATIE gegevens
    line = f.readline().rstrip("\n")
    assert line == "LOCATIE gegevens"

    strt_locatie = f.tell()
    # determine the number of rows
    nrows = -1  # the header does not count
    line = f.readline()
    while line not in ["\n", ""]:
        nrows += 1
        line = f.readline()
    eind_locatie = f.tell()
    # go back to where we were before
    f.seek(strt_locatie)
    # read the location-information
    locatie = pd.read_csv(f, sep="\t", nrows=nrows)
    # there is always only one location (change if this is not the case)
    assert nrows == 1
    locatie = locatie.squeeze()

    # KWALITEIT gegevens
    f.seek(eind_locatie)
    read_quality = False
    while line:
        if line.startswith("KWALITEIT gegevens VAST\n"):
            logger.warning("ignoring 'KWALITEIT gegevens VAST'! ")

        if line.startswith("KWALITEIT gegevens VLOEIBAAR"):
            read_quality = True
            break
        line = f.readline()

    if read_quality:
        strt_locatie = f.tell()
        nrows = -1
        while line not in ["\n", ""]:
            nrows += 1
            line = f.readline()
        eind_locatie = f.tell()
        f.seek(strt_locatie)

        measurements = pd.read_csv(
            f,
            sep="\t",
            parse_dates=["Monster datum", "Analyse datum"],
            dayfirst=True,
            index_col="Monster datum",
            nrows=nrows,
        )
    else:
        measurements = pd.Series()

    f.close()

    meta = {
        "filename": fname,
        "source": "dino",
        "monitoring_well": locatie["NITG-nr"],
        "name": locatie["NITG-nr"],
        "x": locatie["X-coord"],
        "y": locatie["Y-coord"],
    }
    try:
        meta["ground_level"] = locatie["Maaiveldhoogte (m tov NAP)"]
    except KeyError:
        meta["ground_level"] = np.nan

    return measurements, meta


def read_dino_groundwater_csv(
    f: Union[str, Path, FileIO],
    to_mnap: bool = True,
    read_series: bool = True,
    remove_duplicates: bool = False,
    keep_dup: str = "last",
):
    """Read dino groundwater quantity data from a dinoloket csv file.

    Parameters
    ----------
    f : Union[str, Path, TextIOWrapper]
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    remove_duplicates : boolean, optional
        if True duplicate indices are removed. Default is False.
    keep_dup : str, optional
        indicate which duplicate indices should be kept, only used when
        remove_duplicates is True. Default is 'last'

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """

    if isinstance(f, (str, Path)):
        if isinstance(f, str):
            f = Path(f)
        fname = f.stem
        f = f.open("r")
    elif isinstance(f, TextIOWrapper):
        fname = f.name
    else:
        raise TypeError("f should be of type str, Path or TextIOWrapper")

    logger.info("reading -> {}".format(fname))

    # read header
    line, header = _read_dino_groundwater_header(f)
    line = _read_empty(f, line)
    if not header:
        logger.warning(f"could not read header -> {fname}")

    # read reference level
    line, ref = _read_dino_groundwater_referencelvl(f, line)
    line = _read_empty(f, line)
    if not ref:
        logger.warning(f"could not read reference level -> {fname}")

    # read metadata
    line, meta, meta_ts = _read_dino_groundwater_metadata(f, line)
    line = _read_empty(f, line)
    if not meta["metadata_available"]:
        logger.warning(f"could not read metadata -> {fname}")
    meta["filename"] = fname
    meta["source"] = "dino"

    # read measurements
    if read_series:
        line, measurements = _read_dino_groundwater_measurements(f, line)
        if measurements is None:
            logger.warning(f"could not read measurements -> {fname}")
        elif measurements[~measurements.stand_cm_tov_nap.isna()].empty:
            logger.warning(f"no NAP measurements available -> {fname}")
        if to_mnap and measurements is not None:
            measurements.insert(
                0, "stand_m_tov_nap", measurements["stand_cm_tov_nap"] / 100.0
            )
            meta["unit"] = "m NAP"
        elif not to_mnap:
            meta["unit"] = "cm NAP"
        if remove_duplicates and measurements is not None:
            measurements = measurements[~measurements.index.duplicated(keep=keep_dup)]

        # add time variant metadata to measurements
        for s in meta_ts.values():
            if measurements is None:
                measurements = pd.DataFrame(data=s.copy(), columns=[s.name])
            else:
                measurements = measurements.join(s, how="outer")
            measurements.loc[:, s.name] = measurements.loc[:, s.name].ffill()

    else:
        measurements = None

    f.close()

    return measurements, meta


def _read_artdino_groundwater_metadata(f, line):
    metalist = list()
    line = line.strip()
    properties = line.split(",")
    line = f.readline()
    while line not in ["\n", "", "\r\n"]:
        meta = dict()
        line = line.strip()
        values = line.split(",")
        for i, val in enumerate(values):
            meta[properties[i].lower()] = val
        metalist.append(meta)
        line = f.readline()

    meta = {}
    if metalist:
        meta["monitoring_well"] = metalist[-1]["locatie"]
        meta["tube_nr"] = int(float(metalist[-1]["filternummer"]))
        meta["name"] = "-".join([meta["monitoring_well"], metalist[-1]["filternummer"]])
        meta["x"] = float(metalist[-1]["x-coordinaat"])
        meta["y"] = float(metalist[-1]["y-coordinaat"])
        meetpunt = metalist[-1]["meetpunt nap"]
        if meetpunt == "":
            meta["tube_top"] = np.nan
        else:
            meta["tube_top"] = float(meetpunt) / 100.0
        maaiveld = metalist[-1]["maaiveld nap"]
        if maaiveld == "":
            meta["ground_level"] = np.nan
        else:
            meta["ground_level"] = float(maaiveld) / 100
        bovenkant_filter = metalist[-1]["bovenkant filter"]
        if bovenkant_filter == "":
            meta["screen_top"] = np.nan
        else:
            meta["screen_top"] = float(bovenkant_filter) / 100
        onderkant_filter = metalist[-1]["onderkant filter"]
        if onderkant_filter == "":
            meta["screen_bottom"] = np.nan
        else:
            meta["screen_bottom"] = float(onderkant_filter) / 100
        meta["metadata_available"] = True
    else:
        # no metadata
        meta["monitoring_well"] = ""
        meta["tube_nr"] = np.nan
        meta["name"] = "unknown"
        meta["x"] = np.nan
        meta["y"] = np.nan
        meta["tube_top"] = np.nan
        meta["ground_level"] = np.nan
        meta["screen_top"] = np.nan
        meta["screen_bottom"] = np.nan
        meta["metadata_available"] = False

    return line, meta


def _read_artdino_groundwater_measurements(f, line):
    line = line.strip()
    titel = line.split(",")
    while "" in titel:
        titel.remove("")

    if line != "":
        # Validate if titles are valid names
        validator = np.lib._iotools.NameValidator()
        titel = [s.lower() for s in validator(titel)]
        usecols = range(0, len(titel))

        try:
            measurements = pd.read_csv(
                f,
                header=None,
                names=titel,
                parse_dates=["peil_datum_tijd"],
                index_col="peil_datum_tijd",
                dayfirst=True,
                usecols=usecols,
            )
        except pd.errors.ParserError:
            # for now the workflow is to remove the files that cannot be read
            # manually.
            measurements = None
    else:
        measurements = None

    return line, measurements


def read_artdino_groundwater_csv(path, to_mnap=True, read_series=True):
    """Read dino groundwater quantity data from a CSV file as exported by
    ArtDiver.

    Parameters
    ----------
    path : str
        path to csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True

    Returns
    -------
    measurements : pd.DataFrame
    meta : dict
        dictionary with metadata
    """
    logger.info(f"reading -> {os.path.split(path)[-1]}")

    with open(path, "r") as f:
        # read header
        line, header = _read_dino_groundwater_header(f)
        line = _read_empty(f, line)
        if not header:
            logger.warning(f"could not read header -> {path}")

        # read metadata
        line, meta = _read_artdino_groundwater_metadata(f, line)
        line = _read_empty(f, line)

        if not meta["metadata_available"]:
            logger.warning(f"could not read metadata -> {path}")

        meta["filename"] = path
        meta["source"] = "dino"

        # read measurements
        if read_series:
            line, measurements = _read_artdino_groundwater_measurements(f, line)
            if measurements is None:
                logger.warning(f"could not read measurements -> {path}")
            elif measurements[~measurements["stand_cm_nap"].isna()].empty:
                logger.warning(f"no NAP measurements available -> {path}")
            if to_mnap and measurements is not None:
                measurements["stand_m_tov_nap"] = measurements["stand_cm_nap"] / 100.0
                meta["unit"] = "m NAP"
            elif not to_mnap:
                meta["unit"] = "cm NAP"
        else:
            measurements = None

    return measurements, meta


def read_artdino_dir(
    dirname,
    ObsClass=None,
    subdir="csv",
    suffix="1.csv",
    unpackdir=None,
    force_unpack=False,
    preserve_datetime=False,
    keep_all_obs=True,
    **kwargs,
):
    """Read Dino directory with point observations.

    TODO:
    - Evt. nog verbeteren door meteen Dataframe te vullen op het moment dat een
    observatie wordt ingelezen. Nu wordt eerst alles ingelezen in een lijst en
    daar een dataframe van gemaakt.

    Parameters
    ----------
    dirname : str
        directory name, can be a .zip file or the parent directory of subdir
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    subdir : str
        subdirectory of dirname with data files
    suffix : str
        suffix of files in subdir that will be read
    unpackdir : str
        destination directory of the unzipped file
    force_unpack : boolean, optional
        force unpack if dst already exists
    preserve_datetime : boolean, optional
        use date of the zipfile for the destination file
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_artdino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations
    """

    from ..util import unzip_file

    # unzip dir
    if dirname.endswith(".zip"):
        zipf = dirname
        if unpackdir is None:
            dirname = tempfile.TemporaryDirectory().name
        else:
            dirname = unpackdir
        unzip_file(
            zipf, dirname, force=force_unpack, preserve_datetime=preserve_datetime
        )

    # read filenames
    files = os.listdir(os.path.join(dirname, subdir))
    if suffix:
        files = [file for file in files if file.endswith(suffix)]

    if not files:
        raise FileNotFoundError(
            "no files were found in {} that end with {}".format(
                os.path.join(dirname, subdir), suffix
            )
        )

    # read individual files
    obs_list = []
    for _, file in enumerate(files):
        path = os.path.join(dirname, subdir, file)
        obs = ObsClass.from_artdino_file(path=path, **kwargs)
        if obs.metadata_available and (not obs.empty):
            obs_list.append(obs)
        elif keep_all_obs:
            obs_list.append(obs)
        else:
            logging.info(f"not added to collection -> {path}")

    return obs_list


def _read_dino_waterlvl_metadata(f, line):
    """read dino waterlevel metadata

    Parameters
    ----------
    f : text wrapper
    line : str
        line with meta dictionary keys
    meta_dic : dict (optional)
        dictionary with metadata

    Returns
    -------
    meta : dict
        dictionary with metadata

    """
    meta_keys = line.strip().split(",")
    meta_values = f.readline().strip().split(",")
    meta = {}
    for key, value in zip(meta_keys, meta_values):
        key = key.strip()
        if key in ["X-coordinaat", "Y-coordinaat"]:
            if key == "X-coordinaat":
                meta["x"] = float(value)
            elif key == "Y-coordinaat":
                meta["y"] = float(value)
        elif key == "Locatie":
            meta["monitoring_well"] = value
            meta["name"] = value

    return meta


def _read_dino_waterlvl_measurements(f, line):
    """

    Parameters
    ----------
    f : text wrapper
    line: str
        header of dataframe

    Returns
    -------
    measurements : pd.DataFrame
    """

    titel = line.strip().split(",")
    while "" in titel:
        titel.remove("")

    validator = np.lib._iotools.NameValidator()
    titel = [i.lower() for i in validator(titel)]
    usecols = range(0, len(titel))

    measurements = pd.read_csv(
        f,
        header=None,
        names=titel,
        parse_dates=["peildatum"],
        index_col="peildatum",
        dayfirst=True,
        usecols=usecols,
    )

    # measurements['Stand (m t.o.v. NAP)'] = measurements['Stand (cm t.o.v. NAP)'] /100.

    # measurements.set_index('Peildatum', inplace=True)

    return measurements


def read_dino_waterlvl_csv(
    f: Union[str, Path, FileIO], to_mnap: bool = True, read_series: bool = True
):
    """Read dino waterlevel data from a dinoloket csv file.

    Parameters
    ----------
    f : str, Path, FileIO
        path to dino water level csv file
    to_mnap : boolean, optional
        if True a column with 'stand_m_tov_nap' is added to the dataframe
    read_series : boolean, optional
        if False only metadata is read, default is True
    """

    fname = ""
    if isinstance(f, (str, Path)):
        if isinstance(f, str):
            f = Path(f)
        fname = f.stem
        f = f.open("r")

    logger.info("reading -> {}".format(fname))

    p_meta = re.compile(
        "Locatie,Externe aanduiding,X-coordinaat,Y-coordinaat, Startdatum, Einddatum"
    )
    p_data = re.compile(r"Locatie,Peildatum,Stand \(cm t.o.v. NAP\),Bijzonderheid")

    line = f.readline()
    while line != "":
        line = f.readline()
        if p_meta.match(line):
            meta = _read_dino_waterlvl_metadata(f, line)
            if meta:
                meta["metadata_available"] = True
            else:
                meta["metadata_available"] = False
            meta["filename"] = fname
            meta["source"] = "dino"
        elif p_data.match(line):
            if read_series:
                measurements = _read_dino_waterlvl_measurements(f, line)
                if to_mnap and measurements is not None:
                    measurements["stand_m_tov_nap"] = (
                        measurements["stand_cm_tov_nap"] / 100.0
                    )
                    meta["unit"] = "m NAP"
                elif not to_mnap:
                    meta["unit"] = "cm NAP"
            else:
                measurements = None

    f.close()

    return measurements, meta


def read_dino_dir(
    path: Union[str, Path],
    ObsClass,
    subdir: str = "DINO_Grondwaterstanden",
    suffix: str = None,
    keep_all_obs: bool = True,
    **kwargs: dict,
):
    """Read Dino directory with point observations.

    TODO:
    - Evt. nog verbeteren door meteen Dataframe te vullen op het moment dat een
    observatie wordt ingelezen. Nu wordt eerst alles ingelezen in een lijst en
    daar een dataframe van gemaakt.
    - aparte unzip functie maken en toch de juiste tijdelijke directory krijgen.

    Parameters
    ----------
    path : str | Path
        directory name, can be a .zip file or the parent directory of subdir
    ObsClass : type
        class of the observations, e.g. GroundwaterObs or WaterlvlObs
    subdir : str
        subdirectory of dirname with data files. For old school dino zip files this
        should be "Grondwaterstanden_Put". For new style the default value
        'DINO_Grondwaterstanden' is sufficient. The default is 'DINO_Grondwaterstanden'.
    suffix : str or None, optional
        suffix of files in subdir that will be read. For old school dino zip files
        this should be '1.csv'. For new style the default value None is sufficient.
        The default is None.
    keep_all_obs : boolean, optional
        add all observation points to the collection, even without data or
        metadata
    **kwargs: dict, optional
        Extra arguments are passed to ObsClass.from_dino_file()

    Returns
    -------
    obs_df : pd.DataFrame
        collection of multiple point observations
    """

    path = Path(path)

    obs_list = []

    def get_dino_obs(f: Union[str, FileIO]):
        obs = ObsClass.from_dino(f, **kwargs)
        if obs.metadata_available and (not obs.empty):
            return obs
        elif keep_all_obs:
            return obs
        else:
            logging.info(f"not added to collection -> {f.name}")
            return None

    if path.suffix == ".zip":
        with ZipFile(path) as zfile:
            fnames = [x for x in zfile.namelist() if f"{subdir}/" in x]
            if suffix:
                if "0" in suffix:
                    raise Exception(f"Cant read dino files with _{suffix}")
                fnames = [x for x in fnames if suffix in x]
            if len(fnames) == 0:
                raise FileNotFoundError(
                    f"no files were found in {path}, {subdir=} with {suffix=}"
                )
            for fname in fnames:
                with zfile.open(fname) as fo:
                    obs = get_dino_obs(TextIOWrapper(fo))
                    if obs is not None:
                        obs_list.append(obs)
    elif path.is_dir():
        subpath = path / subdir
        if suffix:
            if "0" in suffix:
                raise Exception(f"Cant read dino files with _{suffix}")
            elif "*" not in suffix:
                suffix = f"*{suffix}"
            files = list(subpath.glob(suffix))
        else:
            files = list(subpath.iterdir())

        if len(files) == 0:
            raise FileNotFoundError(f"no files were found in {subpath} with {suffix=}")

        for file in files:
            obs = get_dino_obs(file)
            if obs is not None:
                obs_list.append(obs)
    else:
        raise ValueError("Path must either be a .zip or directory")

    return obs_list
