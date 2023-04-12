import os
from pathlib import Path

import hydropandas as hpd
from hydropandas.io import fews

fews_testdata = Path(__file__).parent / "data/2019-FEWS-test"

dirname, fnames = hpd.util.get_files(
    str(fews_testdata / "WaalenBurg_201810-20190215_prod.zip"), ".xml"
)


def test_fews_highmemory() -> None:
    fews.read_xml_fname(
        os.path.join(dirname, fnames[0]), low_memory=False, ObsClass=hpd.WaterlvlObs
    )


def test_fews_lowmemory() -> None:
    fews.read_xml_fname(
        os.path.join(dirname, fnames[0]), low_memory=True, ObsClass=hpd.WaterlvlObs
    )


def test_obscollection_fews_selection() -> None:
    fews.read_xml_filelist(fnames, hpd.WaterlvlObs, dirname)


def test_fews_pid_wsvv() -> None:
    wsvv_pid = hpd.get_fews_pid("WSVV")
    oc = hpd.ObsCollection.from_fews_xml(
        str(fews_testdata / "test_wsvv_fews.xml"), ObsClass=wsvv_pid
    )
    assert isinstance(oc.iloc[0].obs, hpd.PrecipitationObs)
    assert isinstance(oc.iloc[1].obs, hpd.EvaporationObs)
