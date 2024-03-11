# %%
import logging

from hydropandas import observation
from hydropandas.io import solinst

logging.basicConfig(level=logging.DEBUG)


# %% test observations


def test_read_solinst_file_obs():
    # observation of slugtest, 8 observations per second
    df, _ = solinst.read_solinst_file(
        "./tests/data/2024-solinst-test/WsNoo_dp366_BUB_20231222_slug1m.xle",
    )

    assert len(df) == 587, "Dataframe should have 587 readings"


def test_read_solinst_file_meta_has_location():
    # observation of slugtest, created via iOs app, location included by app
    _, meta = solinst.read_solinst_file(
        "./tests/data/2024-solinst-test/WsNoo_dp366_BUB_20231222_slug1m.xle",
    )

    assert meta["x"] == 52730.58, "x coordinate should be 52730.58"


def test_read_solinst_file_meta_without_location():
    # example observation created via desktop, location not included
    _, meta = solinst.read_solinst_file(
        "./tests/data/2024-solinst-test/example-10min-interval-via-laptop.zip",
    )

    assert meta["x"] is None, "x coordinate not available"


def test_read_solinst_file_with_manual_meta():
    # manual metadata about levels provided
    screen_bottom = -10
    screen_top = -5
    ground_level = -1
    tube_nr = 1
    tube_top = -0.5

    oc = observation.GroundwaterObs.from_solinst(
        "./tests/data/2024-solinst-test/WsNoo_dp366_BUB_20231222_slug1m.xle",
        screen_bottom=screen_bottom,
        screen_top=screen_top,
        ground_level=ground_level,
        tube_nr=tube_nr,
        tube_top=tube_top,
    )

    assert oc.tube_top == tube_top, "tube_top should be in oc"
