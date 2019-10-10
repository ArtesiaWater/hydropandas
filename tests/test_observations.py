import os
from observations import observation as obs
from observations import obs_collection as oc
import numpy as np
import sys
sys.path.insert(1, "..")


TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
sys.path.insert(0, PROJECT_DIR)
os.chdir(TEST_DIR)


plot_dir = r".\data\2019-Dino-test\plots"
dinozip = r'.\data\2019-Dino-test\Dino.zip'


def test_observation_gwq():
    # single observation
    fname = r'.\data\2019-Dino-test\Grondwatersamenstellingen_Put\B52C0057.txt'
    ogq = obs.GroundwaterQualityObs.from_dino_file(fname, verbose=True)
    return ogq


def test_observation_wl():
    fname = r'.\data\2019-Dino-test\Peilschaal\P58A0001.csv'
    wl = obs.WaterlvlObs.from_dino_file(fname, verbose=True)
    return wl


def test_observation_gw():
    fname = r'.\data\2019-Dino-test\Grondwaterstanden_Put\B33F0080001_1.csv'
    gw = obs.GroundwaterObs.from_dino_file(fname=fname, verbose=True)
    return gw


def test_observation_dino_download():
    # download dino
    gw2 = obs.GroundwaterObs.from_dino_server(name="B57F0077", filternr=4.,
                                              tmin="2000-01-01",
                                              tmax="2010-01-01", unit="NAP")
    return gw2


def test_observation_dino_download2():
    # download dino
    gw2 = obs.GroundwaterObs.from_dino_server(name="B57B0069", filternr=1.,
                                              tmin="2000-01-01",
                                              tmax="2010-01-01", unit="NAP")
    return gw2


def test_observation_dino_download3():
    # download dino data from pb without extra metadata. For this pb
    # art.dino_wfs.get_dino_piezometer_metadata() returns an empty list
    gw3 = obs.GroundwaterObs.from_dino_server(name="B45G1147", filternr=1.,
                                              tmin="1900-01-01",
                                              tmax="2020-01-01", unit="NAP")
    return gw3


def test_interactive_plot():
    gw = test_observation_gw()
    gw.to_interactive_plot(savedir=plot_dir, plot_columns=['Stand_m_tov_NAP'],
                           hoover_date_format="{%F}",
                           add_filter_to_legend=True)
    return


def test_obscollection_fieldlogger():
    # collection of observations
    fl = oc.ObsCollection.from_fieldlogger(
        r'.\data\2019-Dino-test\fieldlogger\locations.csv')
    return fl

# read dino directories


def test_obscollection_dinozip_gw():
    # groundwater quantity
    dino_gw = oc.ObsCollection.from_dino_dir(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        keep_all_obs=False,
        verbose=False)
    return dino_gw


def test_obscollection_dinozip_gw_keep_all_obs():
    # do not delete empty dataframes
    dino_gw = oc.ObsCollection.from_dino_dir(
        dirname=dinozip,
        ObsClass=obs.GroundwaterObs,
        subdir='Grondwaterstanden_Put',
        suffix='1.csv',
        keep_all_obs=True,
        verbose=False)
    return dino_gw


def test_obscollection_dinozip_wl():
    # surface water
    dino_ps = oc.ObsCollection.from_dino_dir(
        dirname=dinozip,
        ObsClass=obs.WaterlvlObs,
        subdir='Peilschaal',
        suffix='.csv',
        verbose=True)
    return dino_ps


def test_obscollection_dinozip_gwq():
    # groundwater quality
    dino_gwq = oc.ObsCollection.from_dino_dir(
        dirname=dinozip,
        ObsClass=obs.GroundwaterQualityObs,
        subdir='Grondwatersamenstellingen_Put',
        suffix='.txt',
        verbose=True)
    return dino_gwq


def test_obscollection_dino_download_extent():
    # download DINO from extent
    extent = [120300, 120500, 439000, 441000]  # Schoonhoven zoomed
    dino_gw_extent = oc.ObsCollection.from_dino_server(
        extent=extent, ObsClass=obs.GroundwaterObs, verbose=True)
    return dino_gw_extent


def test_obscollection_dino_download_bbox():
    # download DINO from bbox
    bbox = [120300, 439000, 120500, 441000]  # Schoonhoven zoomed
    bbox = np.array([191608.334, 409880.402, 193072.317, 411477.894])

    dino_gw_bbox = oc.ObsCollection.from_dino_server(
        bbox=bbox, ObsClass=obs.GroundwaterObs, verbose=True)
    return dino_gw_bbox


def test_obscollection_dino_download_bbox_empty():
    # download DINO from bbox
    bbox = [88596.63500000164, 407224.8449999988,
            89623.4149999991, 407804.27800000086]

    dino_gw_bbox = oc.ObsCollection.from_dino_server(
        bbox=bbox, ObsClass=obs.GroundwaterObs, verbose=True)
    return dino_gw_bbox


# collection methods


def test_get_nearest_point():
    dino_gw = test_obscollection_dinozip_gw()
    fl = test_obscollection_fieldlogger()
    dino_gw[['nearest point', 'distance nearest point']
            ] = dino_gw.get_nearest_point(fl)
    return dino_gw


def test_get_filternr():
    dino_gw = test_obscollection_dinozip_gw()
    dino_gw.get_filternr(if_exists='replace')
    return dino_gw


def test_obscollection_to_fieldlogger():
    dino_gw = test_obscollection_dinozip_gw()
    fdf = dino_gw.to_fieldlogger(
        r'.\data\2019-Dino-test\fieldlogger\locations.csv', verbose=True)

    return fdf


def test_within_extent():
    dino_gw = test_obscollection_dinozip_gw()
    extent = [210350, 213300, 473300, 474000]
    dino_gw.within_extent(extent, inplace=True)
    assert dino_gw.shape[0] == 4
    return dino_gw


def test_obscollection_dino_to_map():
    dino_gw = test_obscollection_dinozip_gw()
    dino_gw.get_lat_lon()
    dino_gw.to_interactive_map(plot_dir, plot_columns=['Stand_m_tov_NAP'],
                               fname='imap.html',
                               legend_name='grondwater DINO',
                               add_filter_to_legend=True, hoover_names=['gws'],
                               zoom_start=9,
                               verbose=True)
    return


def test_obscollection_dino_to_mapgraph():
    gw = test_obscollection_dinozip_gw()
    gw.to_mapgraphs(plot_ylim='min_dy')

    return


def test_obscollection_consecutive_obs_years():
    gw = test_obscollection_dinozip_gw_keep_all_obs()
    coy = gw.consecutive_obs_years()

    return coy


def test_obscollection_get_seasonal_stats():
    gw = test_obscollection_dinozip_gw_keep_all_obs()
    st = gw.get_seasonal_stat(stat='mean')

    return st

# read FEWS data


def test_obscollection_fews():
    fews_gw_prod = oc.ObsCollection.from_fews(
        r'.\data\2019-FEWS-test\WaalenBurg_201810-20190215_prod.zip',
        verbose=True,
        to_mnap=False,
        remove_nan=False)
    return fews_gw_prod


def test_obscollection_fews2():
    fews_gw_prod = oc.ObsCollection.from_fews2(
        r'.\data\2019-FEWS-test\WaalenBurg_201810-20190215_prod.zip',
        verbose=True,
        locations=None)
    return fews_gw_prod


def test_obscollection_fews2_selection():
    fews_gw_prod = oc.ObsCollection.from_fews2(
        r'.\data\2019-FEWS-test\WaalenBurg_201810-20190215_prod.zip',
        verbose=True,
        locations=(
            "MPN-N-2",
        ))
    return fews_gw_prod


def test_obscollection_to_map():
    fname = 'texel_fews.html'
    plot_dir = r".\data\2019-FEWS-test\plots"
    fews_gw_prod = test_obscollection_fews()
    ax = fews_gw_prod.to_map()
    fews_gw_prod.to_interactive_map(
        plot_dir,
        plot_columns=['value'],
        fname=fname,
        plot_freq='D',
        legend_name='opp water FEWS',
        map_label='locationId',
        map_label_size=10)
    return ax


# read WISKI data

def test_observation_wiskicsv_gw():
    wiski_gw = obs.GroundwaterObs.from_wiski(
        r".\data\2019-WISKI-test\1016_PBF.csv",
        sep=r'\s+',
        header_sep=':',
        header_identifier=':',
        parse_dates={
            "datetime": [
                0,
                1]},
        index_col=["datetime"],
        translate_dic={
            'name': 'Station Number',
            'x': 'GlobalX',
            'y': 'GlobalY'},
        verbose=True)

    return wiski_gw


def test_obscollection_wiskizip_gw():
    wiski_col = oc.ObsCollection.from_wiski(
        r".\data\2019-WISKI-test\1016_PBF.zip",
        translate_dic={
            'name': 'Station Number',
            'x': 'GlobalX',
            'y': 'GlobalY'},
        sep=r'\s+',
        header_sep=':',
        dayfirst=True,
        header_identifier=':',
        parse_dates={
            "datetime": [
                0,
                1]},
        index_col=["datetime"],
        verbose=True)

    return wiski_col


# Test Pystore

def test_obscollection_to_pystore():
    obsc = test_obscollection_fews()
    obsc.to_pystore("test_pystore", "./data/2019-Pystore-test",
                    groupby="locationId", overwrite=True)


def test_obscollection_from_pystore():
    obsc = oc.ObsCollection.from_pystore(
        "test_pystore", "./data/2019-Pystore-test")
    return obsc


def test_obscollection_pystore_only_metadata():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./data/2019-Pystore-test",
                                         read_series=False)
    return obsc


def test_obscollection_pystore_extent():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./data/2019-Pystore-test",
                                         extent=[115534, 115539, 0, 10000000]
                                         )
    return obsc


def test_obscollection_pystore_item_names():
    obsc = oc.ObsCollection.from_pystore("test_pystore",
                                         "./data/2019-Pystore-test",
                                         item_names=['MPN-N-2']
                                         )
    return obsc


def test_obs_from_pystore_item():
    import pystore
    pystore.set_path("./data/2019-Pystore-test")
    store = pystore.store("test_pystore")
    coll = store.collection(store.collections[0])
    item = coll.item(list(coll.list_items())[0])
    o = obs.GroundwaterObs.from_pystore_item(item)
    return o

# Test KNMI Obs
def test_knmi_obs_from_stn():
    return obs.KnmiObs.from_knmi(829, "RD")


def test_knmi_obs_from_xy():
    return obs.KnmiObs.from_nearest_xy(100000, 350000, "RD")


def test_knmi_obs_from_obs():
    pb = test_observation_gw()
    return obs.KnmiObs.from_obs(pb, "EV24")

# Test Menyanthes (still need a small menyanthes file to do the test)

#def test_obscollection_menyanthes():
#            
#    fname = r'g:\My Drive\m\01projekt\19042016 BRABANT WATER, Uitwerking pompproef Gilze\02 Data\export_from_ADI_20191007.men'
#    obsc = oc.ObsCollection.from_menyanthes(fname, verbose=True)
#    
#    return obsc
