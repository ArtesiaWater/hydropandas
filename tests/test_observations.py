import numpy as np
import os

from art_tools.observations import observation as obs
from art_tools.observations import obs_collection as oc


# %% single observation
fname = r'..\data\2019-Dino-test\Grondwatersamenstellingen_Put\B52C0057.txt'
ogq = obs.GroundwaterQualityObs.from_dino_file(fname, verbose=True)

fname = r'..\data\2019-Dino-test\Peilschaal\P58A0001.csv'
wl = obs.WaterlvlObs.from_dino_file(fname, verbose=True)

fname = r'..\data\2019-Dino-test\Grondwaterstanden_Put\B33F0080001_1.csv'
gw = obs.GroundwaterObs.from_dino_file(fname=fname, verbose=True)

# download dino
gw2 = obs.GroundwaterObs.from_dino_server(name="B57F0077", filternr=4.,
                                          tmin="2000-01-01",
                                          tmax="2010-01-01", unit="NAP")


plot_dir = r"..\data\2019-Dino-test\plots"

gw.to_interactive_plot(savedir=plot_dir, plot_columns=['Stand_m_tov_NAP'],
                       hoover_date_format="{%F}",
                       add_filter_to_legend=True)


# %% collection of observations

fl = oc.ObsCollection.from_fieldlogger(
    r'..\data\2019-Dino-test\fieldlogger\locations.csv')

plot_dir = r"..\data\2019-Dino-test\plots"

# %% read dino directories


dinozip = r'..\data\2019-Dino-test\Dino.zip'


# groundwater quantity
dino_gw = oc.ObsCollection.from_dino_dir(dirname=dinozip, ObsClass=obs.GroundwaterObs,
                                         subdir='Grondwaterstanden_Put',
                                         suffix='1.csv', keep_all_obs=False,
                                         verbose=False)

# surface water
dino_ps = oc.ObsCollection.from_dino_dir(dirname=dinozip, ObsClass=obs.WaterlvlObs,
                                         subdir='Peilschaal', suffix='.csv',
                                         verbose=True)

# groundwater quality
dino_gwq = oc.ObsCollection.from_dino_dir(dirname=dinozip, ObsClass=obs.GroundwaterQualityObs,
                                          subdir='Grondwatersamenstellingen_Put',
                                          suffix='.txt', verbose=True)


# %% downlaod
extent = [120300, 120500, 439000, 441000]  # Schoonhoven zoomed

dino_gw_extent = oc.ObsCollection.from_dino_server(extent=extent,
                                                   ObsClass=obs.GroundwaterObs,
                                                   verbose=True)
bbox = [120300, 439000, 120500, 441000]  # Schoonhoven zoomed
bbox = np.array([191608.334, 409880.402, 193072.317, 411477.894])

dino_gw_bbox = oc.ObsCollection.from_dino_server(bbox=bbox,
                                                 ObsClass=obs.GroundwaterObs,
                                                 verbose=True)

# %% collection methods

dino_gw[['nearest point', 'distance nearest point']
        ] = dino_gw.get_nearest_point(fl)

fdf = dino_gw.to_fieldlogger(
    r'..\data\2019-Dino-test\fieldlogger\locations.csv', verbose=True)

dino_gw.get_lat_lon()
dino_gw.to_interactive_map(plot_dir, plot_columns=['Stand_m_tov_NAP'],
                           fname='imap.html',
                           legend_name='grondwater DINO',
                           add_filter_to_legend=True, hoover_names=['gws'],
                           zoom_start=9,
                           verbose=True)

# %% read FEWS data

plot_dir = r"..\data\2019-FEWS-test\plots"
fname = 'texel_fews.html'

fews_gw_prod = oc.ObsCollection.from_fews(r'..\data\2019-FEWS-test\WaalenBurg_201810-20190215_prod.zip',
                                          verbose=True)

ax = fews_gw_prod.to_map()

fews_gw_prod.to_interactive_map(plot_dir, plot_columns=['value'], fname=fname,
                                plot_freq='D', legend_name='opp water FEWS',
                                map_label='locationId', map_label_size=10)
