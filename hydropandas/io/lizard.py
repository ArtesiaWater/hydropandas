import pandas as pd
import requests 
from pyproj import Transformer 
import logging
import numpy as np
logger = logging.getLogger(__name__)
from shapely.geometry import Polygon
from tqdm import tqdm
import math
import geopandas
import concurrent.futures
#%% check_status_obs
def check_status_obs (metadata, timeseries):
    '''
    checks if a monitoring well is still active 

    Parameters
    ----------
    metadata : pandas.DataFrame
        metadata of the monitoring well 
    timeseries : pandas DataFrame
        timeseries of the monitoring well 
        
    Returns
    -------
    metadata DataFrame including the status of the monitoring well 
    
    '''
    if timeseries.empty:
        metadata["status"] = 'geen tijdsserie beschikbaar'
        return metadata
        
    last_measurement_date = timeseries.last_valid_index()
    today = pd.to_datetime('today').normalize() 

    if today - last_measurement_date < pd.Timedelta(days= 180):
        metadata["status"] = 'active'

    else:
        metadata["status"] = 'niet actief'
    
    
    return metadata
#%% extent_to_polygon
def extent_to_polygon (coordinates):
    '''
    Translates a list of coordinates (xmin, ymin, xmax, ymax) to a polygon with 
    coordinate system WGS84 

    Parameters
    ----------
    coordinates : lst
        list of the modelextent within which the observations are collected 

    Returns
    -------
    polygon of the modelextent with coordinate system WGS84

    '''
    transformer = Transformer.from_crs("EPSG:28992","WGS84")

    lon_min, lat_min = transformer.transform(coordinates[0], coordinates[1])
    lon_max, lat_max = transformer.transform(coordinates[2], coordinates[3])

    poly_T =  Polygon([
        (lat_min, lon_min), (lat_max, lon_min),(lat_max, lon_max), (lat_min,lon_max)])
    
    return poly_T
#%% translate_flag
def translate_flag (timeseries):
    '''
    Translates Vitens Lizard flags from interter to text

    Parameters
    ----------
    timeseries : pandas.DataFrame 
    timeseries of a monitoring well with flags

    Returns
    -------
    timeseries : pandas.DataFrame 
        timeseries with translated quality flags

    '''
    for idx, flag in enumerate(timeseries.loc[:,"flag"]):
        if flag == 0 or flag ==1:
            timeseries.loc[idx,'flag'] = 'betrouwbaar'
        
        elif flag == 3 or flag ==4:
            timeseries.loc[idx,'flag'] = 'onbeslist'
        
        elif flag == 6 or flag ==7:
            timeseries.loc[idx,'flag'] = 'onbetrouwbaar'
            
        elif flag == 99:
            timeseries.loc[idx,'flag'] = 'onongevalideerd'
            
        elif flag == -99:
            timeseries.loc[idx,'flag'] = 'verwijderd'
    
    return timeseries
#%% get_API_results_from_code
def get_API_results_from_code(code, url_lizard): 
    '''
    extracts the Groundwater Station parameters from a monitoring well based 
    on the code of the monitoring well 
    
    Parameters
    ----------
    code : str
        code of the monitoring well 
    url_lizard : str
        location of the LIZARD-API

    Raises
    ------
    ValueError
        if code of the monitoring well is not known

    Returns
    -------
    groundwaterstation_metadata : dict
        dictionary with all available metadata of the monitoring well and its filters

    '''
    lizard_GWS_endpoint = f'{url_lizard}groundwaterstations/'
    url_groundwaterstation_code = f'{lizard_GWS_endpoint}?code={code}'
    
    try:
        groundwaterstation_metadata = requests.get(url_groundwaterstation_code).json()["results"][0]
 
    except IndexError:
        raise ValueError("Code is invalid")  
        
    return groundwaterstation_metadata
#%%
def _prepare_API_input(nr_pages, url_groundwater):
    '''
    prepare API data pages within the defined extent 

    Parameters
    ----------
    nr_pages : int
        number of the pages on which the information is stored
    url_groundwater : str
        location of the used API to extract the data 

    Returns
    -------
    proces_input : list
        list of the page number and the corresponding url 

    '''
    proces_input = []
    for page in range(nr_pages):
        true_page = page+1 # Het echte paginanummer wordt aan de import thread gekoppeld
        url = url_groundwater+'&page={}'.format(true_page)
        item = [true_page,url]
        proces_input += [item]
    return proces_input

def _download_API(data):
    '''
    Function to download the data from the API using the ThreadPoolExecutor 

    Parameters
    ----------
    data : list
        list of the page number and the corresponding url

    Returns
    -------
    None.

    '''
    page, url = data
    try:
        data = requests.get(url = url)
        # succes += 1
        data = data.json()['results']
    except:
        data = []
    return(data)


#%% get_API_results_from_extent
def get_API_results_from_extent(polygon_extent, url_lizard, page_size = 100, nr_threads = 10): 
    '''
    extracts the Groundwater Station parameters from a monitoring well based 
    on the code 
    
    Parameters
    ----------
    code : str
        code of the monitoring well 
    url_lizard : str
        location of the LIZARD-API

    Raises
    ------
    ValueError
        if code of the monitoring well is not known

    Returns
    -------
    groundwaterstation_metadata :dict
        dictionary with all available metadata of the monitoring well and its filters

    '''
    lizard_GWS_endpoint = f'{url_lizard}groundwaterstations/'
    url_groundwaterstation_extent = f'{lizard_GWS_endpoint}?geometry__within={polygon_extent}&page_size={page_size}'
    
    try:
        groundwaterstation_data = requests.get(url_groundwaterstation_extent).json()
        nr_results = groundwaterstation_data['count']
        nr_pages = math.ceil(nr_results/page_size)
        
        print("Number of monitoring wells: {}".format(nr_results))
        print("Number of pages: {}".format(nr_pages))
            
            
        if nr_threads > nr_pages:
            nr_threads = nr_pages
        
        proces_input = _prepare_API_input(nr_pages, url_groundwaterstation_extent)
        groundwaterstation_results = pd.DataFrame()


        with concurrent.futures.ThreadPoolExecutor(max_workers = nr_threads) as executor:
            for result in tqdm(executor.map(_download_API,proces_input),total = nr_pages, desc="Page"):
                groundwaterstation_results = pd.concat([groundwaterstation_results,pd.DataFrame(result)])
             
    except IndexError:
        raise ValueError("Extent is invalid")  
        
    return nr_results, groundwaterstation_results

#%% get_metadata_filter
def get_metadata_filter(API_result, tube_nr, url_lizard):  
    """
    extract the metadata for a specific location from the dict with all 
    groundwater station metadata
    
    Parameters
    ----------
    API_result : dict
        dictionary with all available metadata of the monitoring well and all its filters
    tube_nr : int, optional
        select metadata from a specific tube number
        Default selects tube_nr = 1
    url_lizard : str
        location of the LIZARD-API
        
    Raises
    ------
    ValueError
        if code of the monitoring well is invalid.
    Returns
    -------
    padndas DataFrame of metadata of the specific monitoring well 
    """
    
    metadata = pd.DataFrame(columns= ["monitoring_well", "tube_nr", "name", "x", "y",
                                      "tube_top", "ground_level","screen_top", "screen_bottom","status",
                                      "timeseries_available", "uuid_hand", "start_hand", "uuid_diver", 
                                      'start_diver', 'source', "unit"])

    lon, lat, _ = API_result['geometry']['coordinates']
    transformer = Transformer.from_crs("WGS84","EPSG:28992")
    x,y = transformer.transform(lat,lon)

    for idx, location in enumerate(API_result["filters"]):
            metadata.loc[idx, 'unit'] = 'm NAP'
            name = location['code'].replace("-","")
            metadata.loc[idx, 'x'] = np.round(x,2)
            metadata.loc[idx, 'y'] = np.round(y,2)
            metadata.loc[idx,"name"] = location["code"]
            metadata.loc[idx, "monitoring_well"] = API_result['name']
            metadata.loc[idx,
                         "ground_level"] = API_result['surface_level']
            
            
            metadata.loc[idx, "tube_nr"] = int(name[-3:])
           
            metadata.loc[idx, "tube_top"] = location['top_level']
            metadata.loc[idx,
                         "screen_top"] = location['filter_top_level']
            metadata.loc[idx,
                         "screen_bottom"] = location['filter_bottom_level']
            metadata.loc[idx,
                         "source"] = 'lizard'

            if not location['timeseries']:
                metadata.loc[idx, "timeseries_available"] = 'Nee'

            else:
                timeseries = location['timeseries']
                # metadata.loc[name, "_wezig"] = 'Ja'
                for series in timeseries:
                    series_info = requests.get(series).json()
                    if series_info["name"] == 'WNS9040.hand':
                        metadata.loc[idx, "uuid_hand"] = series_info["uuid"]
                        metadata.loc[idx,
                                     "start_hand"] = series_info["start"]
                    elif series_info["name"] == 'WNS9040':
                        metadata.loc[idx, "uuid_diver"] = series_info["uuid"]
                        metadata.loc[idx,
                                     "start_diver"] = series_info["start"]

            # geen tijdreeksen aanwezig
            if pd.isna(metadata.loc[idx, "start_diver"]) and pd.isna(metadata.loc[idx, "start_hand"]):
                metadata.loc[idx, "timeseries_available"] = 'Nee'
            elif pd.notna(metadata.loc[idx, "start_diver"]) and pd.isna(metadata.loc[idx, "start_hand"]):
                metadata.loc[idx, "timeseries_available"] = 'Diver'
            elif pd.isna(metadata.loc[idx, "start_diver"]) and pd.notna(metadata.loc[idx, "start_hand"]):
                metadata.loc[idx, "timeseries_available"] = 'handpeilingen'
            elif pd.notna(metadata.loc[idx, "start_diver"]) and pd.notna(metadata.loc[idx, "start_hand"]):
                metadata.loc[idx, "timeseries_available"] = 'Diver + hand'
        
    metadata.sort_values(by = ['tube_nr'], inplace = True, ignore_index = True)
    
    tube_nr = [tube_nr] if isinstance(tube_nr, int) else tube_nr
 
    if tube_nr is None:
        metadata = metadata.loc[metadata["tube_nr"] == 1]
    elif tube_nr is not None and tube_nr != 'all':
            metadata =metadata[metadata['tube_nr'].isin(tube_nr)]
            
    metadata = metadata if isinstance(metadata, pd.DataFrame) else pd.DataFrame(metadata).T
    return metadata
#%% get_timeseries
def get_timeseries (uuid, code, tube_nr, tmin , tmax ,url_lizard, page_size = 100000):
    """
    Get the time series of a specific monitoring well 
    ----------
    uuid : str
        Universally Unique Identifier of the monitoring well.
    code : str
        code or name of the monitoring well  
    tube_nr : int, optional
        select specific tube number
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : int YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    url_lizard : str
        location of the LIZARD-API.
    page_size : int, optional
        Query parameter which can extend the response size. The default is 100000.

    Returns
    -------
    pandas DataFrame with the timeseries of the monitoring well 

    """
    
    url_timeseries = url_lizard+'timeseries/{}'.format(uuid)
    
    if tmin != None:
        tmin = pd.to_datetime(tmin).isoformat('T')
        
    if tmax != None:
        tmax = pd.to_datetime(tmax).isoformat('T')
        
        
    params= {'start':tmin, 'end': tmax, 'page_size': page_size}
    url = url_timeseries + '/events/'  

    time_series_events = requests.get(url=url, params=params).json()['results']
    time_series_df = pd.DataFrame(time_series_events) 
    
    if time_series_df.empty:
        # raise ValueError("{} doesn't have measurements in the selected period between {} and {}".format(code, tmin, tmax))
        # print("{} doesn't have measurements in the selected period between {} and {}".format(code, tmin, tmax))
        return pd.DataFrame()
    
    else: 
        time_series_df = translate_flag(time_series_df)
    
        timeseries_sel = time_series_df.loc[:,['time','value', "flag","comment"]]
        timeseries_sel['time'] = pd.to_datetime(timeseries_sel['time'], format = '%Y-%m-%dT%H:%M:%SZ',
                                             errors = 'coerce') + pd.DateOffset(hours = 1)
    
        timeseries_sel = timeseries_sel[~timeseries_sel['time'].isnull()]
        
        timeseries_sel.set_index('time', inplace = True)
        timeseries_sel["name"] = code 
        timeseries_sel["filter_nr"] = tube_nr
        timeseries_sel.index.rename("peil_datum_tijd", inplace = True)
        timeseries_sel = timeseries_sel.loc[:,['name', 'filter_nr','value','flag']]
        timeseries_sel.dropna(inplace = True)
    
    
    return timeseries_sel
#%% merge_timeseries
def merge_timeseries(hand_measurements, diver_measurements):
    """
    merges the timeseries of the hand and diver measurements into one timeserie

    Parameters
    ----------
    hand_measurements : DataFrame
        DataFrame containing the hand measurements of the monitoring well 
    diver_measurements : DataFrame
        DataFrame containing the Diver measurements of the monitoring well 
    
    Returns
    -------
    DataFrame where hand and diver measurements are merged in one timeseries 

    """
    if hand_measurements.empty and diver_measurements.empty:
        measurements = pd.DataFrame()
        
    elif diver_measurements.first_valid_index() == None:
        measurements = hand_measurements
        print("no diver measuremets available for {}".format(hand_measurements.iloc[0]['name']))
    
    else:
          
    
        hand_measurements_sel = hand_measurements.loc[hand_measurements.index < diver_measurements.first_valid_index()]
        measurements = pd.concat([hand_measurements_sel, diver_measurements], axis = 0)
     
    return measurements
#%% combine_timeseries
def combine_timeseries (hand_measurements, diver_measurements):
    """
    combines the timeseries of the hand and diver measurements into one DataFrame

    Parameters
    ----------
    hand_measurements : DataFrame
        DataFrame containing the hand measurements of the monitoring well 
    diver_measurements : DataFrame
        DataFrame containing the Diver measurements of the monitoring well 

    Returns
    -------
    a combined DataFrame with both hand, and diver measurements
        DESCRIPTION.

    """
    hand_measurements.rename(columns =  {"value": "value_hand",  "flag": "flag_hand"}, inplace = True)
    diver_measurements.rename(columns = {"value": "value_diver", "flag": "flag_diver"}, inplace = True)
    
    measurements = pd.concat([hand_measurements, diver_measurements], axis = 1)
    measurements = measurements.loc[:,[ "value_hand","value_diver", "flag_hand","flag_diver"]]
    measurements.loc[:,"name"] = hand_measurements.loc[:,"name"][0]
    measurements.loc[:,"filter_nr"] = hand_measurements.loc[:,"filter_nr"][0]
    
    return measurements 

#%% extract_timeseries_from_API
def extract_timeseries_from_API (metadata_df, tmin, tmax, type_timeseries, url_lizard):
    '''
    extracts timeseries for a specific monitoring well 

    Parameters
    ----------
    metadata_df : pandas DataFrame
        metadata dataframe of the monitoring wel 
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        type of timeseries to;
            hand: returns only hand measurements
            diver: returns only diver measurements  
            merge: the hand and diver measurements into one time series (merge; default) or
            combine: keeps hand and diver measurements separeted       
        The default is merge.
    url_lizard : str
        location of the LIZARD-API.

    Returns
    -------
    measurements : pandas DataFrame
        dataframe with the timeseries of the monitoring well 
    metadata_df : pandas DataFrame
        dataframe with the metadata of the monitoring well 

    '''
   
    metadata_df = metadata_df.squeeze()
    if metadata_df["timeseries_available"] != "Nee":
   
        if metadata_df["timeseries_available"] == 'Diver + hand':
            hand_measurements =  get_timeseries(metadata_df['uuid_hand'], metadata_df["name"],metadata_df["tube_nr"], tmin, tmax, url_lizard)
            diver_measurements = get_timeseries(metadata_df['uuid_diver'],metadata_df["name"],metadata_df["tube_nr"], tmin, tmax, url_lizard)
        
            if type_timeseries == "hand":
                measurements = hand_measurements
            elif type_timeseries == "diver":
                measurements = diver_measurements
                    
            elif type_timeseries == "merge":
                measurements = merge_timeseries(hand_measurements, diver_measurements)
            elif type_timeseries =="combine":
                measurements = combine_timeseries(hand_measurements, diver_measurements)
            
        # Diver
        elif metadata_df ["timeseries_available"] == 'Diver':
             measurements = get_timeseries(metadata_df['uuid_diver'],metadata_df["name"], metadata_df["tube_nr"],
                                          tmin, tmax, url_lizard)
           
        # HAND
        elif metadata_df["timeseries_available"] == 'handpeilingen':
                measurements = get_timeseries(metadata_df['uuid_hand'],metadata_df["name"], metadata_df["tube_nr"], tmin, tmax, url_lizard)
    
    elif metadata_df["timeseries_available"] == "Nee": 
        measurements = pd.DataFrame()
        
        
        
    metadata_df.drop(['uuid_hand','uuid_diver'],  inplace = True)
        
    return measurements, metadata_df 

#%% read_lizard_groundwater
def read_lizard_groundwater_from_code (code, tube_nr=None,
                            tmin = None, tmax = None,
                            type_timeseries = 'merge',
                            url_lizard='https://vitens.lizard.net/api/v4/'):
    """
    extracts the metadata and timeseries of a observation well from a LIZARD-API based on
    the code of a monitoring well 

    Parameters
    ----------
    code : str
        code of the measuring well  
    tube_nr : int, optional
        select specific tube top
        Default selects tube_nr = 1
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements  
        merge: the hand and diver measurements into one time series (merge; default) or
        combine: keeps hand and diver measurements separeted       
        The default is merge.
    url_lizard : str
        location of the LIZARD-API.

    Returns
    -------
    returns a DataFrame with metadata and timeseries 
    """
        
    groundwaterstation_metadata = get_API_results_from_code(code, url_lizard)
     
    obs_df = get_metadata_filter(groundwaterstation_metadata,tube_nr, url_lizard)
    
    if obs_df.empty:
        raise ValueError("{} doesn't have a tube number {}".format(code,tube_nr))
        
    measurements, obs_df = extract_timeseries_from_API (obs_df, tmin, tmax,type_timeseries, url_lizard)
    obs_df = check_status_obs(obs_df,measurements)


    return measurements,obs_df.to_dict()

#%% get_obs_list_from_code
def get_obs_list_from_code (code, tube_nr='all', 
                            tmin = None, tmax = None,
                            type_timeseries = 'merge',
                            url_lizard='https://vitens.lizard.net/api/v4/'):
    """
    get all observations from a list of codes of the monitoring wells and a 
    list of tube numbers 

    Parameters
    ----------
    code : lst of str 
        codes of the monitoring wells 
    tube_nr : lst of str 
        list of tube numbers of the monitoring wells that should be selected.
        By default 'all' available tubes are selected.
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        hand: returns only hand measurements
        diver: returns only diver measurements  
        merge: the hand and diver measurements into one time series (merge; default) or
        combine: keeps hand and diver measurements separeted       
        The default is merge.
    url_lizard : str
        location of the LIZARD-API.
    Returns
    -------
    ObsCollection
        ObsCollection DataFrame with the 'obs' column

    """
    
    obs_col = pd.DataFrame(columns = ['monitoring_well', 'tube_nr', 'name', 'x', 'y', 'tube_top',
                                      'ground_level', 'screen_top', 'screen_bottom', 'status',
                                      'timeseries_available', 'start_hand', 'start_diver', 'source', 
                                      'unit','obs'])
    if not isinstance(code, list):
        raise ValueError("Code should be a list")
    
    if len(code) == 1:
        code = code[0]
        groundwaterstation_metadata = get_API_results_from_code(code, url_lizard)
        obs_df = get_metadata_filter(groundwaterstation_metadata,tube_nr, url_lizard)

        for idx, row in obs_df.iterrows():
            measurements, obs_series = extract_timeseries_from_API (row, tmin, tmax,type_timeseries, url_lizard)
            obs_series['obs'] = measurements.squeeze()
            obs_series = check_status_obs(obs_series, measurements)
            obs_col.loc[idx] = obs_series                
            
                
    else:
        for elem in code: 
            groundwaterstation_metadata = get_API_results_from_code(elem, url_lizard)                 
            obs_df = get_metadata_filter(groundwaterstation_metadata,tube_nr, url_lizard)
            obs_col = pd.concat([obs_df,obs_col], axis = 0, ignore_index = True)
            
            for idx, row in obs_df.iterrows():
                measurements, obs_series = extract_timeseries_from_API (row, tmin, tmax,type_timeseries, url_lizard)
                obs_series['obs'] = measurements.squeeze()
                obs_series = check_status_obs(obs_series, measurements)
                obs_col.loc[idx] = obs_series                
            
    return obs_col


#%% get_obs_list_from_extent 
                       
def get_obs_list_from_extent(extent, extract_timeseries = True,
                            tmin = None, tmax = None,
                            type_timeseries = 'merge',
                            url_lizard='https://vitens.lizard.net/api/v4/'):
    '''
    get all observations within a specified extent 
    Parameters
    ----------
    extent : list or a shapefile 
        get groundwater monitoring wells wihtin this extent [xmin, ymin, xmax, ymax]
        or within a predefined Polygon from a shapefile 
    extract_timeseries : Bool, optional
        Extract timeseries or not, if not only metadata are returned 
    tmin : str YYYY-m-d, optional
        start of the observations, by default the entire serie is returned
    tmax : Ttr YYYY-m-d, optional
        end of the observations, by default the entire serie is returned
    type_timeseries : str, optional
        merge: the hand and diver measurements into one time series (merge; default) or
        combine: keeps hand and diver measurements separeted       
        The default is merge.
    url_lizard : str
        location of the LIZARD-API.

    Returns
    -------
    obs_col : TYPE
        ObsCollection DataFrame with the 'obs' column

    '''
    
    if type(extent) == list:
        polygon_T = extent_to_polygon(extent)
    

    
    elif extent.endswith('.shp'):
        polygon = geopandas.read_file(extent)
        polygon_T = polygon.to_crs("WGS84","EPSG:28992").loc[0,"geometry"]
        
        
    else:
        print("Extent should be a shapefile or a list of coordinates")
        return
   
    r_results, groundwaterstation_info= get_API_results_from_extent(polygon_T, url_lizard)
    
    obs_col = pd.DataFrame(columns = ['monitoring_well', 'tube_nr', 'name', 'x', 'y', 'tube_top',
                                      'ground_level', 'screen_top', 'screen_bottom', 'status',
                                      'timeseries_available', 'start_hand', 'start_diver', 'source', 
                                      'unit','obs'])
    
    groundwaterstation_info = [(series, 'all', url_lizard) for _, series in groundwaterstation_info.iterrows()]
    groundwaterstation_filters = pd.DataFrame()
        
    nr_threads = 10            
    if nr_threads > r_results:
        nr_threads = r_results
    
    with concurrent.futures.ThreadPoolExecutor(max_workers= nr_threads) as executor:
        for result in tqdm(executor.map(lambda args : get_metadata_filter(*args), groundwaterstation_info),
                           total = r_results, desc='Monitoring well'):
            groundwaterstation_filters = pd.concat([groundwaterstation_filters,pd.DataFrame(result)])
    
    if extract_timeseries == True:
        groundwaterstation_filters = [
            (series, tmin,tmax,type_timeseries, url_lizard) for _, series in groundwaterstation_filters.iterrows()]
              
        with concurrent.futures.ThreadPoolExecutor(max_workers= nr_threads) as executor:
            for measurement, obs_series in tqdm(executor.map(lambda args : extract_timeseries_from_API(*args), groundwaterstation_filters),
                               total =len(groundwaterstation_filters), desc='Timeseries'):
                obs_series['obs'] = measurement                
                obs_col = pd.concat([obs_col,pd.DataFrame([obs_series])])
                
    else:
        groundwaterstation_filters.drop(['uuid_hand','uuid_diver'], axis = 1, inplace = True)
        
        obs_col = groundwaterstation_filters
    
    obs_col.reset_index(drop = True, inplace = True)
    
    return obs_col
        
      

        
        
        
        
        
        
