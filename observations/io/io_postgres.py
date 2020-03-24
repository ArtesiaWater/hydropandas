import pandas as pd
from sqlalchemy import (Boolean, Column, DateTime, Float, Integer, MetaData,
                        String, Table, create_engine)
from sqlalchemy.ext.automap import automap_base
# from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import observations as obs

engine = create_engine(
    'postgresql://postgres:Io0bgSd9@127.0.0.1:5432/test_men2', echo=True, paramstyle="format")


ts_table_name = '"StandpipeMeas"'
pipe_id = 581
sql = """ SELECT * FROM {} WHERE "PipeId" = {};""".format(
    ts_table_name, pipe_id)
#sql = """ SELECT * FROM {};""".format('"ObservationWells"')


raw_series = pd.read_sql(
    sql, con=engine, index_col='DateTime', parse_dates=[0])
# %%
Base = automap_base()
Base.prepare(engine, reflect=True)

DBSession = sessionmaker(bind=engine)
session = DBSession()


obs_well_id = session.query(
    Base.classes.ObservationWells.ObsWellId).all()[0][0]

obs_dic = session.query(Base.classes.ObservationWells).get(
    obs_well_id).__dict__
loc_dic = session.query(Base.classes.Locations).get(
    obs_dic['LocationId']).__dict__
sp_dic = session.query(Base.classes.Standpipes).filter(
    Base.classes.Standpipes.ObsWellId == obs_well_id)

obs_dic.pop('_sa_instance_state')

o = obs.GroundwaterObs(x=loc_dic.pop('XCoordinate'), y=loc_dic.pop('YCoordinate'),
                       maaiveld=loc_dic.pop('SurfaceLevel'), )
