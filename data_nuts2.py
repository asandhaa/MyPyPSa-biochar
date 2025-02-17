import os
import pandas as pd
import numpy as np
import math
import logging
logger = logging.getLogger(__name__)
from numpy import pi as pi
from numpy import sin as sin
from numpy import cos as cos 
from numpy import arccos as arccos
import json
import time
import requests

print('enter data network name')
network_name=str(input())


name=network_name[:-3] 

def add_missing_gens(n, correct): 
    """
    Funktion die überprüft, ob es buses ohne conventionelle pp gibt, jedoch laut power plant data, welche geben müsste
    Halte die Funktion so allgemein, dass sie auch für andere konventionelle verwendet werden kann"""
    missing_gens = [item for item in correct.index if item not in n.generators.index]
   
    for gens in missing_gens:
        carrier = correct.loc[gens].carrier
        n.add("Generator",gens,bus=gens[:6]) 
        n.generators.loc[gens] = n.generators[n.generators.carrier == carrier].iloc[0] 
        n.generators.loc[gens, 'p_nom'] = 0
        n.generators.loc[gens, 'bus'] = gens[:6] 
        if carrier == 'ror':
            n.generators_t.p_max_pu[gens] = list(n.generators_t.p_max_pu.loc[:,n.generators_t.p_max_pu.columns.str.contains('ror')].iloc[:,0])
    return missing_gens
        
#Step 1 for power plant clustering
def cluster_by_postcode(data):
    """ Location of Power plants are assigned to NUTS1-cluster according their postcode
    This works only for NUTS1-clustering, for a NUTS3-clustering the CODE and bus names are already the same
    
    INPUT: data 
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'postcode']
    OUTPUT: data_1, data_nan1
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    """
    # a) all power plants without postcode are placed in data_nan1 and removed from data
    data_nan1= data[data['postcode'].isna()==True]
    data.dropna(subset=['postcode'], inplace =True)
    
    data.postcode = data.postcode.astype(int)
    plz_nuts3 = pd.read_csv(r'1_Regionality/plz_nuts3.csv', header = 0, sep=';')
    plz_nuts3['CODE'] =plz_nuts3['CODE'].str.replace("'", "")
    plz_nuts3['bus']= ''
    
    #Translate CODE (which is NUTS3) into NUTS2 bus names:
    for row in plz_nuts3.index:
        plz_nuts3.loc[row,'bus'] = 'DE0 ' + plz_nuts3.iloc[row]['NUTS3'][3:5] 
    plz_nuts1 = plz_nuts3.drop(['NUTS3'], axis=1)
    plz_nuts1.CODE=plz_nuts1.CODE.astype(int)
    data_1 = data.merge(plz_nuts1, left_on="postcode", right_on="CODE", how = 'left')
    
    # b) all power plants without matching NUTS1 are placed in data_nan1 and removed from merged
    data_nan1 = data_nan1.append(data_1[data_1['bus'].isna()==True])
    data_1.dropna(subset = ['bus'], inplace = True)
    
    data_nan1.drop( columns =['postcode', 'CODE'], inplace = True)
    data_1.drop( columns =['postcode', 'CODE'], inplace = True)
    
    return data_1, data_nan1
    

#Step 2 for power plant clustering
def cluster_by_geoshape(data_nan1):
    """INPUT: data_nan1
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    OUTPUT: data_nan2
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    """
    import geopandas as gpd
    from shapely.geometry import Point

    # Load the NUTS shapefile
    nuts_shapefile_path = r'data\NUTS_2013_60M_SH\data\NUTS_RG_60M_2013.shp' #'NUTS_RG_60M_2013.shp' #
    nuts_gdf = gpd.read_file(nuts_shapefile_path)
    nuts_gdf.drop(nuts_gdf[nuts_gdf['STAT_LEVL_'] != 2].index, inplace=True) 
    nuts_gdf=nuts_gdf[nuts_gdf['NUTS_ID'].str.contains('DE')]
    
    buses = []
    for i in data_nan1.index:
        
        latitude = data_nan1.loc[i].latitude
        longitude = data_nan1.loc[i].longitude
        point = Point(longitude, latitude)
        
        matching_nuts = nuts_gdf[nuts_gdf.geometry.contains(point)]['NUTS_ID']
        
        if not matching_nuts.empty:
            buses.append(matching_nuts.iloc[0])
        else:
            buses.append(np.nan)
             
    data_nan1['bus'] = buses
    
    data_nan2 = data_nan1[data_nan1['bus'].isna()==True]
    data_2 = data_nan1[data_nan1['bus'].isna()==False]

    data_2.bus = data_2['bus'].apply(lambda x: 'DE0 ' + x[2:])
    
    return data_2, data_nan2

#Step 3 for power plant clustering
def cluster_by_distance(n, data_nan2):
    """INPUT: data_nan2
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    OUTPUT: data_3
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    """
    a=n.buses.loc[n.generators.bus[n.generators.carrier=='solar'], ['x','y']] 
    a['bus']=a.index
    a.index=range(len(a))
    a.y*=pi/180
    a.x*=pi/180
    data_nan2.longitude*=pi/180
    data_nan2.latitude*=pi/180
    
    
    R=6371 # mean earth radius in km
    buses=[]
    #Great-Circle Distance Formula
    for j in range(len(data_nan2)):
        distances=[]
        for i in range(len(a)):
            distances.append(arccos(sin(a.y.iloc[i]) * sin(data_nan2.latitude.iloc[j]) + cos(a.y.iloc[i]) *
                                    cos(data_nan2.latitude.iloc[j]) * cos(a.x[i] - data_nan2.longitude.iloc[j])) * R)
        df=pd.DataFrame({'bus':a.bus,'distance':distances})
        buses.append(df.bus[np.where(df.distance==min(df.distance))[0][0]])
    
    data_3 = data_nan2.copy()
    data_3['bus']=buses
    
    return data_3 


#Only necessary for NUTS2 --> 37 clusters
def city_into_nuts3(data):
    data.bus = data.bus.str.replace("DE0 50","DE0 93") # Bremen into Lüneburg, Stat. Region (Niedersachsen)


def Biomass_data (n,name):
    if os.path.exists('{}/bio_basic_removal.csv'.format(name)):
        remove=pd.read_csv('{}/bio_basic_removal.csv'.format(name),index_col=0)
    else:
        data=pd.read_csv('data/renewable_power_plants_DE_bio.csv')
        data=data[['commissioning_date','decommissioning_date','technology','electrical_capacity','lat','lon','postcode']]
        data.columns=['year_added','year_removed','carrier','p_nom','latitude','longitude','postcode']
        data.carrier=data.carrier.replace('Sewage and landfill gas', 'biomass')
        data.carrier=data.carrier.replace('Biomass and biogas', 'biomass')
        data.drop(data[data.carrier != 'biomass' ].index, inplace=True)
        data.index=range(len(data))
    
        year_added=[]
        year_removed=[]
        for i in range(len(data)):
            year_added.append( int(data[data.columns[0]][i][:4]) )
            year_removed.append( int(data[data.columns[0]][i][:4]) +30)
        
        data.year_added=year_added
        data.year_removed=year_removed
        
        data.drop(data[data.year_added < 1989 ].index, inplace=True)
        data.drop(data[data.year_added > 2020].index, inplace=True)
    
        data.drop(data[np.isnan(data.longitude) ==True].index, inplace=True)
        data.drop(data[np.isnan(data.latitude) ==True].index, inplace=True)
        data.drop(data[np.isnan(data.p_nom) ==True].index, inplace=True)
        
        data_backup = data.copy()
        
        data = data_backup.copy()
        
        """Assigning coordinates to NUTS1 """ 
        #Step 1) Assigning according to postcode into dataframe merged_data, where there is no matching is displayed in data_nan
        """INPUT: data 
            columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'postcode']
        OUTPUT: data_1, data_nan1
            columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
        """
        data_1, data_nan1 = cluster_by_postcode(data)
        
        #Step 2) Assigning according to point in or nearby next geoshape
        data_2, data_nan2 = cluster_by_geoshape(data_nan1) 
        data_2 = pd.concat([data_1,data_2], ignore_index = True)
        
        #Step 3) Assigning according Great Circle 
        data_3 = cluster_by_distance(n, data_nan2)
        data = pd.concat([data_2,data_3], ignore_index = True)
        ######End#####
        
        
        city_into_nuts3(data)
        data=data[['carrier' ,'p_nom','bus','year_removed' ]]
        data.index=data.bus + ' ' + data.carrier

        ### Removal     
        removal=data[['year_removed','carrier','p_nom','bus']]
        years_res=list(removal.year_removed.value_counts().index)
        years_res.sort()
        remove= pd.DataFrame(columns = removal.columns)    
        for i in years_res:
            dat=removal.loc[removal.year_removed == i]
        
            df = dat.groupby(['year_removed', 'bus','carrier']).sum().sum(
                level=['year_removed', 'bus','carrier']).unstack('year_removed').fillna(0).reset_index()
        
            df['year']=[df.columns[2][1]]*len(df)
            df.columns=['bus','carrier','p_nom','year_removed']
            remove= pd.concat([df, remove], ignore_index=False)
        
        remove.index= remove.bus + ' ' + remove.carrier
        remove.to_csv('{}/bio_basic_removal.csv'.format(name))
    
    correct=remove.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    correct=correct.stack().reset_index()
    correct=correct[['carrier','p_nom','bus']]
    correct.index=correct.bus + ' ' + correct.carrier
    
    
    #1) Änderung AS: Add gens to network
    missing_gens = add_missing_gens(n, correct) #but not necessary
    
    #2) Änderung AS: Drop gens from network 
    missing_gens = [item for item in n.generators[n.generators.carrier == 'biomass'].index if item not in correct.index ]
    for gen in missing_gens:
        n.generators.drop(gen, inplace=True)    
      
            
    for i in correct.index:
        logger.info("Correcting {} capacity: {} MW instead of {} MW".format(i,
                                                                            round(correct.loc[i,'p_nom'],3),
                                                                            round(n.generators.loc[i,'p_nom'],3)))
        n.generators.loc[i,'p_nom']=correct.loc[i,'p_nom']


    addition=remove[['carrier','p_nom','bus']]
    addition.index=addition.bus + ' ' + addition.carrier
    add=addition.set_index(['carrier','bus']).stack().reset_index()
    add.columns=['carrier','bus','remove','p_nom']
    add=add.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    bus=[]
    carrier=[]
    idx=[]
    for i in add.index:
        carrier.append(i[0])
        bus.append(i[1])
        idx.append(i[1] + ' ' + i[0])
    
    add['carrier']=carrier
    add['bus']=bus
    add.index=idx
    add.columns=['p_nom','carrier','bus']
    add.to_csv('{}/biomass_basic_addition.csv'.format(name))


    return remove



def RES_data(n,name):
    
    if os.path.exists('data/renewable_power_plants_DE.csv'):
        data=pd.read_csv('data/renewable_power_plants_DE.csv')
    else:
        URL= 'https://zenodo.org/record/6826283/files/renewable_power_plants_DE.csv?download=1'
        data=pd.read_csv(URL)
        data.to_csv('data/renewable_power_plants_DE.csv')
    data=data[['commissioning_date','decommissioning_date','technology','electrical_capacity','lat','lon', 'postcode']] #Änderung AS
    data.columns=['year_added','year_removed','carrier','p_nom','latitude','longitude', 'postcode']
    data.carrier=data.carrier.replace('Photovoltaics', 'solar')
    data.carrier=data.carrier.replace('Onshore', 'onwind')
    data.carrier=data.carrier.replace('Photovoltaics ground', 'solar')
    data.carrier=data.carrier.replace('Run-of-river', 'ror')
    data.carrier=data.carrier.replace('Sewage and landfill gas', 'biomass')
    data.carrier=data.carrier.replace('Biomass and biogas', 'biomass')
    data.drop(data[data.carrier == 'biomass' ].index, inplace=True)

    """
    Change Offshore data
    """
    a=data[data.carrier=='Offshore']
    vals=pd.DataFrame({'lat_vals':a.latitude.unique()})
    vals.drop(vals[np.isnan(vals.lat_vals) ==True].index, inplace=True)
    a.latitude[np.isnan(a.latitude)]=vals.mean()[0]
    
    vals=pd.DataFrame({'long_vals':a.longitude.unique()})
    vals.drop(vals[np.isnan(vals.long_vals) ==True].index, inplace=True)
    
    a.longitude[np.isnan(a.longitude)] = np.random.choice(vals.long_vals,
                size=len(a.longitude[np.isnan(a.longitude)]))
    j=0
    for i in a.index:
        if a.carrier[i] == 'Offshore':
            j=j+1
            if j % 2 ==0:
                a.carrier[i]=('offwind-ac')
            else:
                a.carrier[i]=('offwind-dc')
    
    for i in a.index:
        data.loc[i,'carrier']=a.loc[i,'carrier']
        data.loc[i,'latitude']=a.loc[i,'latitude']
        data.loc[i,'longitude']=a.loc[i,'longitude']

    """ """

    data.index=range(len(data))
    
    year_added=[]
    year_removed=[]
    for i in range(len(data)):
        year_added.append( int(data[data.columns[0]][i][:4]) )
        year_removed.append( int(data[data.columns[0]][i][:4]) +25)
    
    data.year_added=year_added
    data.year_removed=year_removed
    
    data.drop(data[data.year_added < 1995 ].index, inplace=True)
    data.drop(data[data.year_added > 2018].index, inplace=True)
    data.drop(data[data.carrier == 'Other fossil fuels' ].index, inplace=True)
    data.drop(data[data.carrier == 'Geothermal' ].index, inplace=True)
    data.drop(data[data.carrier == 'Storage' ].index, inplace=True)
    

    data.drop(data[np.isnan(data.longitude) ==True].index, inplace=True)
    data.drop(data[np.isnan(data.latitude) ==True].index, inplace=True)
    data.drop(data[np.isnan(data.p_nom) ==True].index, inplace=True)
    
    data.drop(data[data.carrier=='ror'].index, inplace=True)
    data.drop(data[np.isnan(data.longitude) ==True].index, inplace=True)
    data.drop(data[np.isnan(data.latitude) ==True].index, inplace=True)
    data.index=range(len(data))

    ### 2019 data
    df2=pd.read_csv('data/renewable_power_plants_DE_2019.csv')
    df2=df2[['commissioning_date','decommissioning_date','technology','electrical_capacity','lat','lon','postcode']] #Änderung AS
    df2.columns=['year_added','year_removed','carrier','p_nom','latitude','longitude','postcode']
    df2.carrier=df2.carrier.replace('Photovoltaics', 'solar')
    df2.carrier=df2.carrier.replace('Onshore', 'onwind')
    df2.carrier=df2.carrier.replace('Photovoltaics ground', 'solar')
    df2.carrier=df2.carrier.replace('Run-of-river', 'ror')
    df2.carrier=df2.carrier.replace('Sewage and landfill gas', 'biomass')
    df2.carrier=df2.carrier.replace('Biomass and biogas', 'biomass')
    df2.drop(df2[df2.carrier == 'biomass' ].index, inplace=True)

    # Correct offiwnd
    a_19=df2[df2.carrier=='Offshore']
    vals=pd.DataFrame({'lat_vals':a.latitude.unique()})
    vals.drop(vals[np.isnan(vals.lat_vals) ==True].index, inplace=True)
    a_19.latitude[np.isnan(a_19.latitude)]=vals.mean()[0]
    
    vals=pd.DataFrame({'long_vals':a.longitude.unique()})
    vals.drop(vals[np.isnan(vals.long_vals) ==True].index, inplace=True)
    
    a_19.longitude[np.isnan(a_19.longitude)] = np.random.choice(vals.long_vals,
                size=len(a_19.longitude[np.isnan(a_19.longitude)]))
    j=0
    for i in a_19.index:
        if a_19.carrier[i] == 'Offshore':
            j=j+1
            if j % 2 ==0:
                a_19.carrier[i]=('offwind-ac')
            else:
                a_19.carrier[i]=('offwind-dc')
    
    for i in a_19.index:
        df2.loc[i,'carrier']=a_19.loc[i,'carrier']
        df2.loc[i,'latitude']=a_19.loc[i,'latitude']
        df2.loc[i,'longitude']=a_19.loc[i,'longitude']
        
    
    df2.drop(df2[np.isnan(df2.longitude) ==True].index, inplace=True)
    df2.drop(df2[np.isnan(df2.latitude) ==True].index, inplace=True)
    df2.drop(df2[np.isnan(df2.p_nom) ==True].index, inplace=True)
    df2.drop(df2[df2.carrier == 'ror' ].index, inplace=True)
    df2.drop(df2[df2.carrier == 'Geothermal' ].index, inplace=True)
    
    
    df2.index=range(len(df2))
    
    year_added=[]
    year_removed=[]
    for i in range(len(df2)):
        year_added.append( int(df2[df2.columns[0]][i][:4]) )
        year_removed.append( int(df2[df2.columns[0]][i][:4]) +25)
    
    df2.year_added=year_added
    df2.year_removed=year_removed    
    
    
    data=pd.concat([data,df2])
    data.index=range(len(data))

    
    data_backup = data.copy()
    
    data = data_backup.copy()
    
    """Assigning coordinates to NUTS1 """ 
    #Step 1) Assigning according to postcode into dataframe merged_data, where there is no matching is displayed in data_nan
    """INPUT: data 
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'postcode']
    OUTPUT: data_1, data_nan1
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    """
    data_1, data_nan1 = cluster_by_postcode(data)
    
    #Step 2) Assigning according to point in or nearby next geoshape
    data_2, data_nan2 = cluster_by_geoshape(data_nan1) 
    data_2 = pd.concat([data_1,data_2], ignore_index = True)
    
    #Step 3) Assigning according Great Circle 
    data_3 = cluster_by_distance(n, data_nan2)
    data = pd.concat([data_2,data_3], ignore_index = True) 
    
    city_into_nuts3(data)   
    data=data[['carrier' ,'p_nom','bus','year_removed' ]]
    data.index=data.bus + ' ' + data.carrier
    
   
    addition=data[['carrier','p_nom','bus']]
    add=addition.set_index(['carrier','bus']).stack().reset_index()
    add.columns=['carrier','bus','remove','p_nom']
    add=add.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    bus=[]
    carrier=[]
    idx=[]
    for i in add.index:
        carrier.append(i[0])
        bus.append(i[1])
        idx.append(i[1] + ' ' + i[0])
    
    add['carrier']=carrier
    add['bus']=bus
    add.index=idx
    add.columns=['p_nom','carrier','bus']
    add.to_csv('{}/res_basic_addition.csv'.format(name))
    
    removal=data[['year_removed','carrier','p_nom','bus']]
    years_res=list(removal.year_removed.value_counts().index)
    years_res.sort()
    remove= pd.DataFrame(columns = removal.columns)    
    for i in years_res:
        dat=removal.loc[removal.year_removed == i]
    
        df = dat.groupby(['year_removed', 'bus','carrier']).sum().sum(
            level=['year_removed', 'bus','carrier']).unstack('year_removed').fillna(0).reset_index()
    
        df['year']=[df.columns[2][1]]*len(df)
        df.columns=['bus','carrier','p_nom','year_removed']
        remove= pd.concat([df, remove], ignore_index=False)
    
    remove.index= remove.bus + ' ' + remove.carrier
    remove.to_csv('{}/res_basic_removal.csv'.format(name))

    return add,remove

def Correct_coal(n,name):
    if os.path.exists('{}/coal_basic_removal.csv'.format(name)):   
        data= pd.read_csv('{}/coal_basic_removal.csv'.format(name),index_col=0) 
    else:
        data = pd.read_csv('data/conventional_power_plants_DE.csv')
        data =data [['capacity_net_bnetza','energy_source','commissioned','retrofit','status','lat','lon','postcode']]
        data=data.rename(columns={"capacity_net_bnetza": "p_nom", "energy_source": "carrier"})
        data.drop(data[data.carrier !='Hard coal'].index, inplace=True)
        data.drop(data[data.status =='shutdown'].index, inplace=True)
        data.carrier=data.carrier.replace('Hard coal', 'coal')
        data.index=range(len(data))
        
        for i in range(len(data)):
            if data.retrofit[i] > data.commissioned[i]:
                print(i)
                data.commissioned[i]=data.commissioned[i] + 20
    
        data.drop(data[data.commissioned<1969].index, inplace=True)
        data['year_removed']=data.commissioned + 50
    
        data.index=range(len(data))
        data=data.rename(columns={"commissioned": "year_added", "lat": "latitude", "lon": "longitude"})
        data.at[31,'postcode'] =data.at[31,'postcode'][-5:] 
        data.postcode = data.postcode.astype(int)
          
            

        data_backup = data.copy()
        
        data = data_backup.copy()
        
        """Assigning coordinates to NUTS1 """ 
        #Step 1) Assigning according to postcode into dataframe merged_data, where there is no matching is displayed in data_nan
        """INPUT: data 
            columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'postcode']
        OUTPUT: data_1, data_nan1
            columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
        """
        data_1, data_nan1 = cluster_by_postcode(data)
        
        #Step 2) Assigning according to point in or nearby next geoshape
        data_2, data_nan2 = cluster_by_geoshape(data_nan1) 
        data_2 = pd.concat([data_1,data_2], ignore_index = True)
        
        #Step 3) Assigning according Great Circle
        data_3 = cluster_by_distance(n, data_nan2)
        data = pd.concat([data_2,data_3], ignore_index = True)

        
        city_into_nuts3(data)
        data=data[['carrier' ,'p_nom','bus','year_removed' ]]
        data.index=data.bus + ' ' + data.carrier
        data=data[data.year_removed >=2020]
        data.to_csv('{}/coal_basic_removal.csv'.format(name)) 
        
    correct=data.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    correct=correct.stack().reset_index()
    correct=correct[['carrier','p_nom','bus']]
    correct.index=correct.bus + ' ' + correct.carrier
    
    #1) Änderung AS: Add gens to network
    missing_gens = add_missing_gens(n, correct)
    
    #2) Änderung AS: Set for generators not in correct p_nom to 0 
    missing_gens = [item for item in n.generators[n.generators.carrier == 'coal'].index if item not in correct.index ]
    for gen in missing_gens:
        n.generators.loc[gen, 'p_nom'] = 0
    
    for i in correct.index:
        logger.info("Correcting {} capacity: {} MW instead of {} MW".format(i,
                                                                            round(correct.loc[i,'p_nom'],3),
                                                                            round(n.generators.loc[i,'p_nom'],3)))
        n.generators.loc[i,'p_nom']=correct.loc[i,'p_nom']


    

    return data

def Base_Removal_Data(n):
    data = pd.read_csv('data/ppl.csv')
    df1 = data[data.carrier == 'hydro']
    df2 = data[data.carrier != 'hydro']
    df1.drop(df1[df1.technology !='Run-Of-River'].index, inplace=True)
    data= pd.concat([df1,df2])
    
    data=data[['carrier','p_nom','lat','lon','yeardecommissioning','yearcommissioned','retrofit']] 
    data=data.rename(columns={"yearcommissioned": "year_added","yeardecommissioning": "year_removed", "lat": "latitude", "lon": "longitude"})

    
    data.drop(data[data.carrier =='nuclear'].index, inplace=True)
    data.drop(data[data.carrier =='biomass'].index, inplace=True)
    data.drop(data[data.carrier =='other'].index, inplace=True)
    data.drop(data[data.carrier =='waste'].index, inplace=True)
    data.drop(data[data.carrier =='storage technologies'].index, inplace=True)
    data.carrier=data.carrier.replace('hydro','ror')
    df=pd.DataFrame({'p_nom':data.p_nom,'carrier':data.carrier})

    data.index=range(len(data))

    data_backup = data.copy()
    
    data = data_backup.copy()
    
    """Assigning coordinates to NUTS1 """ 
    #Step 1) Assigning according to postcode into dataframe merged_data, where there is no matching is displayed in data_nan
    """INPUT: data 
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'postcode']
    OUTPUT: data_1, data_nan1
        columns=['year_added', 'year_removed', 'carrier', 'p_nom', 'latitude','longitude', 'bus']
    """
    # data_1, data_nan1 = cluster_by_postcode(data)
    data_nan1 = data.copy()
    #Step 2) Assigning according to point in or nearby next geoshape
    data_2, data_nan2 = cluster_by_geoshape(data_nan1) 
    # data_2 = pd.concat([data_1,data_2], ignore_index = True)
    
    #Step 3) Assigning according Great Circle  
    data_3 = cluster_by_distance(n, data_nan2)
    data = pd.concat([data_2,data_3], ignore_index = True)
    
    data = data[['carrier', 'p_nom', 'year_removed','retrofit', 'bus']]
    data.drop(data[data.carrier =='coal'].index, inplace=True)

    data.drop(data[data.year_removed <=2020].index, inplace=True)

    ### Correcting CCGT & OCGT
    pre_CCGT=round(data.p_nom[data.carrier=='CCGT'].sum(),3)
    pre_OCGT=round(data.p_nom[data.carrier=='OCGT'].sum(),3)

    data.p_nom[data.carrier=='CCGT']=data.p_nom[data.carrier=='CCGT']*1.25921273 
    data.p_nom[data.carrier=='OCGT']=data.p_nom[data.carrier=='OCGT']*1.25921273

    logger.info("Correcting CCGT capacity: {} MW instead of {} MW".format(round(data.p_nom[data.carrier=='CCGT'].sum(),3),pre_CCGT))
    logger.info("Correcting OCGT capacity: {} MW instead of {} MW".format(round(data.p_nom[data.carrier=='OCGT'].sum(),3),pre_OCGT))

    city_into_nuts3(data)
    data.index=data.bus + ' ' + data.carrier

    #1) Änderung AS: Add gens to network
    missing_gens = add_missing_gens(n, data) 

    #2) Änderung AS: Set p_nom for gens not in data to zero
    missing_gens = [item for item in n.generators[n.generators.carrier.isin(['lignite','ror','oil'])].index if item not in data.index]
    for gen in missing_gens:
        n.generators.loc[gen, 'p_nom'] =0


        
    #Correct 'lignite','ror','oil'    
    correct=data.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    correct=correct.stack().reset_index()
    correct=correct[['carrier','p_nom','bus']]
    correct=correct[correct.carrier.isin(['lignite','ror','oil'])]
    correct.index=correct.bus + ' ' + correct.carrier 
    
    for i in correct.index:
        logger.info("Correcting {} capacity: {} MW instead of {} MW".format(i,
                                                                            round(correct.loc[i,'p_nom'],3),
                                                                            round(n.generators.loc[i,'p_nom'],3)))
        n.generators.loc[i,'p_nom']=correct.loc[i,'p_nom']
    
    data.to_csv('{}/conventional_basic_removal.csv'.format(name))
    return data


def resample_gen_profiles(n,net_name,resample_factor):
    gen_profiles=pd.read_csv('{}/gen_profiles.csv'.format(net_name),index_col=0,parse_dates=True)
    gen_profiles=gen_profiles.resample('{}H'.format(resample_factor)).mean()

    for i in n.generators.index[n.generators.carrier=='solar']:
        if (i not in gen_profiles.columns) | (math.isnan(gen_profiles[i][0])):
            continue
        else:
            data= pd.DataFrame ({'electricity':gen_profiles[i]})
            logger.info("Correcting {} generation profile: {} CF instead of {}".format(i,
                                                                        round(data.electricity.mean(),3),
                                                                        round(n.generators_t.p_max_pu[i].mean(),3)))

    for i in n.generators.index[n.generators.carrier=='onwind']:
        if (i not in gen_profiles.columns) | (math.isnan(gen_profiles[i][0])):
            continue
        else:
            data= pd.DataFrame ({'electricity':gen_profiles[i]})
            logger.info("Correcting {} generation profile: {} CF instead of {}".format(i,
                                                                        round(data.electricity.mean(),3),
                                                                        round(n.generators_t.p_max_pu[i].mean(),3)))

def update_availability_profiles(n,name):              
    gen_profiles=pd.read_csv('{}/gen_profiles.csv'.format(name),index_col=0)
    for gen in gen_profiles.columns:
        logger.info("Correcting {} generation profile: {} CF instead of {}"
                    .format(gen,round(gen_profiles[gen].mean(),3),
                            round(n.generators_t.p_max_pu[gen].mean(),3)))
        n.generators_t.p_max_pu[gen]=gen_profiles[gen].values

def update_rens_profiles(n,reference_year,name):       
    
    token = '207fc85e63dc9da68510ae062a90511d239b7497' #'1cacc0c037ec8d41e4ed821bcd9d1ab673766aee'#'daa9bc54f7bc4af9ffd1124602f6dd8a206b60ca '
    api_base = 'https://www.renewables.ninja/api/'
    
    s = requests.session()
    s.headers = {'Authorization': 'Token ' + token}
    
    if os.path.exists('{}/gen_profiles.csv'.format(name)):
        gen_profiles=pd.read_csv('{}/gen_profiles.csv'.format(name),index_col=0)
        gen_profiles.index =n.generators_t.p_max_pu.index
    else:
        gen_profiles= pd.DataFrame(index =n.generators_t.p_max_pu.index,
                                   columns=list(n.generators.index[(n.generators.carrier=='solar') |
                                                                   (n.generators.carrier=='onwind')]))
    args_wind = {
        'date_from': str(reference_year)+'-01-01',
        'date_to': str(reference_year)+'-12-31',
        'capacity': 1.0,
        'height': 101.5, # min:84, max:119, avg:101.5
        'turbine': 'Vestas V112 3000',
        'format': 'json'}
    args_solar = {
        'date_from': str(reference_year)+'-01-01',
        'date_to': str(reference_year)+'-12-31',
        'dataset': 'merra2',
        'capacity': 1.0,
        'system_loss': 0.1,
        'tracking': 0,
        'tilt': 35,
        'azim': 180,
        'format': 'json'
        }
    args = {'solar':args_solar, 'onwind':args_wind}
    urls= {'solar':'pv', 'onwind':'wind'}
    for tech in ['solar','onwind']:
        for i in n.generators.index[n.generators.carrier==tech]:
            if (i not in gen_profiles.columns) | (math.isnan(gen_profiles[i][0])):
                args[tech]['lat']=n.buses.loc[n.generators.loc[i,'bus'],'y']
                args[tech]['lon']=n.buses.loc[n.generators.loc[i,'bus'],'x']
                url = api_base + 'data/{}'.format(urls[tech])
                r = s.get(url, params=args[tech])
                parsed_response = json.loads(r.text)
                data = pd.read_json(json.dumps(parsed_response['data']), orient='index')
                resample_factor=int(8760/len(n.snapshots))
                gen_profiles[i]=data.electricity.resample('{}H'.format(resample_factor)).mean()
                gen_profiles.to_csv('{}/gen_profiles.csv'.format(name))
                time.sleep(15)



