# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 22:18:14 2023

@author: asandhaa
"""




"""Plotting"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import yaml
import json
import pypsa
import matplotlib as mpl
import cartopy.crs as ccrs
from matplotlib.offsetbox import AnchoredText
import re
import data
import glob
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper, FullLoader

import pyrolysis 
import geopandas as gpd


with open(r'config.yaml') as file:
    config= yaml.load(file, Loader=yaml.FullLoader)
    
network_name=data.network_name

name=network_name[:-3]

"""Network"""

pp_list = ['battery','H2','hydro', 'PHS',
           'lignite', 'coal', 'oil', 'OCGT', 'CCGT', 'geothermal',
            'ror', 'biomass', 'offwind-dc', 'offwind-ac','onwind', 'solar',  'pyrolysis']  


new_names = ['batteries', 'hydrogen', 'reservoir & dam','pumped hydro storage',
             'lignite', 'hard coal', 'oil', 'open-cycle gas', 'combined-cycle gas', 'geothermal', 
             'run of river', 'biomass', 'wind offshore (dc)', 'wind offshore (ac)', 'wind onshore', 'solar PV', 'pyrolysis']

new_names_big = ['Batteries', 'Hydrogen', 'Reservoir & Dam','Pumped Hydro Storage',
             'Lignite', 'Hard Coal', 'Oil', 'Open-Cycle Gas', 'Combined-Cycle Gas', 'Geothermal', 
             'Run of River', 'Biomass', 'Wind Offshore (DC)', 'Wind Offshore (AC)', 'Wind Onshore', 'Solar PV', 'Pyrolysis']

neue_namen = ['Batterien', 'Wasserstoff', 'Reservoir & Dämme','Pumpspeicher',
             'Braunkohle', 'Steinkohle', 'Mineralöl', 'Gasturbinen', 'GuD Kraftwerk', 'Geothermal', 
             'Laufwasser', 'Biomasse', 'Wind Offshore (DC)', 'Wind Offshore (AC)', 'Wind Onshore', 'Solar PV', 'Pyrolyse']

neue_namen_gross = ['Batterien', 'Wasserstoff', 'Reservoir & Dämme','Pumpspeicher',
             'Braunkohle', 'Steinkohle', 'Mineralöl', 'Gasturbinen', 'GuD Kraftwerk', 'Geothermal', 
             'Laufwasser', 'Biomasse', 'Wind Offshore (DC)', 'Wind Offshore (AC)', 'Wind Onshore', 'Solar PV', 'Pyrolyse']

translate_carriers = pd.DataFrame(new_names, pp_list, columns =['new_names'])
translate_carriers['new_names_big'] = new_names_big
uebersaetze_carriers = pd.DataFrame(neue_namen, pp_list, columns =['neue_namen'])

def new_labels(labels):
    new_labels = []
    for label in labels:
        new_label = translate_carriers.loc[label, 'new_names'] if label in translate_carriers.index else label
        new_labels.append(new_label)
    return new_labels

def new_labels_big(labels):
    new_labels = []
    for label in labels:
        new_label = translate_carriers.loc[label, 'new_names_big'] if label in translate_carriers.index else label
        new_labels.append(new_label)
    return new_labels

#deutsche Namen
def neue_labels(labels):
    new_labels = []
    for label in labels:
        new_label = uebersaetze_carriers.loc[label, 'neue_namen'] if label in translate_carriers.index else label
        new_labels.append(new_label)
    return new_labels

def colors_map(data, **kwargs):

    tech_colors=pd.DataFrame.from_dict(config['plotting']['tech_colors'], orient='index')
    df = pd.DataFrame({'carrier': tech_colors.index, 'color': tech_colors.iloc[:,0]})
                       
    c=[]
    for i in range(len(data)):
        c.append(list(df.color[df.carrier==data.index[i]])[0])
    return c

def pie_exp(data):
    df=pd.DataFrame({'carrier':['CCGT','OCGT','biomass','coal','lignite','nuclear','offwind-ac',
                                    'offwind-dc','oil','onwind','ror','solar','pyrolysis'], 
                         'exp':[0,0,0,0,0,0,0.1,0.1,0,0.1,0.1,0.1,0.2]}) 
    d=[]
    for i in range(len(data)):
        d.append(list(df.exp[df.carrier==data.index[i]])[0])
    return d

def pie_chart(n,year, PYROLYSIS): #Generation share for one year
    p_by_carrier = pd.DataFrame(index=n.generators_t.p[n.generators.index[n.generators.carrier!='load']].columns, columns=['p_by_carrier'])
    p_by_carrier.p_by_carrier= n.generators_t.p.sum()
    
    if PYROLYSIS==1: p_by_carrier=pyrolysis.pie_chart_pyro(n, p_by_carrier) 
    
    p_by_carrier.p_by_carrier*=int(8760/len(n.snapshots))

    val=[]
    typ=[]
    for i in range(len(p_by_carrier)):
        val.append(p_by_carrier.p_by_carrier[i])
        typ.append(re.split('\s+', p_by_carrier.index[i])[len(re.split('\s+', p_by_carrier.index[i]))-1])
    
    nnam=list(dict.fromkeys(typ))
    nnam = [i for i in pp_list if i in nnam] 
    
    new_gen= []
    for i in range(len(nnam)):  
        x=0
        for j in range(len(val)):
            if typ[j]==nnam[i]:
                x=x+val[j]
        new_gen.append(x)

    new_gen=pd.DataFrame(new_gen, columns=['gens'],index=nnam)
    
    colors = colors_map(new_gen)
    explode = pie_exp(new_gen)
    perc=(new_gen.loc['solar']+new_gen.loc['onwind']+new_gen.loc['offwind-ac']+
          new_gen.loc['offwind-dc']+new_gen.loc['ror'])/new_gen.sum()
    
    fig=plt.figure(figsize=(40, 10))
    plt.title('Generation shares year {} \n Renewables={}%'.format(year,round(perc*100,3)[0]),fontsize=20)
    patches, texts, junk = plt.pie(list(new_gen.gens),explode=explode, colors=colors,
                             autopct='%1.1f%%', counterclock=False, shadow=True)
    labels = new_labels(new_gen.index)
    plt.legend(patches, labels, loc="best")
    plt.savefig('{}/All Generation year {}'.format(name,year), pi=1600, bbox_inches='tight')

    return round(perc*100,3)[0]


def installed_capacities(n,year, PYROLYSIS,elec_eff_pyro):
    df=pd.DataFrame({'p_nom':n.generators.p_nom[n.generators.carrier!='load'],'carrier':n.generators.carrier[n.generators.carrier!='load']})
    
    if PYROLYSIS==1: df=pyrolysis.installed_capacities_pyro(n, df, elec_eff_pyro)   
    
    summation=df.groupby('carrier').sum()
    summation = summation.reindex([i for i in pp_list if i in summation.index]) 
    shares=int(summation.loc['solar']+summation.loc['onwind']+
               summation.loc['offwind-ac']+summation.loc['offwind-dc']+summation.loc['ror'])/sum(summation.p_nom)
    
    colors = colors_map(summation)
    explode = pie_exp(summation)
    fig=plt.figure(figsize=(40, 10))
    plt.title('Installed Capacity year {}\n Total= {} GW, Renewables={} %'.format(year,
              round(sum(summation.p_nom)/1000,3), round(shares*100,3)),fontsize=20)
    patches, texts,junk = plt.pie(list(summation.p_nom),explode=explode, colors=colors,
                             autopct='%1.1f%%', shadow=True, counterclock=False, textprops={'color':"black"})
    labels = new_labels(summation.index)
    plt.legend(patches, labels, loc="best")
    plt.savefig('{}/Installation Shares year {}'.format(name,year), pi=1600, bbox_inches='tight')
    plt.show()
    
    return round(shares*100,3)

def Bar_to_PNG(n,bar,name,bar_type): #Bar_to_PNG(n,inst_bar,name,'Installation')
    
    cols = list(bar.columns.values) 
    cols.pop(cols.index('solar')) 
    cols.pop(cols.index('offwind-dc'))
    cols.pop(cols.index('offwind-ac'))
    cols.pop(cols.index('onwind'))
    cols.pop(cols.index('ror'))
    cols.pop(cols.index('oil'))
    cols.pop(cols.index('coal'))
    cols.pop(cols.index('lignite'))
    
    

    if bar_type =='Generation':
        cols.pop(cols.index('load')) 
        unit='TWh'
        bar*=int(8760/len(n.snapshots))
        bar= bar[ ['lignite','coal', 'oil'] + cols + ['ror','offwind-ac','offwind-dc','onwind','solar','load']]  
        top = 900
        colors = colors_map(pd.DataFrame(index=bar.columns))

    if bar_type =='Installation':
        unit='GW'
        bar= bar [ ['lignite','coal', 'oil'] + cols + ['ror','offwind-ac','offwind-dc','onwind','solar']]
        colors = colors_map(pd.DataFrame(index=bar.columns))
        top = 600
    
    ax = bar.plot(figsize=(40, 10),kind='bar', stacked=True,color=colors,fontsize=16)
    ax.set_xlabel("Years",fontsize=20)
    ax.set_ylabel("{} in {}".format(bar_type,unit),fontsize=20)
    ax.set_ylim(top=top)
    ax.grid()
    labels = new_labels(bar.columns)
    ax.legend(reversed(plt.legend().legendHandles), reversed(labels), loc="upper left", fontsize = 14) 
    plt.show()
    plt.savefig('{}/{} Bar Plot'.format(name,bar_type), bbox_inches='tight') 
    bar.to_excel('{}/{}_Bar.xlsx'.format(name,bar_type))

    

def Country_Map(n,year, PYROLYSIS,elec_eff_pyro):


    shapefile_path= "C:\\Users\\asandhaa\\Desktop\\pypsa-eur\\MyPyPSA-biochar\\1_Regionality\\nuts3_shapes.geojson"
    gdf= gpd.read_file(shapefile_path)
    gdf.drop(gdf[gdf.CNTR_CODE != 'DE'].index,
             inplace=True) 
    
    nuts1 = gdf.drop(gdf[gdf.LEVL_CODE != 2].index)

    nuts1= nuts1.drop(nuts1[nuts1.NUTS_ID == 'DE3'].index)
    nuts1= nuts1.drop(nuts1[nuts1.NUTS_ID == 'DE5'].index)
    nuts1= nuts1.drop(nuts1[nuts1.NUTS_ID == 'DE6'].index)
        
    opts = config['plotting']
    map_figsize = [10,10]
    map_boundaries = opts['map']['boundaries']
    to_rgba = mpl.colors.colorConverter.to_rgba

    line_colors = {'cur': "purple",
                   'exp': mpl.colors.rgb2hex(to_rgba("red", 0.7), True)}
    tech_colors = opts['tech_colors']
    
    bus_sizes = (n.generators.query('carrier != "load"').groupby(['bus', 'carrier']).p_nom.sum())
    
    if PYROLYSIS==1: bus_sizes = pyrolysis.Country_Map_pyro(n,bus_sizes,elec_eff_pyro)  
    
    #Adding storage_units to Country Map 
    bussizes = n.storage_units.query('carrier == "battery"').groupby(['bus', 'carrier']).p_nom.sum()
    bus_sizes = bus_sizes.append(bussizes)
    bussizes = n.storage_units.query('carrier == "H2"').groupby(['bus', 'carrier']).p_nom.sum()
    bus_sizes = bus_sizes.append(bussizes)    
    bussizes = n.storage_units.query('carrier == "PHS"').groupby(['bus', 'carrier']).p_nom.sum()
    bus_sizes = bus_sizes.append(bussizes)    
    bussizes = n.storage_units.query('carrier == "hydro"').groupby(['bus', 'carrier']).p_nom.sum()
    bus_sizes = bus_sizes.append(bussizes)  

    
    line_widths_exp =  n.lines.s_nom_opt

    attribute='p_nom'
    linewidth_factor = opts['map'][attribute]['linewidth_factor']
    bus_size_factor  = opts['map'][attribute]['bus_size_factor']
    
    bus=[]
    car=[]
    for i in bus_sizes.index:
        bus.append(i[0])
        car.append(i[1])
    txt=pd.DataFrame({'Bus':bus,'Carrier':car,'P':bus_sizes})
    txt.index=range(len(bus))
    

    shares=txt.groupby('Carrier').sum()
    res=(shares.P['solar']+shares.P['offwind-ac']+shares.P['offwind-dc']+shares.P['onwind']+shares.P['ror']+shares.P['biomass'])/shares.sum()[0] #Änderung AS: +shares.P['biomass']
    

    map_boundaries=[5,15,46,55]
    n_plo=n.copy()
    
    # Plot the map
    fig, ax = plt.subplots(figsize=map_figsize, subplot_kw={"projection": ccrs.PlateCarree()})
    
    # Plot NUTS1 boundaries with dark grey edges and light grey face color
    nuts1.plot(ax=ax, edgecolor='darkgrey', linewidth=1, facecolor='lightgrey')
    
    # Plot the network data over the NUTS1 boundaries
    n_plo.plot(line_widths=line_widths_exp / linewidth_factor,
               title='Installation Distribution, RES = {}%'.format(round(res * 100, 3)),
               line_colors=pd.Series(line_colors['exp'], n_plo.lines.index),
               bus_sizes=bus_sizes / (2.5 * bus_size_factor),
               bus_colors=tech_colors,
               boundaries=map_boundaries,
               geomap=True,
               color_geomap=True,
               ax=ax)
    ax.add_artist(AnchoredText("{}".format(year), loc=2))
    ax.set_aspect('equal')
    ax.axis('off')
    
    #Änderung AS: include legend
    new_order = [idx for idx in pp_list if idx in shares.index]
    shares = shares.reindex(new_order)
    legend_labels = shares.index
    handles = [mpatches.Patch(label=label, color=colors_map(shares)[i]) for i, label in enumerate(legend_labels)]
    # legend_labels = new_labels(legend_labels) #neue_labels
    ax.legend = ax.legend(handles=reversed(handles), labels=reversed(legend_labels), loc='lower right')#, handlelength=0.8, handletextpad=0.5)
   
       
    # plt.savefig('{}/Installation Map {}'.format(name,year), pi=1600, bbox_inches='tight')
    plt.show()


    return txt

def Generation_BarChart(name):
    networks_names=(glob.glob("{}/*.nc".format(name)))
    n=pypsa.Network(networks_names[0]) ### TODO: try a better way
    carriers=list(dict.fromkeys(n.generators.carrier))
    
    years=[]
    for i in networks_names:
        years.append(i[-7:-3])
    
    bar= pd.DataFrame(index = years,columns=carriers)
    
    for k in range(len(networks_names)):
        print(networks_names[k])
        print(years[k])
    
        n=pypsa.Network(networks_names[k]) 
        p_by_carrier = n.generators_t.p.sum()
        val=[]
        typ=[]
        for i in range(len(p_by_carrier)):
            val.append(p_by_carrier[i])
            typ.append(re.split('\s+', p_by_carrier.index[i])[len(re.split('\s+', p_by_carrier.index[i]))-1])
        nnam=list(dict.fromkeys(typ))
        new_gen= []
        
        for i in range(len(nnam)):
            x=0
            for j in range(len(val)):
                if typ[j]==nnam[i]:
                    x=x+val[j]
            new_gen.append(x)
        
        new_gen=pd.DataFrame(new_gen, index=nnam)
        
        for p in range(len(new_gen.index)):
            bar.loc[years[k],new_gen.index[p]]=new_gen.loc[new_gen.index[p]][0]/10**6
    
    bar=bar.sort_index(axis=0,ascending=True)
    colors = colors_map(pd.DataFrame(index=bar.columns))

    ax = bar.plot(figsize=(40, 10),kind='bar', stacked=True,color=colors,fontsize=15)
    ax.set_xlabel("Years",fontsize=15)
    ax.set_ylabel("Generation [TWh]",fontsize=15)
    ax.grid()
    ax.set_title('Generation shares {}'.format(name),fontsize=30)
    plt.savefig('{}/Generation Bar Plot {}'.format(str(len(n.buses))+'//'+scenario,name), pi=1600, bbox_inches='tight')
    bar.to_csv('{}/Generation_Bar.csv'.format(str(len(n.buses))+'//'+scenario))

    return bar

def Installation_BarChart(name):
    networks_names=(glob.glob("{}/*.nc".format(name)))
    n=pypsa.Network(networks_names[0]) ### TODO: try a better way
    years=[]
    for i in networks_names:
        years.append(i[-7:-3])
    
    df=pd.DataFrame({'p_nom':n.generators.p_nom,'carrier':n.generators.carrier})
    summation=df.groupby('carrier').sum()

    bar= pd.DataFrame(index = years,columns=summation.index)
    
    for k in range(len(networks_names)):
        print(networks_names[k])
        print(years[k])
    
        n=pypsa.Network(networks_names[k]) 
        df=pd.DataFrame({'p_nom':n.generators.p_nom,'carrier':n.generators.carrier})
        summation=df.groupby('carrier').sum()
        for p in range(len(summation.index)):
            bar.loc[years[k],summation.index[p]]=summation.loc[summation.index[p]][0]/10**3

    bar=bar.sort_index(axis=0,ascending=True)
    colors = colors_map(pd.DataFrame(index=bar.columns))
    ax = bar.plot(figsize=(40, 10),kind='bar', stacked=True,color=colors,fontsize=15)
    ax.set_xlabel("Years",fontsize=15)
    ax.set_ylabel("Installation [GW]",fontsize=15)
    ax.grid()
    ax.set_title('Installation shares {}'.format(name),fontsize=30)
    plt.savefig('{}/Installation Bar Plot {}'.format(name,name), pi=1600, bbox_inches='tight')
    
    bar.to_csv('{}/Installation_Bar.csv'.format(str(len(n.buses))+'//'+scenario))

    return bar
#%% Änderungen AS

def p_per_carrier(n, PYROLYSIS):
    columns = [i for i in pp_list if i in list(n.generators.carrier)] #Lists carrier which are used in the correct order
    p_per_carrier = pd.DataFrame(index = n.generators_t.p.index, columns = columns)
    for carrier in p_per_carrier.columns:
        p_per_carrier_type = n.generators_t.p[n.generators_t.p.columns[n.generators_t.p.columns.str.contains(carrier)]].sum(axis=1)
        p_per_carrier[carrier] = p_per_carrier_type
    
    if PYROLYSIS==1: p_per_carrier = pyrolysis.p_per_carrier_pyrolysis(n, p_per_carrier)
    return p_per_carrier

def unit_comitment_plot(n, i, day, minmax, PYROLYSIS):
    # Settings
    data = p_per_carrier(n,PYROLYSIS)
    colors= colors_map(pd.DataFrame(index=data.columns)) 
    
    ax = data[day:day].plot(kind = 'area', stacked = True, color = colors, figsize = (9,6), legend = 'reverse')
    ax = (n.loads_t.p[day:day].sum(axis=1)).plot(marker ='o', color = '0.5', label = 'load')  #plot load just to check
    plt.fill_between(data[day:day].index, data[day:day].sum(axis=1),n.loads_t.p[day:day].sum(axis=1), where=data[day:day].sum(axis=1)>n.loads_t.p[day:day].sum(axis=1), alpha=0, hatch='\\')
    plt.fill_between(data[day:day].index, data[day:day].sum(axis=1),n.loads_t.p[day:day].sum(axis=1), where=data[day:day].sum(axis=1)<n.loads_t.p[day:day].sum(axis=1), alpha=0, hatch='/')
    charging = mpatches.Patch(alpha = 0, hatch = r'\\\\', label ='charging storages')
    discharging = mpatches.Patch(alpha = 0, hatch = r'////', label ='discharging storages')
    plt.gca().add_artist(plt.legend(handles =[charging, discharging], loc= 'upper left'))
    handles, labels = ax.get_legend_handles_labels()
    labels = new_labels(labels)
    ax.legend(reversed(handles), reversed(labels), loc = 'upper right')
    ax.set_ylim(bottom =0, top=120000)
    ax.set_title('Unit commitment of %d%s'%(i,day[4:]), fontsize = 18)
    ax.set_xlabel('Time', fontsize = 16)
    ax.set_ylabel('Generation in MW', fontsize = 16)
    if i == 2050: plt.savefig('{}/Unit commitment year {}_{}'.format(name,i,minmax), pi=1600, bbox_inches='tight')
    # plt.show()
    
def unit_commitment(n,i,PYROLYSIS):

    
    day_min = '2013-06-02' #day with min electricity demand
    day_max = '2013-12-11' #day with max electricity demand
    #plot unit commitment, simple plotting with pandas
    
    unit_comitment_plot(n, i, day_min, 'min', PYROLYSIS)
    unit_comitment_plot(n, i, day_max, 'max', PYROLYSIS)


def co2_emission(n,i,PYROLYSIS):
    # CO2 sequestration
    if PYROLYSIS == 1:
        co2_seq_sum = n.stores_t.e['co2 atmosphere']*(n.carriers.loc['co2','co2_emissions']) #Alternativ
    else:
        co2_seq_sum = 0

    # CO2 emission by conventionals
    efficiency_by_gen = n.generators.efficiency.groupby(n.generators.carrier).mean()[['oil', 'OCGT', 'CCGT', 'coal', 'lignite']] #efficiency generators
    emi_carriers_gen = n.carriers.co2_emissions[['oil', 'OCGT', 'CCGT', 'coal', 'lignite']] #2030
    emissions_gen = (p_per_carrier(n, PYROLYSIS)[['oil', 'OCGT', 'CCGT', 'coal', 'lignite']]/efficiency_by_gen) *emi_carriers_gen   #emissions by snapshot each generator
    total_emi_gen = emissions_gen.sum(axis=1)       #total emissions (sum of all generators)
    accu_emi_gen = emissions_gen.cumsum(axis=0)  #accumulated emissions from generators
    total_accu_emi_gen = accu_emi_gen.sum(axis=1)   #accumulated emissions generators
    total_accu_emi_gen = pd.DataFrame(total_accu_emi_gen)
    total_accu_emi_gen.columns= ["emissions"] # co2 emission in tonnes
    
    # CO2 total
    total_co2 = total_accu_emi_gen.emissions - co2_seq_sum
    
    # Plotting ------------------------------------
    # fig, ax = plt.subplots(3,1, figsize=(10,10), dpi = 500, sharex = True)
    # ax[0].fill_between(total_accu_emi_gen.index, total_accu_emi_gen.emissions)
    # ax[0].plot(total_accu_emi_gen.index, total_accu_emi_gen.emissions)
    # ax[0].set_ylabel("CO2eq. in t")
    # ax[0].set_title("Accumulative emissions generators in year {}".format(i))
    # ax[1].plot(co2_seq_sum.index, co2_seq_sum)
    # ax[1].fill_between(co2_seq_sum.index, co2_seq_sum)
    # ax[1].set_title("CO2 sequestrated")
    # ax[1].set_ylabel(" CO2eq. in t") 
    # ax[2].plot(total_co2, color = 'orange')
    # ax[2].fill_between(total_co2.index, total_co2, color = 'orange')
    # ax[2].set_title("CO2 total balance")
    # ax[2].set_ylabel(" CO2eq. in t") 
    # ax[2].set_xlabel("Time") 
    # if i == 2050: plt.savefig('{}/CO2 emission in year {}'.format(name,i), pi=1600, bbox_inches='tight')
    # plt.show()
    
    return total_co2[-1]
    
    # print('Mt Co2 emission')
    # print((total_emi_gen.sum()/10**6).round(4))
    # print('Mt Co2 sequestration')
    # print((co2_seq_p.sum()/10**6).round(4))
    # print('Mt Co2 in emission - sequestration')
    # print(((total_emi_gen.sum()+co2_seq_p.sum())/10**6).round(4))
    
    # return ((total_emi_gen.sum()+co2_seq_p.sum())/10**6).round(4)
    
#  ---------------------------------------------



