# -*- coding: utf-8 -*-
"""
Created on Sun May  10 02:18:04 2020

@author: anasi
          name      color
0         CCGT        red
1         OCGT     salmon
2      biomass  darkgreen
3         coal     sienna
4      lignite       gray
5      nuclear    fuchsia
6   offwind-ac     indigo
7   offwind-dc       blue
8          oil          k
9       onwind     violet
10         ror       cyan
11       solar     orange
"""


import pandas as pd
import matplotlib.pyplot as plt
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper, FullLoader

import plotting 



def convert_opt_storage_to_conv(n,year):
    idx=[]
    bus=[]
    opt_p=[]
    for i in range(len(n.storage_units[n.storage_units.p_nom_extendable==True].index)):
        idx.append(n.storage_units[n.storage_units.p_nom_extendable==True].index[i])
        bus.append(n.storage_units[n.storage_units.p_nom_extendable==True].bus[i])
        opt_p.append(n.storage_units[n.storage_units.p_nom_extendable==True].p_nom_opt[i])
    df=pd.DataFrame({'bus':bus,'p_nom_opt':opt_p,'year_added':[year]*len(opt_p),'year_removed':[year+15]*len(opt_p)})
    df.index=idx
    
    for i in range(len(df.index)):
        print('Storage Extendable: This is executed for {}'.format(idx[i]))

        g=n.storage_units.loc[n.storage_units.index == idx[i]].copy()
        Value=df.p_nom_opt[idx[i]]
        n.add("StorageUnit","Fixed " + idx[i],
              bus=g.bus[0],p_nom=Value,p_nom_opt=0,marginal_cost=g.marginal_cost[0],
              capital_cost=0, max_hours=g.max_hours[0], carrier=g.carrier[0],p_nom_extendable=False,
              p_nom_max=g.p_nom_max[0],control=g.control[0], p_min_pu=g.p_min_pu[0],
              efficiency_dispatch=g.efficiency_dispatch[0] , efficiency_store=g.efficiency_store[0],
              cyclic_state_of_charge=g.cyclic_state_of_charge[0] )
        
        n.storage_units_t.p_store["Fixed " + idx[i]]=n.storage_units_t.p_store[ idx[i]]*0
        n.storage_units_t.p_dispatch["Fixed " + idx[i]]=n.storage_units_t.p_dispatch[ idx[i]]*0


    for i in n.storage_units.index:
        print('Setting p_nom_opt to zero for {}'.format(i))
        n.storage_units.loc[i,'p_nom_opt']=0
    return df

def update_const_storage(n):
    for i in range(len(n.storage_units.index)):
        if n.storage_units.index[i][:5] =='Fixed':
            n.storage_units.p_nom[i]+=n.storage_units.p_nom_opt[n.storage_units.index[i][6:]]
            n.storage_units.p_nom_opt[i]=n.storage_units.p_nom[i]
    return 

def append_storages(n,year,df):
    idx=[]
    bus=[]
    opt_p=[]
    for i in range(len(n.storage_units[n.storage_units.p_nom_extendable==True].index)):
        idx.append(n.storage_units[n.storage_units.p_nom_extendable==True].index[i])
        bus.append(n.storage_units[n.storage_units.p_nom_extendable==True].bus[i])
        opt_p.append(n.storage_units[n.storage_units.p_nom_extendable==True].p_nom_opt[i])
    temp=pd.DataFrame({'bus':bus,'p_nom_opt':opt_p,'year_added':[year]*len(opt_p),'year_removed':[year+15]*len(opt_p)})
    temp.index=idx
    df = pd.concat([df, temp], ignore_index=False)
    return df

def initial_storage(n):   
    # Last state-of-charge is set as state_of_charge_initial
    initial=n.storage_units.state_of_charge_initial.copy()
    last_state=n.storage_units_t.state_of_charge[n.storage_units_t.state_of_charge.index==n.storage_units_t.state_of_charge.index[-1]].copy()
    
    for i in last_state.columns:
        print(i)
        initial.loc[i]=last_state[i][0]

    n.storage_units.state_of_charge_initial=initial

def delete_storage(n,year,df_stor):
    wanted=df_stor[{'bus','p_nom_opt'}][df_stor.year_removed==year]
    for i in range(len(wanted.index)):
        n.storage_units.loc['Fixed ' + wanted.index[i], 'p_nom']-=wanted.loc[wanted.index[i],'p_nom_opt']


def storage_installation(n,store_bar,year):
    stores=pd.DataFrame({'p_nom':n.storage_units.p_nom,'carrier':n.storage_units.carrier})
    stores.p_nom*=(n.storage_units.max_hours/1e3) 
    summation=stores.groupby('carrier').sum()
        
    for p in range(len(summation.index)):
        store_bar.loc[year,summation.index[p]]=summation.loc[summation.index[p]][0] ## in MW

    store_bar=store_bar.sort_index(axis=0,ascending=True)

    return store_bar

def storage_installation_power(n,store_bar,year):
    stores=pd.DataFrame({'p_nom':n.storage_units.p_nom,'carrier':n.storage_units.carrier})
    
    summation=stores.groupby('carrier').sum()
        
    for p in range(len(summation.index)):
        store_bar.loc[year,summation.index[p]]=summation.loc[summation.index[p]][0] ## in MW

    store_bar=store_bar.sort_index(axis=0,ascending=True)

    return store_bar

def Storage_Bar(store_bar,name):
    colors = plotting.colors_map(pd.DataFrame(index=store_bar.columns)) 
    
    ax = store_bar.plot(figsize=(40, 10),kind='bar', stacked=True,color=colors,fontsize=15)
    ax.set_xlabel("Years",fontsize=15)
    ax.set_ylabel("Storage capacity in GWh",fontsize=15) 
    ax.grid()
    labels = plotting.new_labels(store_bar.columns)
    ax.legend(reversed(plt.legend().legendHandles), reversed(labels), loc="upper left", fontsize = 14) 
    plt.savefig('{}/Storage Bar Plot '.format(name), pi=1600, bbox_inches='tight')
    store_bar.to_excel('{}/Storage_Bar.xlsx'.format(name))



def Storage_Bar_Power(store_bar,name):
    colors = plotting.colors_map(pd.DataFrame(index=store_bar.columns))
    
    ax = store_bar.plot(figsize=(40, 10),kind='bar', stacked=True,color=colors,fontsize=15)
    ax.set_xlabel("Years",fontsize=15)
    ax.set_ylabel("Storage installation in GW",fontsize=15) 
    ax.grid()
    labels = plotting.new_labels(store_bar.columns)
    ax.legend(reversed(plt.legend().legendHandles), reversed(labels), loc="upper left", fontsize = 14)
    plt.savefig('{}/Storage Bar Power Plot '.format(name), pi=1600, bbox_inches='tight')
    store_bar.to_excel('{}/Storage_Bar_Power.xlsx'.format(name))

