# -*- coding: utf-8 -*-
"""
Created on Sun May  10 02:18:04 2020

@author: asandhaa
"""
import os
   
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import pypsa
import pandas as pd
import numpy as np
# from prepare_network import set_line_s_max_pu
import yaml
import logging
logger = logging.getLogger(__name__)
from solve_network import (solve_network, prepare_network)

import base
import Myopic 
import plotting
from override_components import override_component_attrs 
import storage 
import data_nuts2 as data


#%% Input part
network_name=data.network_name # network_name = elec_s_37_ecb_lcopt_Co2L-1H-Ep.nc

"""
Upload start network and config
"""

with open(r'config.yaml') as file:
    config= yaml.load(file, Loader=yaml.FullLoader)


tmpdir = config['solving'].get('tmpdir')
opts = config['scenario']['opts'][0].split('-')
solve_opts = config['solving']['options']

S_N=1
scenario='S'+str(int(S_N))

regional_potential=Myopic.regional_potential # regional_potential= 2500
name=network_name[:-3] 
base.createFolder(network_name[:-3])

overrides = override_component_attrs(r'scripts/override_component_attrs') 
n = pypsa.Network(network_name, override_component_attrs=overrides) 
n.generators['p_nom_min']=0  

"""Start to establish Network for 2020"""
data.update_rens_profiles(n,2013, name) 
data.update_availability_profiles(n,name)

co2lims = pd.read_excel('data/co2limits.xlsx') 
cost_factors= pd.read_csv('data/Cost_Factor.csv', index_col=0, header=0) 
fuel_cost= pd.read_csv('data/fuel_cost.csv', index_col=0, header=0) 
co2price = pd.read_csv('data/co2_price.csv', index_col=0)


"""----------------------------------------------------------------------------
Pyrolysis BOX


Pyrolysis plant""" 
elec_eff_pyro = config['not-electricity']['biochar']['electrical_eff'] 
bc_eff_pyro = config['not-electricity']['biochar']['biochar_eff']
regional_potential_pyrolysis = regional_potential/elec_eff_pyro 

if 'co2' in n.carriers.index:
    PYROLYSIS = 1 
    import pyrolysis

    """1) Pyrolysis plant"""
    n.carriers.loc['co2','co2_emissions'] = -2.96

    """2) Biomass""" 
    cluster_number = len(n.buses[n.buses.v_nom == 380])
    nuts1_pyrolysis = pyrolysis.biomass_nutsx(n,name,config, cluster_number)      

    """3) Costs"""
    df_capex_pyro = pd.read_excel('data/capex_pyro.xlsx', index_col=0)
    capex_pyro = df_capex_pyro.loc[2020,'costs']
    opex_pyro = 0
    n.links.loc[n.links.carrier =='pyrolysis','capital_cost'] = capex_pyro
    n.links.loc[n.links.carrier =='pyrolysis','marginal_cost'] = opex_pyro
    n.stores.marginal_cost = 0
else: PYROLYSIS = 0 

"""----------------------------------------------------------------------------"""
#%% Data part

"""
Logger
"""
agg_p_nom_minmax = pd.read_csv(config['electricity'].get('agg_p_nom_limits'), index_col=1)
co2030=(1-(co2lims[co2lims.year==2030].co2limit.values[0]/(460*10**6)))*100
co2040=(1-(co2lims[co2lims.year==2040].co2limit.values[0]/(460*10**6)))*100
co2045=(1-(co2lims[co2lims.year==2045].co2limit.values[0]/(460*10**6)))*100
co2050=(1-(co2lims[co2lims.year==2050].co2limit.values[0]/(460*10**6)))*100

var=['regional potential',
     'CO2-reduction-2030','CO2-reduction-2040',
     'CO2-reduction-2045','CO2-reduction-2050']
var.extend(list(agg_p_nom_minmax.index))

val=[regional_potential,co2030,co2040,co2045,co2050]
val.extend(list(agg_p_nom_minmax['max']))

txt=pd.DataFrame({'val':val}, index=var)
txt.to_excel("{}/Scenario_settings.xlsx".format(name), index = True)

""""""""""""

Bio_data=data.Biomass_data(n,name)

def read_data(n,name):
    if os.path.exists('{}/conventional_basic_removal.csv'.format(name)):
        conventional_base=pd.read_csv('{}/conventional_basic_removal.csv'.format(name),index_col=0)
    else:
        conventional_base=data.Base_Removal_Data(n)

    if os.path.exists('{}/res_basic_removal.csv'.format(name)):
        RES_base_remove=pd.read_csv('{}/res_basic_removal.csv'.format(name),index_col=0)
    else:
        _,RES_base_remove=data.RES_data(n,name)

    if os.path.exists('{}/res_basic_addition.csv'.format(name)):
        RES_base_addition=pd.read_csv('{}/res_basic_addition.csv'.format(name),index_col=0)
    else:
        RES_base_addition,_=data.RES_data(n,name)
    
    return conventional_base, RES_base_remove, RES_base_addition


conventional_base, RES_base_remove, RES_base_addition =read_data(n,network_name[:-3]) 
coal_data=data.Correct_coal(n,network_name[:-3]) 

conventional_base=conventional_base.append(coal_data, ignore_index=False)       
conventional_base=conventional_base.append(Bio_data, ignore_index=False)        
removal_data = conventional_base.append(RES_base_remove, ignore_index=False)    
removal_data.to_csv('{}/All_removal_data.csv'.format(network_name[:-3]))        


renewables=removal_data[removal_data.carrier=='ror']                            
renewables = renewables.append([removal_data[removal_data.carrier=='solar'],
                                removal_data[removal_data.carrier=='onwind'],
                                removal_data[removal_data.carrier=='offwind-dc'],
                                removal_data[removal_data.carrier=='offwind-ac'],], ignore_index=False)
renewables.drop(renewables[renewables.carrier == 'ror'].index, inplace=True)   



saved_potential=n.generators.p_nom_max[n.generators.p_nom_extendable==True] 
phase_out_removal,yearly_phase_out=Myopic.Phase_out(n,'coal',2036)
phase_out_removal_lignite,yearly_phase_out_lignite=Myopic.Phase_out(n,'lignite',2036) 
if PYROLYSIS==1: pyrolysis.create_fixed_pyro(n) 
df=base.convert_opt_to_conv(n,2020,RES_base_addition, name) 
df_stor = storage.convert_opt_storage_to_conv(n,2020) 

Myopic.update_load(n,545/(n.loads_t.p_set.sum().sum()/1000000))

n.lines.s_max_pu=1.
Myopic.update_co2limit(n,int(co2lims.co2limit[co2lims.year==2020]))
Myopic.update_co2price(n,year=2020,co2price=co2price)

obj=[]
ren_perc=[]

n.storage_units.state_of_charge_initial=0 



gen_bar= pd.DataFrame(index =list(range(2020,2051)),columns=list(n.generators.carrier.unique()))
if PYROLYSIS == 1: gen_bar['pyrolysis'] = np.nan 
inst_bar= pd.DataFrame(index = list(range(2020,2051)) , columns=list(n.generators.carrier[n.generators.carrier!='load'].unique()))
if PYROLYSIS == 1: inst_bar ['pyrolysis'] = np.nan
store_bar= pd.DataFrame(index =list(range(2020,2051)),columns=list(n.storage_units.carrier.unique())) 
saved_potential = Myopic.Yearly_potential(n,saved_potential, regional_potential) 
if PYROLYSIS==1: pyrolysis.Yearly_potential_pyro(n, regional_potential_pyrolysis)

n = prepare_network(n, solve_opts) 
if not os.path.exists('{}/prepared.nc'.format(name)):
    n.export_to_netcdf('{}/prepared.nc'.format(name)) 

#%% Solving test_back.nc
if os.path.exists('{}/test_back.nc'.format(name)):
    n = pypsa.Network('{}/test_back.nc'.format(name), override_component_attrs=overrides) 
else:
    n = pypsa.Network('{}/prepared.nc'.format(name), override_component_attrs=overrides) 
    n = solve_network(n, config=config, solver_dir=tmpdir,opts=opts, year=2020, PYROLYSIS = PYROLYSIS, elec_eff_pyro =elec_eff_pyro) 
    n.export_to_netcdf('{}/test_back.nc'.format(name))


storage.initial_storage(n) 

n.generators.sign[n.generators.carrier=='load']=1
n.generators.marginal_cost[n.generators.carrier=='load']*=10

years_cols=list(range(2020,2051))
Potentials_over_years=pd.DataFrame({'2020':saved_potential})
Potentials_over_years=Potentials_over_years.reindex(columns=Potentials_over_years.columns.tolist() + years_cols) 


obj.append(n.objective - n.generators_t.p[n.generators.index[n.generators.carrier=='load']].sum().sum()*n.generators.marginal_cost.max())
networks=[]

#%% Myopic 
for i in range(2021, 2051):
    Myopic.remove_Phase_out (n,phase_out_removal, yearly_phase_out)
    Myopic.remove_Phase_out (n,phase_out_removal_lignite, yearly_phase_out_lignite) 
    if i == 2036:
        n.generators.p_nom[n.generators.carrier=='coal']=0
        n.generators.p_nom[n.generators.carrier=='lignite']=0   
    Myopic.update_const_lines(n) 	
    Myopic.update_const_gens(n) 
    if PYROLYSIS==1: pyrolysis.update_const_gens_pyro(n)  
    storage.update_const_storage(n)  

    gen_bar=Myopic.Gen_Bar(n,gen_bar,i-1, PYROLYSIS)  
    inst_bar=Myopic.Inst_Bar(n,inst_bar,i-1, PYROLYSIS,elec_eff_pyro) 
    store_bar=storage.storage_installation(n,store_bar,i-1) 
        
    networks.append(n)    
    n.export_to_netcdf('{}/{}.nc'.format(name,i-1)) 

    #from here the year i
    Myopic.update_cost(n,i,cost_factors,fuel_cost=fuel_cost)
    if PYROLYSIS==1: pyrolysis.update_cost_pyro(n, i, df_capex_pyro)
    Myopic.update_load(n,1.01)
    
    Myopic.delete_gens(n,i,df,saved_potential) 
    Myopic.delete_original_RES(n,i,renewables,saved_potential, regional_potential) 
    Myopic.delete_old_gens(n,i,conventional_base) 

    storage.delete_storage(n, i, df_stor) 
    
    n.lines.s_max_pu=1.
    
    Myopic.update_co2limit(n,int(co2lims.co2limit[co2lims.year==i])) 
    Myopic.update_co2price(n,year=i,co2price=co2price)  
      
    n = solve_network(n, config=config, solver_dir=tmpdir,opts=opts, year=i, PYROLYSIS = PYROLYSIS, elec_eff_pyro = elec_eff_pyro)

    storage.initial_storage(n) 
    saved_potential=Myopic.Yearly_potential(n,saved_potential, regional_potential)
    if PYROLYSIS==1: pyrolysis.Yearly_potential_pyro(n, regional_potential_pyrolysis)  
    
    #Logging
    Potentials_over_years.loc[:,i]=saved_potential
    df=Myopic.append_gens(n,i,df)  
    if PYROLYSIS==1: df=pyrolysis.append_gens_pyro(n,i,df, elec_eff_pyro) 
    df_stor = storage.append_storages(n,i,df_stor) 
    obj.append(n.objective - n.generators_t.p[n.generators.index[n.generators.carrier=='load']].sum().sum()*n.generators.marginal_cost.max())


#%% year = 2050
    
Myopic.update_const_lines(n)
Myopic.update_const_gens(n)
storage.update_const_storage(n) 

n.export_to_netcdf('{}/{}.nc'.format(name,i))
with pd.ExcelWriter("{}/addition.xlsx".format(name)) as writer:
    df.to_excel(writer, index = True, sheet_name='Sheet1') 
    df_stor.to_excel(writer, index = True, sheet_name='Sheet2')

objective=pd.DataFrame({'year':list(range(2020,2051)) , 'objective':obj})
objective.to_excel("{}/objective.xlsx".format(name), index = False)
Potentials_over_years.to_excel("{}/Potentials_over_years.xlsx".format(name), index = True)


gen_bar=Myopic.Gen_Bar(n,gen_bar,i, PYROLYSIS) 
inst_bar=Myopic.Inst_Bar(n,inst_bar,i, PYROLYSIS,elec_eff_pyro)
storage.storage_installation(n,store_bar,i), storage.Storage_Bar(store_bar,name) 
networks.append(n) 

 


