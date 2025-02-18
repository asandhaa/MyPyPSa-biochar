# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 13:59:23 2023

@author: asandhaa

- Functions with the same names and suffix _pyro

"""

import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pypsa.linopf import (get_var, define_constraints, linexpr, join_exprs,
                          network_lopf, ilopf)


"""--------
MODEL
---------"""
def biomass_nutsx(n,name,config, cluster_number): 
    nuts_pyrolysis = pd.read_excel(f'1_Regionality/biomass_distribution_{cluster_number}_clusters.xlsx', index_col=0, header=0)
    biomass_straw = config['not-electricity']['biochar']['straw_potential']
    straw_absolut = nuts_pyrolysis['straw_share']*biomass_straw
    nuts_pyrolysis['straw_absolut'] = straw_absolut
    
    biomass_wood = config['not-electricity']['biochar']['wood_potential']
    wood_absolut = nuts_pyrolysis['wood_share']*biomass_wood #in MWh
    nuts_pyrolysis['wood_absolut'] = wood_absolut
    
    
    elec_eff_straw = config['not-electricity']['biochar']['electrical_eff_straw']
    elec_eff_wood = config['not-electricity']['biochar']['electrical_eff_wood']
    bc_eff_straw = config['not-electricity']['biochar']['biochar_eff_straw']
    bc_eff_wood = config['not-electricity']['biochar']['biochar_eff_wood']
    co2_straw = config['not-electricity']['biochar']['co2_seq_factor_straw']
    co2_wood = config['not-electricity']['biochar']['co2_seq_factor_wood']
    def efficiency(efficiency_straw, efficiency_wood):
        return ((nuts_pyrolysis.straw_absolut*efficiency_straw + nuts_pyrolysis.wood_absolut*efficiency_wood)/(straw_absolut+wood_absolut)).round(2)
    nuts_pyrolysis['elec_eff'] = efficiency(elec_eff_straw, elec_eff_wood)
    nuts_pyrolysis['bc_eff'] = efficiency(bc_eff_straw, bc_eff_wood)
    nuts_pyrolysis['co2_seq'] = efficiency(co2_straw, co2_wood)*nuts_pyrolysis['bc_eff']*-1
    nuts_pyrolysis.to_excel(f'{name}/nuts_pyrolysis_{cluster_number}.xlsx')
    biomass_absolut = nuts_pyrolysis.straw_absolut + nuts_pyrolysis.wood_absolut
    n.stores.loc[n.stores.index.str.contains('pyrolysis biomass storage'), 'e_initial'] = biomass_absolut.values
    n.stores.loc[n.stores.index.str.contains('pyrolysis biomass storage'), 'e_nom'] = biomass_absolut.values #AS new 20.12.24
    n.links.loc[n.links.index.str.contains('pyrolysis'), 'efficiency'] = nuts_pyrolysis['elec_eff'].values
    n.links.loc[n.links.index.str.contains('pyrolysis'), 'efficiency2'] = nuts_pyrolysis['bc_eff'].values
    n.links.loc[n.links.index.str.contains('pyrolysis'), 'efficiency3'] = nuts_pyrolysis['co2_seq'].values
    
    return nuts_pyrolysis
"""--------
BASE
---------"""
"""In the original convert_opt_to_conv the existing power plants are set by setting the values for p_nom fixed
For pyrolysis we assume there are zero existing plants, therefore the first p_nom_max = regional_potential, but keep the option if we need to include the 50 in Germany. This is then
done here as:
SOLL-Zustand	p_nom_max   	p_nom  	p_nom_opt
pyrolyse	    8500	        0	     0
Fixed pyrolyse  inf (default)   0        0

"""
def convert_opt_to_conv_pyro(n, regional_potential_pyrolysis):
    for i in range(len(n.links.index)):
        if 'Fixed' not in n.links.index[i] and 'pyrolysis' in n.links.index[i]:
            n.links.loc[i,'p_nom_max'] = regional_potential_pyrolysis


def create_fixed_pyro(n):
            
    for i in n.links.index:
        if 'pyrolysis' in i:
            print(i)
            a=n.links.loc[n.links.index == str(i)]
            n.add("Link",name='Fixed {}'.format(i), bus0=a.bus0[0],bus1=a.bus1[0],carrier= a.carrier[0],
            length=a.length[0], p_nom_extendable=False, capital_cost = 0, p_nom=a.p_nom[0])
            n.links.loc[n.links.index==str(i),'p_nom_opt']=0
            n.links.loc[n.links.index==str(i),'p_nom']=0
            n.links.loc[n.links.index=='Fixed '+ str(i),'capital_cost_lc']=0
            n.links.loc[n.links.index=='Fixed '+ str(i),'capital_cost']=0
            n.links.loc[n.links.index=='Fixed '+ str(i),'carrier']='DC'
            n.links.loc[n.links.index=='Fixed '+ str(i),'geometry']=a.geometry[0]
            n.links.loc[n.links.index=='Fixed '+ str(i),'tags']=a.tags[0]
            n.links.loc[n.links.index=='Fixed '+ str(i),'p_min_pu']=a.p_min_pu[0]
            n.links.loc[n.links.index=='Fixed '+ str(i),'type']=a.type[0]
            n.links.loc[n.links.index=='Fixed '+ str(i),'efficiency']=a.efficiency[0]
            n.links.loc[n.links.index=='Fixed '+ str(i),'p_nom_min']=0
            
            n.links.loc[n.links.index=='Fixed '+ str(i),'marginal_cost']=a.marginal_cost[0]                 
            n.links.loc[n.links.index=='Fixed '+ str(i),'bus2']=a.bus2[0]                                  
            n.links.loc[n.links.index=='Fixed '+ str(i),'bus3']=a.bus3[0]                                  
            n.links.loc[n.links.index=='Fixed '+ str(i),'efficiency2']=a.efficiency2[0]                    
            n.links.loc[n.links.index=='Fixed '+ str(i),'efficiency3']=a.efficiency3[0]                     
        
    for i in n.links.index:
        if 'Fixed' not in i and 'pyrolysis' in i:
            print(i)
            n.links.loc[i,'p_nom_extendable']=True
        if 'Fixed' in i and 'pyrolysis' in i:                 
            n.links.loc[i,'carrier'] = 'pyrolysis'
            print(i)
            
    


"""--------
MYOPIC
---------"""


"""
In the originial update_cost it is done for ['OCGT', 'CCGT', 'biomass', 'coal', 'lignite', 'oil'] and are fixed.
However the marginal cost for pyrolysis are adapted via a yearly factors marginal[year]=marginal[year-1]*factor[year]
"""
def update_cost_pyro(n, i, df_capex_pyro):   
    new_capex = df_capex_pyro.loc[i,'costs']
    n.links.capital_cost[n.links.carrier =='pyrolysis'] = new_capex   


def update_const_gens_pyro(n):
    for i in range(len(n.links.index)):
        if 'Fixed' in n.links.index[i] and 'pyrolysis' in n.links.index[i]:
            print(i)
            n.links.p_nom[i]+=n.links.p_nom_opt[n.links.index[i][6:]]
            n.links.p_nom_opt[i]=n.links.p_nom[i]


def Yearly_potential_pyro(n, regional_potential_pyrolysis):  
    """'In loop: The regional potential per cluster is used to limit the installation of each technology per year. 
            For pyrolysis we are not limited by a p_max like solar and wind, because the installed capacitiy is not directly 
            proportional to the available biomass. So at each year the p_nom_max = 80000 is set for the extendable pyrolysis.
            Which is not a limit at all because with 10 Mio. t Biomass you can install only 2,2 GWel for full load"""
    for i in n.links.index:
        if 'Fixed' not in i and 'pyrolysis' in i:
            print(i)
            n.links.loc[i,'p_nom_max'] = regional_potential_pyrolysis
    
    return

def append_gens_pyro(n,year,df,elec_eff_pyro):
    idx=[]
    bus=[]
    opt_p=[]
    life=[]

    lifetime=pd.DataFrame({'carrier':['pyrolysis'],
                         'life':[100]})                   
    for i in n.links[n.links.p_nom_extendable==True].index:
        if 'pyrolysis' in i:
            print(i)
            idx.append(i)
            bus.append(n.links.loc[i,'bus1'])
            opt_p.append(n.links.loc[i,'p_nom_opt']*elec_eff_pyro)
            life.append(year + lifetime.loc[lifetime.carrier==idx[-1].split()[-1], 'life'].item())
            
    temp=pd.DataFrame({'bus':bus,'p_nom':opt_p,'year_added':[year]*len(opt_p),'year_removed':life})
    temp.index=idx
    df = pd.concat([df, temp], ignore_index=False)
    return df

      



"""--------
PLOTTING
---------"""

def pie_chart_pyro(n, p_by_carrier):
    p_by_pyrolysis = pd.DataFrame(index=n.links_t.p1.columns[n.links_t.p1.columns.str.contains('pyrolysis')], columns = ['p_by_carrier']) #AS
    p_by_pyrolysis.p_by_carrier = n.links_t.p1.sum()*(-1) 
    p_by_carrier = p_by_carrier.append(p_by_pyrolysis) 
    return p_by_carrier


def installed_capacities_pyro(n,df,elec_eff_pyro):
	df_pyrolysis = pd.DataFrame({'p_nom':n.links.p_nom_opt[n.links.p_nom_opt.index.str.contains('pyrolysis')]*elec_eff_pyro,'carrier':'pyrolysis'}) #Änderung AS
	df=df.append(df_pyrolysis)
	return df

def Gen_Bar_pyro(n, p_by_carrier):
	p_by_pyrolysis = n.links_t.p1.sum()[n.links_t.p1.columns.str.contains('pyrolysis')]*(-1)
	p_by_carrier = p_by_carrier.append(p_by_pyrolysis)
	return p_by_carrier
    
def Inst_Bar_pyro(n,df, elec_eff_pyro):
    df_pyro=pd.DataFrame({'p_nom':n.links.p_nom[n.links.carrier=='pyrolysis']*elec_eff_pyro,'carrier':n.links.carrier[n.links.carrier=='pyrolysis']})
    df = df.append(df_pyro)
    return df
    
def Country_Map_pyro(n, bus_sizes, elec_eff_pyro):
    bussizes_pyro = (n.links.query('carrier == "pyrolysis"').groupby(['bus1', 'carrier']).p_nom.sum()*elec_eff_pyro)
    bus_sizes = bus_sizes.append(bussizes_pyro) 
    return bus_sizes


def p_per_carrier_pyrolysis(n,p_per_carrier):
    p_per_carrier = p_per_carrier.rename(columns={'co2':'pyrolysis'})
    p_per_carrier['pyrolysis'] = n.links_t.p1[n.links_t.p1.columns[n.links_t.p1.columns.str.contains('pyrolysis')]].sum(axis=1)*-1
    return p_per_carrier


def biomass_used(n,i, name):
    stores_t=n.stores_t.e[n.stores_t.e.columns[n.stores_t.e.columns.str.contains('pyrolysis')]]
    stores_t.drop(labels=stores_t.index[1:-1], axis=0,inplace=True)

    # Plotting ------------------------------------
    ax = stores_t.plot(kind='bar')
    ax.set_ylabel('MWh')
    ax.set_title('Biomass used in year %d'%i)
    ax.set_xticklabels(['Start', 'End'])
    ax.legend().set_visible(False)
    plt.show()
    if i == 2050: plt.savefig('{}/biomass used {}'.format(name,i), bbox_inches='tight')
    return stores_t.iloc[1,:].sum().round(0)

def flh_pyro(n, elec_eff_pyro):
    generation_pyro = (n.links_t.p1.sum()[n.links_t.p1.columns.str.contains('pyrolysis')]*(-1)).sum()
    installation_pyro = (pd.DataFrame({'p_nom':n.links.p_nom_opt[n.links.p_nom_opt.index.str.contains('pyrolysis')]*elec_eff_pyro,'carrier':'pyrolysis'})).p_nom.sum()
    FLH_pyro = (generation_pyro/installation_pyro).round(0)
    return FLH_pyro


"""--------
SOLVER_NETWORK
---------"""  
def add_CCL_constraints_pyro(n,agg_p_nom_minmax, p_nom_per_cc, year, elec_eff_pyro):
    'The yearly installation of pyrolysis is constrained by yearly values saved in an excel file'
    'The values describe p_nom at the link which is the biomass not the electricity'
    'The values need to be multiplied by the efficiency, then its the maximal installable electrical capacity'

    agg_p_nom_minmax_pyro_year = pd.read_excel('data/agg_p_nom_minmax_pyro.xlsx', index_col=list(range(2)))
    agg_p_nom_minmax_pyro_year['max'] = agg_p_nom_minmax_pyro_year['max']/elec_eff_pyro
    agg_p_nom_minmax_pyro = agg_p_nom_minmax_pyro_year.query(f'year == {year}').drop('year', axis=1)
    
    #Method: Endogen value of yearly limit. For pyrolysis: 70 % of the installed capacity of this year
                # can be installed for the next year. agg_p_nom_minmax_pyro = 0.7 * p_nom
                # Starting value is 8 MWth (1 Mwel), as this is the actual capacity in 2020               
    inst_pyro = round(n.links[n.links.carrier == 'pyrolysis'].p_nom.sum(),0)
    if inst_pyro >= 8:
        agg_p_nom_minmax_pyro['max'] = 0.7*n.links[n.links.carrier == 'pyrolysis'].p_nom.sum()
    else:
        agg_p_nom_minmax_pyro['max'] = 8
    #Method end
    
    agg_p_nom_minmax = agg_p_nom_minmax.append(agg_p_nom_minmax_pyro)
    
    link_country = n.links.bus0.map(n.buses.country)
    1.95090912e+10
    p_nom_per_cc_p = (pd.DataFrame(
                    {'p_nom': linexpr((1, get_var(n, 'Link', 'p_nom'))),
                    'country': link_country, 'carrier': 'pyrolysis'})
                    .dropna(subset=['p_nom'])
                    .groupby(['country', 'carrier']).p_nom
                    .apply(join_exprs))  
    p_nom_per_cc = pd.concat([p_nom_per_cc, p_nom_per_cc_p], ignore_index=False)
    
    return agg_p_nom_minmax, p_nom_per_cc   


