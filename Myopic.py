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
import pypsa
import data
import re
network_name=data.network_name

print('enter the regional potential')
regional_potential=int(input())
import pyrolysis 

def update_load(n,factor):
    for i in range(len(n.loads_t.p_set.columns)):
        n.loads_t.p_set[n.loads_t.p_set.columns[i]]*=factor

def update_cost(n,year,cost_factors,fuel_cost):
    
    for tech in cost_factors.columns:
        n.generators.capital_cost[n.generators.carrier==tech]*= cost_factors[tech].loc[year]
        n.generators.marginal_cost[n.generators.carrier==tech]*= cost_factors[tech].loc[year]
        n.storage_units.capital_cost[n.storage_units.carrier==tech]*= cost_factors[tech].loc[year]
        n.storage_units.marginal_cost[n.storage_units.carrier==tech]*= cost_factors[tech].loc[year]

    
    fuels=list(fuel_cost.columns)

    for tech in fuels:
            print('correcting {} marginal cost at {}'.format(tech,year))
            eff=n.generators.efficiency[n.generators.carrier==tech].unique().item()
            old_cost=fuel_cost.loc[year-1,tech]/eff
            new_cost=fuel_cost.loc[year,tech]/eff
            n.generators.marginal_cost[n.generators.carrier==tech]+= new_cost-old_cost
    


def update_co2price(n,year,co2price):
    dif= co2price.loc[year].item() - co2price.loc[year-1].item()
    for carrier in list(n.carriers.index[n.carriers.co2_emissions >0]):
        for gen in list(n.generators.index[n.generators.carrier==carrier]):
            print(gen)
            val=dif* n.carriers.loc[carrier].co2_emissions / n.generators.loc[gen,'efficiency']
            n.generators.loc[gen,'marginal_cost']+=val

def update_co2limit(n, new_lim):

    n.global_constraints.constant=float(new_lim)

def delete_old_gens(n,year,base):    
    wanted=base[{'bus','p_nom'}][base.year_removed==year].sum(axis=1)
    for i in range(len(wanted.index)):
        if wanted.index[i] in n.generators.index:
            if wanted.index[i].split()[-1] in ['biomass','ror']:
                n.generators.loc[n.generators.index == wanted.index[i],'p_nom_max'] =+ wanted[i]
                val= n.generators.loc['Fixed ' + wanted.index[i],'p_nom'] - wanted[i]
                if val >= 0 :
                    n.generators.loc['Fixed ' + wanted.index[i],'p_nom']= val
                else:
                    n.generators.loc['Fixed ' + wanted.index[i],'p_nom'] = 0
            else:

                val= n.generators.loc[wanted.index[i], 'p_nom'] - wanted[i]
                if val >= 0 :               
                    n.generators.loc[wanted.index[i], 'p_nom']= val
                else:
                    n.generators.loc[wanted.index[i], 'p_nom'] = 0



def Phase_out(n,carrier, phase_year):
    a=n.generators[n.generators.carrier==carrier]
    Total=a.p_nom.sum()
    Yearly = Total / (phase_year - 2020)
    dist = []
    for i in a.index:
        dist.append( Yearly* a.p_nom.loc[i] / Total)
    return a, Yearly


def remove_Phase_out (n,removal, yearly_value):
    for i in removal.index:
        remove=yearly_value* removal.p_nom.loc[i] / removal.p_nom.sum()
        val= n.generators.loc[i, 'p_nom'] - remove
        if val >= 1 :
            n.generators.loc[[i], 'p_nom']= val
        else:
            n.generators.loc[[i], 'p_nom'] = 0


def update_const_gens(n):      
    for i in range(len(n.generators.index)):
        if n.generators.index[i][:5] =='Fixed':
            n.generators.p_nom[i]+=n.generators.p_nom_opt[n.generators.index[i][6:]]
            n.generators.p_nom_opt[i]=n.generators.p_nom[i]
            if n.generators.index[i].split()[-1] not in ['CCGT', 'OCGT']:
                if n.generators.index[i].split()[-1] in ['ror','biomass']:
                    print(n.generators.index[i])
                    n.generators.p_nom_max[n.generators.index[i][6:]] -= n.generators.p_nom_opt[n.generators.index[i][6:]]
                    
                else:                        
                    n.generators.p_nom_max[n.generators.index[i]]-=n.generators.p_nom_opt[n.generators.index[i][6:]]
                if n.generators.p_nom_max[n.generators.index[i][6:]] < 0:
                    n.generators.loc[n.generators.index == n.generators.index[i][6:], 'p_nom_max']=0
    return


def update_const_lines(n): 
   
    for line in n.lines.index:
        if n.lines.loc[line,'s_nom_opt']>n.lines.loc[line,'s_nom']:
            print(line)
            n.lines.loc[line,'s_nom']=n.lines.loc[line,'s_nom_opt']
    for links in n.links.index:
        if n.links.loc[links,'p_nom_opt']>n.links.loc[links,'p_nom'] and 'pyrolysis' not in links: 
            print(links)
            n.links.loc[links,'p_nom']=n.links.loc[links,'p_nom_opt']   
    



def Yearly_potential(n,saved_potential,regional_potential):        
    for i in n.generators.index[n.generators.p_nom_extendable==True]:
        if i.split()[-1] not in ['biomass','ror']:
            print(i)
            saved_potential[i]-=n.generators.p_nom_opt[i]
            if saved_potential[i] >= regional_potential:
                n.generators.loc[i,'p_nom_max']=regional_potential
            else:
                n.generators.loc[i,'p_nom_max']=saved_potential[i]
    for i in saved_potential.index:
        if saved_potential[i] <= 1:
            n.generators.loc[i,'p_nom_max']=0
            saved_potential[i]=0
    return saved_potential
    

def append_gens(n,year,df):
    idx=[]
    bus=[]
    opt_p=[]
    life=[]
    carrier = [] 
    lifetime=pd.DataFrame({'carrier':['CCGT','OCGT','offwind-ac',
                                    'offwind-dc','onwind','solar','ror','biomass'],
                         'life':[30,30,25,25,25,25,80,30]})    
    for i in n.generators[n.generators.p_nom_extendable==True].index:
        if i.split()[-1] not in ['biomass','ror']:
            print(i)
            idx.append(i)
            bus.append(n.generators.loc[i,'bus'])
            opt_p.append(n.generators.loc[i,'p_nom_opt'])
            life.append(year + lifetime.loc[lifetime.carrier==idx[-1].split()[-1], 'life'].item())
            carrier.append(n.generators.loc[i,'carrier'])  
    temp=pd.DataFrame({'bus':bus,'p_nom':opt_p,'year_added':[year]*len(opt_p),'year_removed':life, 'carrier': carrier})  #Änderung AS
    temp.index=idx
    df = pd.concat([df, temp], ignore_index=False)
    return df


def delete_gens(n,year,df,saved_potential):             
    wanted=df[{'bus','p_nom'}][df.year_removed==year]
    wanted=wanted.groupby(level=0).sum()
    for i in range(len(wanted.index)):
        if 'biomass' not in wanted.index[2]:
            n.generators.loc['Fixed ' + wanted.index[i], 'p_nom']-=wanted.loc[wanted.index[i],'p_nom']
            n.generators.loc['Fixed ' + wanted.index[i], 'p_nom_max']+=wanted.loc[wanted.index[i],'p_nom']
            saved_potential[wanted.index[i]]+=wanted.loc[wanted.index[i],'p_nom']
            if n.generators.loc[wanted.index[i], 'p_nom_max'] + wanted.loc[wanted.index[i],'p_nom'] < regional_potential:
                n.generators.loc[wanted.index[i], 'p_nom_max']+=wanted.loc[wanted.index[i],'p_nom']
            if n.generators.loc['Fixed ' + wanted.index[i], 'p_nom']<=0:
                n.generators.loc['Fixed ' + wanted.index[i], 'p_nom'] = 0


def delete_original_RES(n,year,df,saved_potential,regional_potential):     
    wanted=df[{'bus','p_nom'}][df.year_removed==year]
    for i in range(len(wanted.index)):
        n.generators.loc['Fixed ' + wanted.index[i], 'p_nom']-=wanted.loc[wanted.index[i],'p_nom']
        n.generators.loc['Fixed ' + wanted.index[i], 'p_nom_max']+=wanted.loc[wanted.index[i],'p_nom']
        saved_potential[wanted.index[i]]+=wanted.loc[wanted.index[i],'p_nom']
        if n.generators.loc[wanted.index[i], 'p_nom_max'] + wanted.loc[wanted.index[i],'p_nom'] < regional_potential:
            n.generators.loc[wanted.index[i], 'p_nom_max']+=wanted.loc[wanted.index[i],'p_nom']
            
        if n.generators.loc['Fixed ' + wanted.index[i], 'p_nom']<=0:
            n.generators.loc['Fixed ' + wanted.index[i], 'p_nom'] = 0


 
def remove_phased_out (n):
    for i in n.generators.index[n.generators.p_nom_extendable==False]:
        if 'Fixed ' not in i:
            print(i)
            if n.generators.p_nom[i] == 0:
                n.remove('Generator', i)

def Gen_Bar(n,gen_bar,year,PYROLYSIS): 
    p_by_carrier = n.generators_t.p.sum()
	
    if PYROLYSIS==1: p_by_carrier=pyrolysis.Gen_Bar_pyro(n, p_by_carrier)
    
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
        gen_bar.loc[year,new_gen.index[p]]=new_gen.loc[new_gen.index[p]][0]/10**6

    return gen_bar
    
def Inst_Bar(n,bar,year, PYROLYSIS,elec_eff_pyro):
    
    df=pd.DataFrame({'p_nom':n.generators.p_nom[n.generators.carrier!='load'],'carrier':n.generators.carrier[n.generators.carrier!='load']})
    
    if PYROLYSIS==1: df=pyrolysis.Inst_Bar_pyro(n,df,elec_eff_pyro)   
    
    summation=df.groupby('carrier').sum()
    
    for p in range(len(summation.index)):
        bar.loc[year,summation.index[p]]=summation.loc[summation.index[p]][0]/10**3

    bar=bar.sort_index(axis=0,ascending=True)

    return bar



