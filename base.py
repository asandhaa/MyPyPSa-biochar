import os
import pandas as pd
from vresutils.costdata import annuity
import data

network_name=data.network_name
name=network_name[:-3]

def annuity(n,r):
    """Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6"""

    if isinstance(r, pd.Series):
        return pd.Series(1/n, index=r.index).where(r == 0, r/(1. - 1./(1.+r)**n))
    elif r > 0:
        return r/(1. - 1./(1.+r)**n)
    else:
        return 1/n
    
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)    




def convert_opt_to_conv(n,year,RES_base_addition,name):   
    idx=[]
    bus=[]
    opt_p=[]
    carrier = [] 
            
        
        
    for i in range(len(n.generators[n.generators.p_nom_extendable==True].index)):
        idx.append(n.generators[n.generators.p_nom_extendable==True].index[i])
        bus.append(n.generators[n.generators.p_nom_extendable==True].bus[i])
        carrier.append(n.generators[n.generators.p_nom_extendable==True].carrier[i])  
        opt_p.append(0)
    df=pd.DataFrame({'bus':bus,'p_nom':opt_p,'year_added':[year]*len(opt_p),'year_removed':[year+25]*len(opt_p),'carrier': carrier}) #Änderung AS: column carrier dazu
    df.index=idx
    
    c=pd.read_csv('{}/conventional_basic_removal.csv'.format(name), index_col=0)
    c=c[['carrier','p_nom','bus']]
    c.drop(c[c.carrier =='biomass'].index, inplace=True)
    c.drop(c[c.carrier =='coal'].index, inplace=True)
    c.drop(c[c.carrier =='lignite'].index, inplace=True)
    c.drop(c[c.carrier =='ror'].index, inplace=True)
    c.drop(c[c.carrier =='oil'].index, inplace=True)
    c=c.groupby(['carrier','bus']).agg({'p_nom':['sum']})
    c=c.stack().reset_index()
    c=c[['carrier','p_nom','bus']]
    c.index=c.bus + ' ' + c.carrier
    c.to_csv('{}/extendable_base_addition.csv'.format(name))
    p_max=0
    Value=0
    for i in range(len(df.index)):
        if idx[i] in n.generators_t.p_max_pu.columns:
            print('Renewable Extendable: This is executed for {}'.format(idx[i]))
            y2=list(n.generators_t.p_max_pu[idx[i]])
    
            g=n.generators.loc[n.generators.index == idx[i]].copy()
            if g.carrier[0] != 'ror':
                if idx[i] in RES_base_addition.index:
                    Value=RES_base_addition.p_nom[idx[i]]
                    p_max=g.p_nom_max[0]
                    if p_max <0:
                        p_max=0
                else:
                    Value=0
#                y1=y2*Value
                if Value >=p_max:
                    Value=p_max
                n.generators.loc[idx[i], 'p_nom_max']-= Value
                n.add("Generator","Fixed " + idx[i],
                      bus=g.bus[0],p_nom=Value,p_nom_opt=0,marginal_cost=g.marginal_cost[0],
                      capital_cost=0, carrier=g.carrier[0],p_nom_extendable=False,
                      p_nom_max=p_max,control=g.control[0],
                      efficiency=g.efficiency[0], p_min_pu=0, p_max_pu=y2)     
                n.generators.loc["Fixed " + idx[i],'weight']=g.weight[0]
                n.generators_t.p["Fixed " + idx[i]]=0
        else:
            print('Conventional Extendable: This is executed for {}'.format(idx[i]))
            g=n.generators.loc[n.generators.index == idx[i]].copy()
            if idx[i] in c.index:
                Value=c.p_nom[idx[i]]
            else:
                Value=0

#            y1=y2*Value
            n.add("Generator","Fixed " + idx[i],
                  bus=g.bus[0],p_nom=Value,p_nom_opt=0,marginal_cost=g.marginal_cost[0],
                  capital_cost=0, carrier=g.carrier[0],p_nom_extendable=False,control='',
                  p_nom_max=g.p_nom_max[0],
                  efficiency=g.efficiency[0])
            
            n.generators.loc["Fixed " + idx[i],'weight']=g.weight[0]
            n.generators_t.p["Fixed " + idx[i]]=0

    for i in n.generators.index:
        print('Setting p_nom_opt to zero for {}'.format(i))
        n.generators.loc[i,'p_nom_opt']=0
    
    
    # Biomass & ROR
    cap_bio= ((annuity(30, 0.07) +
                             3.6/100.) *
                             2350*1e3 * 1)
    cap_ror= ((annuity(80, 0.07) +
                             2/100.) *
                             2500*1e3 * 1)

    for idx in n.generators.index[(n.generators.carrier=='biomass') | (n.generators.carrier=='ror')]:
        print(idx)
        n.generators.loc[n.generators.index == idx,'capital_cost']= cap_bio if idx.split()[-1] == 'biomass' else cap_ror
        g=n.generators.loc[n.generators.index == idx].copy()
        n.add("Generator","Fixed " + idx,
              bus=g.bus[0],p_nom=g.p_nom[0],p_nom_opt=0,marginal_cost=g.marginal_cost[0],
              capital_cost=0, carrier=g.carrier[0],p_nom_extendable=False,
              p_nom_max=0,efficiency=g.efficiency[0])
        n.generators_t.p["Fixed " + idx]=0
        n.generators.loc[n.generators.index == idx,'p_nom_max']=0
        n.generators.loc[n.generators.index == idx,'p_nom_extendable']=True
        n.generators.loc[n.generators.index == idx,'p_nom']=0
        n.generators.loc[n.generators.index == idx,'p_nom_opt']=0

        if g.carrier[0] == 'ror':
            y2=list(n.generators_t.p_max_pu[g.index[0]])
            n.generators_t.p_max_pu['Fixed '  + idx] = y2
    return df
