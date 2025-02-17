# -*- coding: utf-8 -*-
"""
@author: asandhaa
"""

## snakemake
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from os.path import normpath, exists
from shutil import copyfile
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")

configfile: "config.yaml"

COSTS="data/costs.csv"
ATLITE_NPROCESSES = config['atlite'].get('nprocesses', 4)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*"


rule cluster_all_networks:
    input: expand("networks/elec_s{simpl}_{clusters}.nc", **config['scenario'])

#modified
rule extra_components_all_networks:
    input: expand("networks/elec_s{simpl}_{clusters}_ecb.nc", **config['scenario'])

#modified
rule prepare_all_networks:
    input: expand("networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc", **config['scenario'])

#modified
rule solve_all_networks:
    input: expand("results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc", **config['scenario'])


if config['enable'].get('prepare_links_p_nom', False):
    rule prepare_links_p_nom:
        output: 'data/links_p_nom.csv'
        log: 'logs/prepare_links_p_nom.log'
        threads: 1
        resources: mem=500
        script: 'scripts/prepare_links_p_nom.py'


datafiles = ['ch_cantons.csv', 'je-e-21.03.02.xls',
            'eez/World_EEZ_v8_2014.shp', 'EIA_hydro_generation_2000_2014.csv',
            'hydro_capacities.csv', 'naturalearth/ne_10m_admin_0_countries.shp',
            'NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp', 'nama_10r_3popgdp.tsv.gz',
            'nama_10r_3gdp.tsv.gz', 'corine/g250_clc06_V18_5.tif']


if not config.get('tutorial', False):
    datafiles.extend(["natura/Natura2000_end2015.shp", "GEBCO_2014_2D.nc"])


if config['enable'].get('retrieve_databundle', True):
    rule retrieve_databundle:
        output: expand('data/bundle/{file}', file=datafiles)
        log: "logs/retrieve_databundle.log"
        script: 'scripts/retrieve_databundle.py'


rule retrieve_load_data:
    input: HTTP.remote("data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv", keep_local=True, static=True)
    output: "data/load_raw.csv"
    shell: "mv {input} {output}"


rule build_load_data:
    input: "data/load_raw.csv"
    output: "resources/load.csv"
    log: "logs/build_load_data.log"
    script: 'scripts/build_load_data.py'


rule build_powerplants:
    input:
        base_network="networks/base.nc",
        custom_powerplants="data/custom_powerplants.csv"
    output: "resources/powerplants.csv"
    log: "logs/build_powerplants.log"
    threads: 1
    resources: mem=500
    script: "scripts/build_powerplants.py"


rule base_network:
    input:
        eg_buses='data/entsoegridkit/buses.csv',
        eg_lines='data/entsoegridkit/lines.csv',
        eg_links='data/entsoegridkit/links.csv',
        eg_converters='data/entsoegridkit/converters.csv',
        eg_transformers='data/entsoegridkit/transformers.csv',
        parameter_corrections='data/parameter_corrections.yaml',
        links_p_nom='data/links_p_nom.csv',
        links_tyndp='data/links_tyndp.csv',
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        europe_shape='resources/europe_shape.geojson'
    output: "networks/base.nc"
    log: "logs/base_network.log"
    benchmark: "benchmarks/base_network"
    threads: 1
    resources: mem=500
    script: "scripts/base_network.py"


rule build_shapes:
    input:
        naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
        eez='data/bundle/eez/World_EEZ_v8_2014.shp',
        nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
        nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
        nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
        ch_cantons='data/bundle/ch_cantons.csv',
        ch_popgdp='data/bundle/je-e-21.03.02.xls'
    output:
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        europe_shape='resources/europe_shape.geojson',
        nuts3_shapes='resources/nuts3_shapes.geojson'
    log: "logs/build_shapes.log"
    threads: 1
    resources: mem=500
    script: "scripts/build_shapes.py"


rule build_bus_regions:
    input:
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        base_network="networks/base.nc"
    output:
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    log: "logs/build_bus_regions.log"
    threads: 1
    resources: mem=1000
    script: "scripts/build_bus_regions.py"

if config['enable'].get('build_cutout', False):
    rule build_cutout:
        input:
            regions_onshore="resources/regions_onshore.geojson",
            regions_offshore="resources/regions_offshore.geojson"
        output: "cutouts/{cutout}.nc"
        log: "logs/build_cutout/{cutout}.log"
        benchmark: "benchmarks/build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources: mem=ATLITE_NPROCESSES * 1000
        script: "scripts/build_cutout.py"


if config['enable'].get('retrieve_cutout', True):
    rule retrieve_cutout:
        input: HTTP.remote("zenodo.org/record/4709858/files/{cutout}.nc", keep_local=True, static=True)
        output: "cutouts/{cutout}.nc"
        shell: "mv {input} {output}"


if config['enable'].get('build_natura_raster', False):
    rule build_natura_raster:
        input:
            natura="data/bundle/natura/Natura2000_end2015.shp",
            cutouts=expand("cutouts/{cutouts}.nc", **config['atlite'])
        output: "resources/natura.tiff"
        log: "logs/build_natura_raster.log"
        script: "scripts/build_natura_raster.py"


if config['enable'].get('retrieve_natura_raster', True):
    rule retrieve_natura_raster:
        input: HTTP.remote("zenodo.org/record/4706686/files/natura.tiff", keep_local=True, static=True)
        output: "resources/natura.tiff"
        shell: "mv {input} {output}"


rule build_renewable_profiles:
    input:
        base_network="networks/base.nc",
        corine="data/bundle/corine/g250_clc06_V18_5.tif",
        natura="resources/natura.tiff",
        gebco=lambda w: ("data/bundle/GEBCO_2014_2D.nc"
                         if "max_depth" in config["renewable"][w.technology].keys()
                         else []),
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        regions=lambda w: ("resources/regions_onshore.geojson"
                                   if w.technology in ('onwind', 'solar')
                                   else "resources/regions_offshore.geojson"),
        cutout=lambda w: "cutouts/" + config["renewable"][w.technology]['cutout'] + ".nc"
    output: profile="resources/profile_{technology}.nc",
    log: "logs/build_renewable_profile_{technology}.log"
    benchmark: "benchmarks/build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources: mem=ATLITE_NPROCESSES * 5000
    script: "scripts/build_renewable_profiles.py"


if 'hydro' in config['renewable'].keys():
    rule build_hydro_profile:
        input:
            country_shapes='resources/country_shapes.geojson',
            eia_hydro_generation='data/bundle/EIA_hydro_generation_2000_2014.csv',
            cutout="cutouts/" + config["renewable"]['hydro']['cutout'] + ".nc"
        output: 'resources/profile_hydro.nc'
        log: "logs/build_hydro_profile.log"
        resources: mem=5000
        script: 'scripts/build_hydro_profile.py'


rule add_electricity:
    input:
        base_network='networks/base.nc',
        tech_costs=COSTS,
        regions="resources/regions_onshore.geojson",
        powerplants='resources/powerplants.csv',
        hydro_capacities='data/bundle/hydro_capacities.csv',
        geth_hydro_capacities='data/geth2015_hydro_capacities.csv',
        load='resources/load.csv',
        nuts3_shapes='resources/nuts3_shapes.geojson',
        **{f"profile_{tech}": f"resources/profile_{tech}.nc"
           for tech in config['renewable']}
    output: "networks/elec.nc",
        costs_updated = "resources/costs.csv",
    log: "logs/add_electricity.log"
    benchmark: "benchmarks/add_electricity"
    threads: 1
    resources: mem=5000
    script: "scripts/add_electricity.py"


rule simplify_network:
    input:
        network='networks/elec.nc',
        tech_costs=COSTS,
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    output:
        network='networks/elec_s{simpl}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}.geojson",
        busmap='resources/busmap_elec_s{simpl}.csv',
        connection_costs='resources/connection_costs_s{simpl}.csv'
    log: "logs/simplify_network/elec_s{simpl}.log"
    benchmark: "benchmarks/simplify_network/elec_s{simpl}"
    threads: 1
    resources: mem=4000
    script: "scripts/simplify_network.py"


rule cluster_network:
    input:
        network='networks/elec_s{simpl}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}.geojson",
        busmap=ancient('resources/busmap_elec_s{simpl}.csv'),
        custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
                       if config["enable"].get("custom_busmap", False) else []),
        tech_costs=COSTS
    output:
        network='networks/elec_s{simpl}_{clusters}.nc',
        regions_onshore="resources/regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/regions_offshore_elec_s{simpl}_{clusters}.geojson",
        busmap="resources/busmap_elec_s{simpl}_{clusters}.csv",
        linemap="resources/linemap_elec_s{simpl}_{clusters}.csv"
    log: "logs/cluster_network/elec_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/cluster_network/elec_s{simpl}_{clusters}"
    threads: 1
    resources: mem=6000
    script: "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        tech_costs=COSTS,
    output: 'networks/elec_s{simpl}_{clusters}_ec.nc'
    log: "logs/add_extra_components/elec_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources: mem=3000
    script: "scripts/add_extra_components.py"

#created
rule add_biochar:
    input:
        network='networks/elec_s{simpl}_{clusters}_ec.nc',
        overrides="data/override_component_attrs",  #file with the overrides for the links
        #tech_costs=COSTS,
    output:'networks/elec_s{simpl}_{clusters}_ecb.nc'
    log: "logs/add_biochar/elec_s{simpl}_{clusters}_ecb.log"
    benchmark: "benchmarks/add_biochar/elec_s{simpl}_{clusters}_ecb"
    threads: 1
    resources: mem=3000
    script: "scripts/add_biochar.py"

#modified
rule prepare_network:
    input:
        network='networks/elec_s{simpl}_{clusters}_ecb.nc', tech_costs=COSTS,
        overrides="data/override_component_attrs",  #file with the overrides for the links
    output: 'networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc'
    log: "logs/prepare_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.log"
    benchmark: "benchmarks/prepare_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}"
    threads: 1
    resources: mem=4000
    script: "scripts/prepare_network.py"


def memory(w):
    factor = 3.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)seg$', o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))

#modified
rule solve_network:
    input:
        network="networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc",
        overrides="data/override_component_attrs",  #file with the overrides for the links
    output: "results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc"
    log:
        solver=normpath("logs/solve_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_solver.log"),
        python="logs/solve_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_python.log",
        memory="logs/solve_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_memory.log"
    benchmark: "benchmarks/solve_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}"
    threads: 4
    resources: mem=memory
    shadow: "shallow"
    script: "scripts/solve_network.py"

#modified
rule solve_operations_network:
    input:
        unprepared="networks/elec_s{simpl}_{clusters}_ecb.nc",
        optimized="results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc"
    output: "results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_op.nc"
    log:
        solver=normpath("logs/solve_operations_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_op_solver.log"),
        python="logs/solve_operations_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_op_python.log",
        memory="logs/solve_operations_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_op_memory.log"
    benchmark: "benchmarks/solve_operations_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}"
    threads: 4
    resources: mem=(lambda w: 5000 + 372 * int(w.clusters))
    shadow: "shallow"
    script: "scripts/solve_operations_network.py"

#modified
rule plot_network:
    input:
        network="results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc",
        tech_costs=COSTS
    output:
        only_map="results/plots/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{attr}.{ext}",
        ext="results/plots/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{attr}_ext.{ext}"
    log: "logs/plot_network/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{attr}_{ext}.log"
    script: "scripts/plot_network.py"

#modified
def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return ([COSTS] +
            expand("results/networks/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}.nc",
                   ll=ll,
                   **{k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
                      for k in ["simpl", "clusters", "opts"]}))

#modified
rule make_summary:
    input: input_make_summary
    output: directory("results/summaries/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{country}")
    log: "logs/make_summary/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{country}.log",
    script: "scripts/make_summary.py"

#modified
rule plot_summary:
    input: "results/summaries/elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{country}"
    output: "results/plots/summary_{summary}_elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{country}.{ext}"
    log: "logs/plot_summary/{summary}_elec_s{simpl}_{clusters}_ecb_l{ll}_{opts}_{country}_{ext}.log"
    script: "scripts/plot_summary.py"


def input_plot_p_nom_max(w):
    return [("networks/elec_s{simpl}{maybe_cluster}.nc"
             .format(maybe_cluster=('' if c == 'full' else ('_' + c)), **w))
            for c in w.clusts.split(",")]


rule plot_p_nom_max:
    input: input_plot_p_nom_max
    output: "results/plots/elec_s{simpl}_cum_p_nom_max_{clusts}_{techs}_{country}.{ext}"
    log: "logs/plot_p_nom_max/elec_s{simpl}_{clusts}_{techs}_{country}_{ext}.log"
    script: "scripts/plot_p_nom_max.py"

## config file

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: CC0-1.0

version: 0.4.0
tutorial: false

logging:
  level: INFO
  format: '%(levelname)s:%(name)s:%(message)s'

summary_dir: results

scenario:
  simpl: ['']
  ll: ['copt']
  clusters: [16]
  opts: [Co2L-3H]

countries: ['DE']

clustering:
  simplify:
    to_substations: false # network is simplified to nodes with positive or negative power injection (i.e. substations or offwind connections)

snapshots:
  start: "2013-01-01"
  end: "2014-01-01"
  closed: 'left' # end is not inclusive

enable:
  prepare_links_p_nom: false
  retrieve_databundle: true
  build_cutout: false
  retrieve_cutout: true
  build_natura_raster: false
  retrieve_natura_raster: true
  custom_busmap: false

#new (NOT YET IMPLEMENTED. Only for option: delete generators + replace with store and links)
#pypsa_eur:
#  Bus:
#    - AC
#  Link:
#    - DC
#  Generator:
#    - onwind
#    - offwind-ac
#    - offwind-dc
#    - solar
#    - ror
#  StorageUnit:
#    - PHS
#    - hydro
#  Store: []

# new
not-electricity:
  biochar:
    Store: [co2-tracking] #co2-tracking
    Link: [pyrolysis] #pyrolysis
    co2_seq: -0.21473   #t Co2eq/MWh
    biomass_potential: 7.5e+8 #MWh/year
    capital_cost: 1300    #€/MW
    marginal_cost: 26.38   #€/MWh
    electrial_eff: 0.32   #%
    biochar_eff: 0.36    #%
    lifetime: 20    #years


electricity:
  voltages: [220., 300., 380.]
  co2limit: 7.75e+7 # 0.05 * 3.1e9*0.5
  co2base: 1.487e+9
  agg_p_nom_limits: data/agg_p_nom_minmax.csv

  extendable_carriers:
    Generator: []
    StorageUnit: [] # battery, H2
    Store: [battery, H2]
    Link: []

  max_hours:
    battery: 6
    H2: 168

  powerplants_filter: false # use pandas query strings here, e.g. Country not in ['Germany']
  custom_powerplants: false # use pandas query strings here, e.g. Country in ['Germany']
  conventional_carriers: [oil, OCGT, CCGT, coal, lignite, geothermal, biomass] #nuclear, oil, OCGT, CCGT, coal, lignite, geothermal, biomass
  renewable_capacities_from_OPSD: [] # onwind, offwind, solar

  # estimate_renewable_capacities_from_capacity_stats:
  #   # Wind is the Fueltype in ppm.data.Capacity_stats, onwind, offwind-{ac,dc} the carrier in PyPSA-Eur
  #   Wind: [onwind, offwind-ac, offwind-dc]
  #   Solar: [solar]

atlite:
  nprocesses: 4
  cutouts:
    # use 'base' to determine geographical bounds and time span from config
    # base:
      # module: era5
    europe-2013-era5:
      module: era5 # in priority order
      x: [-12., 35.]
      y: [33., 72]
      dx: 0.3
      dy: 0.3
      time: ['2013', '2013']
    europe-2013-sarah:
      module: [sarah, era5] # in priority order
      x: [-12., 45.]
      y: [33., 65]
      dx: 0.2
      dy: 0.2
      time: ['2013', '2013']
      sarah_interpolate: false
      sarah_dir:
      features: [influx, temperature]


renewable:
  onwind:
    cutout: europe-2013-era5
    resource:
      method: wind
      turbine: Vestas_V112_3MW
    capacity_per_sqkm: 3 # ScholzPhd Tab 4.3.1: 10MW/km^2
    # correction_factor: 0.93
    corine:
      # Scholz, Y. (2012). Renewable energy based electricity supply at low costs:
      #  development of the REMix model and application for Europe. ( p.42 / p.28)
      grid_codes: [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                   24, 25, 26, 27, 28, 29, 31, 32]
      distance: 1000
      distance_grid_codes: [1, 2, 3, 4, 5, 6]
    natura: true
    potential: simple # or conservative
    clip_p_max_pu: 1.e-2
  offwind-ac:
    cutout: europe-2013-era5
    resource:
      method: wind
      turbine: NREL_ReferenceTurbine_5MW_offshore
    capacity_per_sqkm: 2
    correction_factor: 0.8855
    # proxy for wake losses
    # from 10.1016/j.energy.2018.08.153
    # until done more rigorously in #153
    corine: [44, 255]
    natura: true
    max_depth: 50
    max_shore_distance: 30000
    potential: simple # or conservative
    clip_p_max_pu: 1.e-2
  offwind-dc:
    cutout: europe-2013-era5
    resource:
      method: wind
      turbine: NREL_ReferenceTurbine_5MW_offshore
    # ScholzPhd Tab 4.3.1: 10MW/km^2
    capacity_per_sqkm: 2
    correction_factor: 0.8855
    # proxy for wake losses
    # from 10.1016/j.energy.2018.08.153
    # until done more rigorously in #153
    corine: [44, 255]
    natura: true
    max_depth: 50
    min_shore_distance: 30000
    potential: simple # or conservative
    clip_p_max_pu: 1.e-2
  solar:
    cutout: europe-2013-sarah
    resource:
      method: pv
      panel: CSi
      orientation:
        slope: 35.
        azimuth: 180.
    capacity_per_sqkm: 1.7 # ScholzPhd Tab 4.3.1: 170 MW/km^2
    # Correction factor determined by comparing uncorrected area-weighted full-load hours to those
    # published in Supplementary Data to
    # Pietzcker, Robert Carl, et al. "Using the sun to decarbonize the power
    # sector: The economic potential of photovoltaics and concentrating solar
    # power." Applied Energy 135 (2014): 704-720.
    # This correction factor of 0.854337 may be in order if using reanalysis data.
    # for discussion refer to https://github.com/PyPSA/pypsa-eur/pull/304
    # correction_factor: 0.854337
    corine: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
             14, 15, 16, 17, 18, 19, 20, 26, 31, 32]
    natura: true
    potential: simple # or conservative
    clip_p_max_pu: 1.e-2
  hydro:
    cutout: europe-2013-era5
    carriers: [ror, PHS, hydro]
    PHS_max_hours: 6
    hydro_max_hours: "energy_capacity_totals_by_country" # one of energy_capacity_totals_by_country, estimate_by_large_installations or a float
    clip_min_inflow: 1.0

lines:
  types:
    220.: "Al/St 240/40 2-bundle 220.0"
    300.: "Al/St 240/40 3-bundle 300.0"
    380.: "Al/St 240/40 4-bundle 380.0"
  s_max_pu: 0.7
  s_nom_max: .inf
  length_factor: 1.25
  under_construction: 'zero' # 'zero': set capacity to zero, 'remove': remove, 'keep': with full capacity

links:
  p_max_pu: 1.0
  p_nom_max: .inf
  include_tyndp: true
  under_construction: 'zero' # 'zero': set capacity to zero, 'remove': remove, 'keep': with full capacity

transformers:
  x: 0.1
  s_nom: 2000.
  type: ''

load:
  url: https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv
  power_statistics: True # only for files from <2019; set false in order to get ENTSOE transparency data
  interpolate_limit: 3 # data gaps up until this size are interpolated linearly
  time_shift_for_large_gaps: 1w # data gaps up until this size are copied by copying from
  manual_adjustments: true # false
  scaling_factor: 1.0

costs:
  year: 2030
  discountrate: 0.07 # From a Lion Hirth paper, also reflects average of Noothout et al 2016
  USD2013_to_EUR2013: 0.7532 # [EUR/USD] ECB: https://www.ecb.europa.eu/stats/exchange/eurofxref/html/eurofxref-graph-usd.en.html
  marginal_cost: # EUR/MWh
    solar: 0.01
    onwind: 0.015
    offwind: 0.015
    hydro: 0.
    H2: 0.
    electrolysis: 0.
    fuel cell: 0.
    battery: 0.
    battery inverter: 0.
  emission_prices: # in currency per tonne emission, only used with the option Ep
    co2: 0.

solving:
  options:
    formulation: kirchhoff
    load_shedding: false
    noisy_costs: true
    min_iterations: 4
    max_iterations: 6
    clip_p_max_pu: 0.01
    skip_iterations: false
    track_iterations: false
    #nhours: 10
  solver:
    name: gurobi
    threads: 4
    method: 2 # barrier
    crossover: 0
    BarConvTol: 1.e-5
    FeasibilityTol: 1.e-6
    AggFill: 0
    PreDual: 0
    GURO_PAR_BARDENSETHRESH: 200
  # solver:
  #   name: cplex
  #   threads: 4
  #   lpmethod: 4 # barrier
  #   solutiontype: 2 # non basic solution, ie no crossover
  #   barrier.convergetol: 1.e-5
  #   feasopt.tolerance: 1.e-6

plotting:
  map:
    figsize: [7, 7]
    boundaries: [-10.2, 29, 35,  72]
    p_nom:
      bus_size_factor: 5.e+4
      linewidth_factor: 3.e+3

  costs_max: 800
  costs_threshold: 1

  energy_max: 15000.
  energy_min: -10000.
  energy_threshold: 50.

  vre_techs: ["onwind", "offwind-ac", "offwind-dc", "solar", "ror"]
  conv_techs: ["OCGT", "CCGT", "Nuclear", "Coal"]
  storage_techs: ["hydro+PHS", "battery", "H2"]
  load_carriers: ["AC load"]
  AC_carriers: ["AC line", "AC transformer"]
  link_carriers: ["DC line", "Converter AC-DC"]
  tech_colors:
    "onwind" : "#235ebc"
    "onshore wind" : "#235ebc"
    'offwind' : "#6895dd"
    'offwind-ac' : "#6895dd"
    'offshore wind' : "#6895dd"
    'offshore wind ac' : "#6895dd"
    'offwind-dc' : "#74c6f2"
    'offshore wind dc' : "#74c6f2"
    "hydro" : "#08ad97"
    "hydro+PHS" : "#08ad97"
    "PHS" : "#08ad97"
    "hydro reservoir" : "#08ad97"
    'hydroelectricity' : '#08ad97'
    "ror" : "#4adbc8"
    "run of river" : "#4adbc8"
    'solar' : "#f9d002"
    'solar PV' : "#f9d002"
    'solar thermal' : '#ffef60'
    'biomass' : '#0c6013'
    'solid biomass' : '#06540d'
    'biogas' : '#23932d'
    'waste' : '#68896b'
    'geothermal' : '#ba91b1'
    "OCGT" : "#d35050"
    "gas" : "#d35050"
    "natural gas" : "#d35050"
    "CCGT" : "#b20101"
    "nuclear" : "#ff9000"
    "coal" : "#707070"
    "lignite" : "#9e5a01"
    "oil" : "#262626"
    "H2" : "#ea048a"
    "hydrogen storage" : "#ea048a"
    "battery" : "#b8ea04"
    "Electric load" : "#f9d002"
    "electricity" : "#f9d002"
    "lines" : "#70af1d"
    "transmission lines" : "#70af1d"
    "AC-AC" : "#70af1d"
    "AC line" : "#70af1d"
    "links" : "#8a1caf"
    "HVDC links" : "#8a1caf"
    "DC-DC" : "#8a1caf"
    "DC link" : "#8a1caf"
  nice_names:
    OCGT: "Open-Cycle Gas"
    CCGT: "Combined-Cycle Gas"
    offwind-ac: "Offshore Wind (AC)"
    offwind-dc: "Offshore Wind (DC)"
    onwind: "Onshore Wind"
    solar: "Solar"
    PHS: "Pumped Hydro Storage"
    hydro: "Reservoir & Dam"
    battery: "Battery Storage"
    H2: "Hydrogen Storage"
    lines: "Transmission Lines"
    ror: "Run of River"




## add_ biochar
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 17:13:14 2022

@author: mdomenech


Biochar script

step 1. add co2 tracking: co2 atmosphere store which tracks the co2 emitted to the atmosphere
step 2. add co2 stored: tracks the co2 stored (captured) in this case this is equivalent to the amount of biochar
step 3. add override components to the networks (for the links)
step 4. add pyrolysis link from the biomass store to the electrical buses and from the pyrolysis to the co2 stored (biochar)


"""
#librariesfrom _helpers python
import logging  #to log messages that you want to see
from _helpers import configure_logging
from override_components import override_component_attrs

#from add_electricity import (add_nice_carrier_names,
#                             _add_missing_carriers_from_costs)


import pypsa
import pandas as pd
import numpy as np




def add_co2_tracking(n):

    elec_opts = snakemake.config['not-electricity']
    carriers = elec_opts['biochar']['Store']

    #_add_missing_carriers_from_costs(n, carriers)

    #buses_i = n.buses.index
    #bus_sub_dict = {k: n.buses[k].values for k in ['x', 'y', 'country']}

    if 'co2-tracking' in carriers:
    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
        n.add("Carrier", "co2",

              co2_emissions= snakemake.config["not-electricity"]["biochar"]["co2_seq"],  #sequestration potential biochar
              nice_name = "CO2",
              color = "#235ebc")

        # this tracks CO2 in the atmosphere
        n.add("Bus",
            "co2 atmosphere",
            location="DE",
            carrier="co2")


        # can also be negative
        n.add("Store", "co2 atmosphere",
            #e_nom = 1000,
            e_nom_extendable=True,
            e_min_pu=-1,
            carrier="co2",
            bus="co2 atmosphere",
            #capital_cost=0,
        )

        #this tracks co2 stored (biochar in this case)
        n.add("Bus","co2 stored",
              location= "DE",
              carrier= "co2 stored")

        n.add("Store","co2 stored",
              #e_nom = 1000,
              #e_min_pu=-1,
              e_nom_extendable=True,
              #e_nom_max=np.inf,
              #capital_cost=0.,  #capital cost for biochar store
              carrier="co2 stored",
              bus="co2 stored"
              )


def add_pyrolysis(n):

    print("adding pyrolysis")

    elec_opts = snakemake.config['not-electricity']
    tech = elec_opts['biochar']['Link']


    buses_i = n.buses.index[n.buses.carrier == "AC"] #this way co2 tracking is not considered
    j = n.buses[n.buses.carrier == "AC"]
    j = j.shape[0]
    bus_sub_dict =  {k: n.buses[k].values[0:j,] for k in ['x', 'y', 'country']} #do not select co2 tracking buses n.buses.values


    if 'pyrolysis' in tech:
        print("adding pyrolysis")

        pyrolysis_buses_i = n.madd("Bus", buses_i + " pyrolysis", carrier="pyrolysis", **bus_sub_dict)#bus_sub_dict)   #n.buses.index[n.buses.carrier == "AC"]

        #this is the wood store
        n.madd("Store", pyrolysis_buses_i + " biomass storage",
               bus= pyrolysis_buses_i,
               carrier="pyrolysis",
               e_nom = snakemake.config["not-electricity"]["biochar"]["biomass_potential"],  #capacity pyrolysis (available biomass)
               e_initial = snakemake.config["not-electricity"]["biochar"]["biomass_potential"],
               #e_cyclic=True,
               #e_nom_extendable=True,
               capital_cost = 0.) #capital cost biomass store


        n.madd("Link",pyrolysis_buses_i, #this is what is written in name
               bus0= pyrolysis_buses_i,
               bus1= buses_i,
               bus2="co2 stored",
               bus3="co2 atmosphere",
               #carrier= 'pyrolysis biochar',
               p_nom_extendable=True,
               #lifetime = snakemake.config["not-electricity"]["biochar"]["lifetime"],,
               efficiency = snakemake.config["not-electricity"]["biochar"]["electrical_eff"],  #electrical efficiency
               efficiency2 = snakemake.config["not-electricity"]["biochar"]["biochar_eff"],    #biochar efficiency
               efficiency3 = -snakemake.config["not-electricity"]["biochar"]["biochar_eff"], #co2 reduction from atmosphere
               capital_cost = snakemake.config["not-electricity"]["biochar"]["capital_cost"],  #1300000.0,    #55250.0, #capital cost pyrolysis
               marginal_cost = snakemake.config["not-electricity"]["biochar"]["marginal_cost"] #9.5  #26.38 #marginal cost pyrolysis
               )

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('add_biochar', network='elec',
                                   simpl='', clusters=5)

    configure_logging(snakemake)

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    add_co2_tracking(n)
    add_pyrolysis(n)

    n.export_to_netcdf(snakemake.output[0])

## add_electricity
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Adds electrical generators and existing hydro storage units to a base network.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        USD2013_to_EUR2013:
        dicountrate:
        emission_prices:

    electricity:
        max_hours:
        marginal_cost:
        capital_cost:
        conventional_carriers:
        co2limit:
        extendable_carriers:
        include_renewable_capacities_from_OPSD:
        estimate_renewable_capacities_from_capacity_stats:

    load:
        scaling_factor:

    renewable:
        hydro:
            carriers:
            hydro_max_hours:
            hydro_capital_cost:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`, :ref:`load_cf`, :ref:`renewable_cf`, :ref:`lines_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``data/bundle/hydro_capacities.csv``: Hydropower plant store/discharge power capacities, energy storage capacity, and average hourly inflow by country.

    .. image:: ../img/hydrocapacities.png
        :scale: 34 %

- ``data/geth2015_hydro_capacities.csv``: alternative to capacities above; not currently used!
- ``resources/opsd_load.csv`` Hourly per-country load profiles.
- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/nuts3_shapes.geojson``: confer :ref:`shapes`
- ``resources/powerplants.csv``: confer :ref:`powerplants`
- ``resources/profile_{}.nc``: all technologies in ``config["renewables"].keys()``, confer :ref:`renewableprofiles`.
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``networks/elec.nc``:

    .. image:: ../img/elec.png
            :scale: 33 %

Description
-----------

The rule :mod:`add_electricity` ties all the different data inputs from the preceding rules together into a detailed PyPSA network that is stored in ``networks/elec.nc``. It includes:

- today's transmission topology and transfer capacities (optionally including lines which are under construction according to the config settings ``lines: under_construction`` and ``links: under_construction``),
- today's thermal and hydro power generation capacities (for the technologies listed in the config setting ``electricity: conventional_carriers``), and
- today's load time-series (upsampled in a top-down approach according to population and gross domestic product)

It further adds extendable ``generators`` with **zero** capacity for

- photovoltaic, onshore and AC- as well as DC-connected offshore wind installations with today's locational, hourly wind and solar capacity factors (but **no** current capacities),
- additional open- and combined-cycle gas turbines (if ``OCGT`` and/or ``CCGT`` is listed in the config setting ``electricity: extendable_carriers``)
"""

import logging
from _helpers import configure_logging, update_p_nom_max

import pypsa
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import powerplantmatching as pm
from powerplantmatching.export import map_country_bus

from vresutils.costdata import annuity
from vresutils.load import timeseries_opsd
from vresutils import transfer as vtransfer

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(s): return s/s.sum()


def _add_missing_carriers_from_costs(n, costs, carriers):
    missing_carriers = pd.Index(carriers).difference(n.carriers.index)
    if missing_carriers.empty: return

    emissions_cols = costs.columns.to_series()\
                           .loc[lambda s: s.str.endswith('_emissions')].values
    suptechs = missing_carriers.str.split('-').str[0]
    emissions = costs.loc[suptechs, emissions_cols].fillna(0.)
    emissions.index = missing_carriers
    n.import_components_from_dataframe(emissions, 'Carrier')


def load_costs(Nyears=1., tech_costs=None, config=None, elec_config=None):
    if tech_costs is None:
        tech_costs = snakemake.input.tech_costs

    if config is None:
        config = snakemake.config['costs']

    # set all asset costs and other parameters
    costs = pd.read_csv(tech_costs, index_col=list(range(3))).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"),"value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"),"value"] *= config['USD2013_to_EUR2013']

    costs = (costs.loc[idx[:,config['year'],:], "value"]
             .unstack(level=2).groupby("technology").sum(min_count=1))

    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                           "discount rate" : config['discountrate'],
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : 25})

    costs["capital_cost"] = ((annuity(costs["lifetime"], costs["discount rate"]) +
                             costs["FOM"]/100.) *
                             costs["investment"] * Nyears)

    costs.at['OCGT', 'fuel'] = costs.at['gas', 'fuel']
    costs.at['CCGT', 'fuel'] = costs.at['gas', 'fuel']

    costs['marginal_cost'] = costs['VOM'] + costs['fuel'] / costs['efficiency']

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})

    costs.at['OCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']
    costs.at['CCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']

    costs.at['solar', 'capital_cost'] = 0.5*(costs.at['solar-rooftop', 'capital_cost'] +
                                             costs.at['solar-utility', 'capital_cost'])

    def costs_for_storage(store, link1, link2=None, max_hours=1.):
        capital_cost = link1['capital_cost'] + max_hours * store['capital_cost']
        if link2 is not None:
            capital_cost += link2['capital_cost']
        return pd.Series(dict(capital_cost=capital_cost,
                              marginal_cost=0.,
                              co2_emissions=0.))

    if elec_config is None:
        elec_config = snakemake.config['electricity']
    max_hours = elec_config['max_hours']
    costs.loc["battery"] = \
        costs_for_storage(costs.loc["battery storage"], costs.loc["battery inverter"],
                          max_hours=max_hours['battery'])
    costs.loc["H2"] = \
        costs_for_storage(costs.loc["hydrogen storage"], costs.loc["fuel cell"],
                          costs.loc["electrolysis"], max_hours=max_hours['H2'])

    for attr in ('marginal_cost', 'capital_cost'):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites

    return costs


def load_powerplants(ppl_fn=None):
    if ppl_fn is None:
        ppl_fn = snakemake.input.powerplants
    carrier_dict = {'ocgt': 'OCGT', 'ccgt': 'CCGT', 'bioenergy': 'biomass',
                    'ccgt, thermal': 'CCGT', 'hard coal': 'coal'}
    return (pd.read_csv(ppl_fn, index_col=0, dtype={'bus': 'str'})
            .powerplant.to_pypsa_names()
            .rename(columns=str.lower).drop(columns=['efficiency'])
            .replace({'carrier': carrier_dict}))


def attach_load(n):
    substation_lv_i = n.buses.index[n.buses['substation_lv']]
    regions = (gpd.read_file(snakemake.input.regions).set_index('name')
               .reindex(substation_lv_i))
    opsd_load = (pd.read_csv(snakemake.input.load, index_col=0, parse_dates=True)
                .filter(items=snakemake.config['countries']))

    scaling = snakemake.config.get('load', {}).get('scaling_factor', 1.0)
    logger.info(f"Load data scaled with scalling factor {scaling}.")
    opsd_load *= scaling

    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index('index')

    def upsample(cntry, group):
        l = opsd_load[cntry]
        if len(group) == 1:
            return pd.DataFrame({group.index[0]: l})
        else:
            nuts3_cntry = nuts3.loc[nuts3.country == cntry]
            transfer = vtransfer.Shapes2Shapes(group, nuts3_cntry.geometry,
                                               normed=False).T.tocsr()
            gdp_n = pd.Series(transfer.dot(nuts3_cntry['gdp'].fillna(1.).values),
                              index=group.index)
            pop_n = pd.Series(transfer.dot(nuts3_cntry['pop'].fillna(1.).values),
                              index=group.index)

            # relative factors 0.6 and 0.4 have been determined from a linear
            # regression on the country to continent load data
            # (refer to vresutils.load._upsampling_weights)
            factors = normed(0.6 * normed(gdp_n) + 0.4 * normed(pop_n))
            return pd.DataFrame(factors.values * l.values[:,np.newaxis],
                                index=l.index, columns=factors.index)

    load = pd.concat([upsample(cntry, group) for cntry, group
                      in regions.geometry.groupby(regions.country)], axis=1)

    n.madd("Load", substation_lv_i, bus=substation_lv_i, p_set=load)


def update_transmission_costs(n, costs, length_factor=1.0, simple_hvdc_costs=False):
    n.lines['capital_cost'] = (n.lines['length'] * length_factor *
                               costs.at['HVAC overhead', 'capital_cost'])

    if n.links.empty: return

    dc_b = n.links.carrier == 'DC'

    # If there are no dc links, then the 'underwater_fraction' column
    # may be missing. Therefore we have to return here.
    if n.links.loc[dc_b].empty: return

    if simple_hvdc_costs:
        costs = (n.links.loc[dc_b, 'length'] * length_factor *
                 costs.at['HVDC overhead', 'capital_cost'])
    else:
        costs = (n.links.loc[dc_b, 'length'] * length_factor *
                ((1. - n.links.loc[dc_b, 'underwater_fraction']) *
                costs.at['HVDC overhead', 'capital_cost'] +
                n.links.loc[dc_b, 'underwater_fraction'] *
                costs.at['HVDC submarine', 'capital_cost']) +
                costs.at['HVDC inverter pair', 'capital_cost'])
    n.links.loc[dc_b, 'capital_cost'] = costs


def attach_wind_and_solar(n, costs):
    for tech in snakemake.config['renewable']:
        if tech == 'hydro': continue

        n.add("Carrier", name=tech)
        with xr.open_dataset(getattr(snakemake.input, 'profile_' + tech)) as ds:
            if ds.indexes['bus'].empty: continue

            suptech = tech.split('-', 2)[0]
            if suptech == 'offwind':
                underwater_fraction = ds['underwater_fraction'].to_pandas()
                connection_cost = (snakemake.config['lines']['length_factor'] *
                                   ds['average_distance'].to_pandas() *
                                   (underwater_fraction *
                                    costs.at[tech + '-connection-submarine', 'capital_cost'] +
                                    (1. - underwater_fraction) *
                                    costs.at[tech + '-connection-underground', 'capital_cost']))
                capital_cost = (costs.at['offwind', 'capital_cost'] +
                                costs.at[tech + '-station', 'capital_cost'] +
                                connection_cost)
                logger.info("Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}"
                            .format(connection_cost.min(), connection_cost.max(), tech))
            else:
                capital_cost = costs.at[tech, 'capital_cost']

            n.madd("Generator", ds.indexes['bus'], ' ' + tech,
                   bus=ds.indexes['bus'],
                   carrier=tech,
                   p_nom_extendable=True,
                   p_nom_max=ds['p_nom_max'].to_pandas(),
                   weight=ds['weight'].to_pandas(),
                   marginal_cost=costs.at[suptech, 'marginal_cost'],
                   capital_cost=capital_cost,
                   efficiency=costs.at[suptech, 'efficiency'],
                   p_max_pu=ds['profile'].transpose('time', 'bus').to_pandas())


def attach_conventional_generators(n, costs, ppl):
    carriers = snakemake.config['electricity']['conventional_carriers']

    _add_missing_carriers_from_costs(n, costs, carriers)

    ppl = (ppl.query('carrier in @carriers').join(costs, on='carrier')
           .rename(index=lambda s: 'C' + str(s)))

    logger.info('Adding {} generators with capacities [MW] \n{}'
                .format(len(ppl), ppl.groupby('carrier').p_nom.sum()))

    n.madd("Generator", ppl.index,
           carrier=ppl.carrier,
           bus=ppl.bus,
           p_nom=ppl.p_nom,
           efficiency=ppl.efficiency,
           marginal_cost=ppl.marginal_cost,
           capital_cost=0)

    logger.warning(f'Capital costs for conventional generators put to 0 EUR/MW.')


def attach_hydro(n, costs, ppl):
    if 'hydro' not in snakemake.config['renewable']: return
    c = snakemake.config['renewable']['hydro']
    carriers = c.get('carriers', ['ror', 'PHS', 'hydro'])

    _add_missing_carriers_from_costs(n, costs, carriers)

    ppl = ppl.query('carrier == "hydro"').reset_index(drop=True)\
             .rename(index=lambda s: str(s) + ' hydro')
    ror = ppl.query('technology == "Run-Of-River"')
    phs = ppl.query('technology == "Pumped Storage"')
    hydro = ppl.query('technology == "Reservoir"')

    country = ppl['bus'].map(n.buses.country).rename("country")

    inflow_idx = ror.index.union(hydro.index)
    if not inflow_idx.empty:
        dist_key = ppl.loc[inflow_idx, 'p_nom'].groupby(country).transform(normed)

        with xr.open_dataarray(snakemake.input.profile_hydro) as inflow:
            inflow_countries = pd.Index(country[inflow_idx])
            missing_c = (inflow_countries.unique()
                         .difference(inflow.indexes['countries']))
            assert missing_c.empty, (f"'{snakemake.input.profile_hydro}' is missing "
                f"inflow time-series for at least one country: {', '.join(missing_c)}")

            inflow_t = (inflow.sel(countries=inflow_countries)
                        .rename({'countries': 'name'})
                        .assign_coords(name=inflow_idx)
                        .transpose('time', 'name')
                        .to_pandas()
                        .multiply(dist_key, axis=1))

    if 'ror' in carriers and not ror.empty:
        n.madd("Generator", ror.index,
               carrier='ror',
               bus=ror['bus'],
               p_nom=ror['p_nom'],
               efficiency=costs.at['ror', 'efficiency'],
               capital_cost=costs.at['ror', 'capital_cost'],
               weight=ror['p_nom'],
               p_max_pu=(inflow_t[ror.index]
                         .divide(ror['p_nom'], axis=1)
                         .where(lambda df: df<=1., other=1.)))

    if 'PHS' in carriers and not phs.empty:
        # fill missing max hours to config value and
        # assume no natural inflow due to lack of data
        phs = phs.replace({'max_hours': {0: c['PHS_max_hours']}})
        n.madd('StorageUnit', phs.index,
               carrier='PHS',
               bus=phs['bus'],
               p_nom=phs['p_nom'],
               capital_cost=costs.at['PHS', 'capital_cost'],
               max_hours=phs['max_hours'],
               efficiency_store=np.sqrt(costs.at['PHS','efficiency']),
               efficiency_dispatch=np.sqrt(costs.at['PHS','efficiency']),
               cyclic_state_of_charge=True)

    if 'hydro' in carriers and not hydro.empty:
        hydro_max_hours = c.get('hydro_max_hours')
        hydro_stats = pd.read_csv(snakemake.input.hydro_capacities,
                                   comment="#", na_values='-', index_col=0)
        e_target = hydro_stats["E_store[TWh]"].clip(lower=0.2) * 1e6
        e_installed = hydro.eval('p_nom * max_hours').groupby(hydro.country).sum()
        e_missing = e_target - e_installed
        missing_mh_i = hydro.query('max_hours == 0').index

        if hydro_max_hours == 'energy_capacity_totals_by_country':
            # watch out some p_nom values like IE's are totally underrepresented
            max_hours_country = e_missing / \
                                hydro.loc[missing_mh_i].groupby('country').p_nom.sum()

        elif hydro_max_hours == 'estimate_by_large_installations':
            max_hours_country = hydro_stats['E_store[TWh]'] * 1e3 / \
                                hydro_stats['p_nom_discharge[GW]']

        missing_countries = (pd.Index(hydro['country'].unique())
                             .difference(max_hours_country.dropna().index))
        if not missing_countries.empty:
            logger.warning("Assuming max_hours=6 for hydro reservoirs in the countries: {}"
                           .format(", ".join(missing_countries)))
        hydro_max_hours = hydro.max_hours.where(hydro.max_hours > 0,
                                hydro.country.map(max_hours_country)).fillna(6)

        n.madd('StorageUnit', hydro.index, carrier='hydro',
               bus=hydro['bus'],
               p_nom=hydro['p_nom'],
               max_hours=hydro_max_hours,
               capital_cost=(costs.at['hydro', 'capital_cost']
                             if c.get('hydro_capital_cost') else 0.),
               marginal_cost=costs.at['hydro', 'marginal_cost'],
               p_max_pu=1.,  # dispatch
               p_min_pu=0.,  # store
               efficiency_dispatch=costs.at['hydro', 'efficiency'],
               efficiency_store=0.,
               cyclic_state_of_charge=True,
               inflow=inflow_t.loc[:, hydro.index])


def attach_extendable_generators(n, costs, ppl):
    elec_opts = snakemake.config['electricity']
    carriers = pd.Index(elec_opts['extendable_carriers']['Generator'])

    _add_missing_carriers_from_costs(n, costs, carriers)

    for tech in carriers:
        if tech.startswith('OCGT'):
            ocgt = ppl.query("carrier in ['OCGT', 'CCGT']").groupby('bus', as_index=False).first()
            n.madd('Generator', ocgt.index,
                   suffix=' OCGT',
                   bus=ocgt['bus'],
                   carrier=tech,
                   p_nom_extendable=True,
                   p_nom=0.,
                   capital_cost=costs.at['OCGT', 'capital_cost'],
                   marginal_cost=costs.at['OCGT', 'marginal_cost'],
                   efficiency=costs.at['OCGT', 'efficiency'])

        elif tech.startswith('CCGT'):
            ccgt = ppl.query("carrier in ['OCGT', 'CCGT']").groupby('bus', as_index=False).first()
            n.madd('Generator', ccgt.index,
                   suffix=' CCGT',
                   bus=ccgt['bus'],
                   carrier=tech,
                   p_nom_extendable=True,
                   p_nom=0.,
                   capital_cost=costs.at['CCGT', 'capital_cost'],
                   marginal_cost=costs.at['CCGT', 'marginal_cost'],
                   efficiency=costs.at['CCGT', 'efficiency'])

        elif tech.startswith('nuclear'):
            nuclear = ppl.query("carrier == 'nuclear'").groupby('bus', as_index=False).first()
            n.madd('Generator', nuclear.index,
                suffix=' nuclear',
                bus=nuclear['bus'],
                carrier=tech,
                p_nom_extendable=True,
                p_nom=0.,
                capital_cost=costs.at['nuclear', 'capital_cost'],
                marginal_cost=costs.at['nuclear', 'marginal_cost'],
                efficiency=costs.at['nuclear', 'efficiency'])

        else:
            raise NotImplementedError(f"Adding extendable generators for carrier "
                                      "'{tech}' is not implemented, yet. "
                                      "Only OCGT, CCGT and nuclear are allowed at the moment.")



def attach_OPSD_renewables(n):

    available = ['DE', 'FR', 'PL', 'CH', 'DK', 'CZ', 'SE', 'GB']
    tech_map = {'Onshore': 'onwind', 'Offshore': 'offwind', 'Solar': 'solar'}
    countries = set(available) & set(n.buses.country)
    techs = snakemake.config['electricity'].get('renewable_capacities_from_OPSD', [])
    tech_map = {k: v for k, v in tech_map.items() if v in techs}

    if not tech_map:
        return

    logger.info(f'Using OPSD renewable capacities in {", ".join(countries)} '
                f'for technologies {", ".join(tech_map.values())}.')

    df = pd.concat([pm.data.OPSD_VRE_country(c) for c in countries])
    technology_b = ~df.Technology.isin(['Onshore', 'Offshore'])
    df['Fueltype'] = df.Fueltype.where(technology_b, df.Technology)
    df = df.query('Fueltype in @tech_map').powerplant.convert_country_to_alpha2()

    for fueltype, carrier_like in tech_map.items():
        gens = n.generators[lambda df: df.carrier.str.contains(carrier_like)]
        buses = n.buses.loc[gens.bus.unique()]
        gens_per_bus = gens.groupby('bus').p_nom.count()

        caps = map_country_bus(df.query('Fueltype == @fueltype'), buses)
        caps = caps.groupby(['bus']).Capacity.sum()
        caps = caps / gens_per_bus.reindex(caps.index, fill_value=1)

        n.generators.p_nom.update(gens.bus.map(caps).dropna())
        n.generators.p_nom_min.update(gens.bus.map(caps).dropna())



def estimate_renewable_capacities(n, tech_map=None):
    if tech_map is None:
        tech_map = (snakemake.config['electricity']
                    .get('estimate_renewable_capacities_from_capacity_stats', {}))

    if len(tech_map) == 0: return

    capacities = (pm.data.Capacity_stats().powerplant.convert_country_to_alpha2()
                  [lambda df: df.Energy_Source_Level_2]
                  .set_index(['Fueltype', 'Country']).sort_index())

    countries = n.buses.country.unique()

    if len(countries) == 0: return

    logger.info('heuristics applied to distribute renewable capacities [MW] \n{}'
                .format(capacities.query('Fueltype in @tech_map.keys() and Capacity >= 0.1')
                        .groupby('Country').agg({'Capacity': 'sum'})))

    for ppm_fueltype, techs in tech_map.items():
        tech_capacities = capacities.loc[ppm_fueltype, 'Capacity']\
                                    .reindex(countries, fill_value=0.)
        #tech_i = n.generators.query('carrier in @techs').index
        tech_i = (n.generators.query('carrier in @techs')
                  [n.generators.query('carrier in @techs')
                   .bus.map(n.buses.country).isin(countries)].index)
        n.generators.loc[tech_i, 'p_nom'] = (
            (n.generators_t.p_max_pu[tech_i].mean() *
             n.generators.loc[tech_i, 'p_nom_max']) # maximal yearly generation
             .groupby(n.generators.bus.map(n.buses.country))
             .transform(lambda s: normed(s) * tech_capacities.at[s.name])
             .where(lambda s: s>0.1, 0.))  # only capacities above 100kW
        n.generators.loc[tech_i, 'p_nom_min'] = n.generators.loc[tech_i, 'p_nom']


def add_nice_carrier_names(n, config=None):
    if config is None: config = snakemake.config
    carrier_i = n.carriers.index
    nice_names = (pd.Series(config['plotting']['nice_names'])
                  .reindex(carrier_i).fillna(carrier_i.to_series().str.title()))
    n.carriers['nice_name'] = nice_names
    colors = pd.Series(config['plotting']['tech_colors']).reindex(carrier_i)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(f'tech_colors for carriers {missing_i} not defined '
                       'in config.')
    n.carriers['color'] = colors


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('add_electricity')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.

    costs = load_costs(Nyears)
    ppl = load_powerplants()

    attach_load(n)

    update_transmission_costs(n, costs)

    attach_conventional_generators(n, costs, ppl)
    attach_wind_and_solar(n, costs)
    attach_hydro(n, costs, ppl)
    attach_extendable_generators(n, costs, ppl)

    estimate_renewable_capacities(n)
    attach_OPSD_renewables(n)
    update_p_nom_max(n)

    add_nice_carrier_names(n)

    n.export_to_netcdf(snakemake.output[0])
    costs.to_csv(snakemake.output.costs_updated)

#add_extra_components
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Adds extra extendable components to the clustered and simplified network.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        USD2013_to_EUR2013:
        dicountrate:
        emission_prices:

    electricity:
        max_hours:
        marginal_cost:
        capital_cost:
        extendable_carriers:
            StorageUnit:
            Store:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.

Outputs
-------

- ``networks/elec_s{simpl}_{clusters}_ec.nc``:


Description
-----------

The rule :mod:`add_extra_components` attaches additional extendable components to the clustered and simplified network. These can be configured in the ``config.yaml`` at ``electricity: extendable_carriers:``.
It processes ``networks/elec_s{simpl}_{clusters}.nc`` to build ``networks/elec_s{simpl}_{clusters}_ec.nc``, which in contrast to the former (depending on the configuration) contain with **zero** initial capacity

- ``StorageUnits`` of carrier 'H2' and/or 'battery'. If this option is chosen, every bus is given an extendable ``StorageUnit`` of the corresponding carrier. The energy and power capacities are linked through
a parameter that specifies the energy capacity as maximum hours at full dispatch power and is configured in ``electricity: max_hours:``. This linkage leads to one investment variable per storage unit.
The default ``max_hours`` lead to long-term hydrogen and short-term battery storage units.

- ``Stores`` of carrier 'H2' and/or 'battery' in combination with ``Links``. If this option is chosen, the script adds extra buses with corresponding carrier where energy ``Stores`` are attached and which are
connected to the corresponding power buses via two links, one each for charging and discharging. This leads to three investment variables for the energy capacity, charging and discharging capacity of the storage unit.
"""
import logging
from _helpers import configure_logging

import pypsa
import pandas as pd
import numpy as np

from add_electricity import (load_costs, add_nice_carrier_names,
                             _add_missing_carriers_from_costs)

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def attach_storageunits(n, costs):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['StorageUnit']
    max_hours = elec_opts['max_hours']

    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index

    lookup_store = {"H2": "electrolysis", "battery": "battery inverter"}
    lookup_dispatch = {"H2": "fuel cell", "battery": "battery inverter"}

    for carrier in carriers:
        n.madd("StorageUnit", buses_i, ' ' + carrier,
               bus=buses_i,
               carrier=carrier,
               p_nom_extendable=True,
               capital_cost=costs.at[carrier, 'capital_cost'],
               marginal_cost=costs.at[carrier, 'marginal_cost'],
               efficiency_store=costs.at[lookup_store[carrier], 'efficiency'],
               efficiency_dispatch=costs.at[lookup_dispatch[carrier], 'efficiency'],
               max_hours=max_hours[carrier],
               cyclic_state_of_charge=True)


def attach_stores(n, costs):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['Store']

    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index
    bus_sub_dict = {k: n.buses[k].values for k in ['x', 'y', 'country']}

    if 'H2' in carriers:
        h2_buses_i = n.madd("Bus", buses_i + " H2", carrier="H2", **bus_sub_dict)

        n.madd("Store", h2_buses_i,
               bus=h2_buses_i,
               carrier='H2',
               e_nom_extendable=True,
               e_cyclic=True,
               capital_cost=costs.at["hydrogen storage", "capital_cost"])

        n.madd("Link", h2_buses_i + " Electrolysis",
               bus0=buses_i,
               bus1=h2_buses_i,
               carrier='H2 electrolysis',
               p_nom_extendable=True,
               efficiency=costs.at["electrolysis", "efficiency"],
               capital_cost=costs.at["electrolysis", "capital_cost"],
               marginal_cost=costs.at["electrolysis", "marginal_cost"])

        n.madd("Link", h2_buses_i + " Fuel Cell",
               bus0=h2_buses_i,
               bus1=buses_i,
               carrier='H2 fuel cell',
               p_nom_extendable=True,
               efficiency=costs.at["fuel cell", "efficiency"],
               #NB: fixed cost is per MWel
               capital_cost=costs.at["fuel cell", "capital_cost"] * costs.at["fuel cell", "efficiency"],
               marginal_cost=costs.at["fuel cell", "marginal_cost"])

    if 'battery' in carriers:
        b_buses_i = n.madd("Bus", buses_i + " battery", carrier="battery", **bus_sub_dict)

        n.madd("Store", b_buses_i,
               bus=b_buses_i,
               carrier='battery',
               e_cyclic=True,
               e_nom_extendable=True,
               capital_cost=costs.at['battery storage', 'capital_cost'],
               marginal_cost=costs.at["battery", "marginal_cost"])

        n.madd("Link", b_buses_i + " charger",
               bus0=buses_i,
               bus1=b_buses_i,
               carrier='battery charger',
               efficiency=costs.at['battery inverter', 'efficiency'],
               capital_cost=costs.at['battery inverter', 'capital_cost'],
               p_nom_extendable=True,
               marginal_cost=costs.at["battery inverter", "marginal_cost"])

        n.madd("Link", b_buses_i + " discharger",
               bus0=b_buses_i,
               bus1=buses_i,
               carrier='battery discharger',
               efficiency=costs.at['battery inverter','efficiency'],
               p_nom_extendable=True,
               marginal_cost=costs.at["battery inverter", "marginal_cost"])


def attach_hydrogen_pipelines(n, costs):
    elec_opts = snakemake.config['electricity']
    ext_carriers = elec_opts['extendable_carriers']
    as_stores = ext_carriers.get('Store', [])

    if 'H2 pipeline' not in ext_carriers.get('Link',[]): return

    assert 'H2' in as_stores, ("Attaching hydrogen pipelines requires hydrogen "
            "storage to be modelled as Store-Link-Bus combination. See "
            "`config.yaml` at `electricity: extendable_carriers: Store:`.")

    # determine bus pairs
    attrs = ["bus0","bus1","length"]
    candidates = pd.concat([n.lines[attrs], n.links.query('carrier=="DC"')[attrs]])\
                    .reset_index(drop=True)

    # remove bus pair duplicates regardless of order of bus0 and bus1
    h2_links = candidates[~pd.DataFrame(np.sort(candidates[['bus0', 'bus1']])).duplicated()]
    h2_links.index = h2_links.apply(lambda c: f"H2 pipeline {c.bus0}-{c.bus1}", axis=1)

    # add pipelines
    n.madd("Link",
           h2_links.index,
           bus0=h2_links.bus0.values + " H2",
           bus1=h2_links.bus1.values + " H2",
           p_min_pu=-1,
           p_nom_extendable=True,
           length=h2_links.length.values,
           capital_cost=costs.at['H2 pipeline','capital_cost']*h2_links.length,
           efficiency=costs.at['H2 pipeline','efficiency'],
           carrier="H2 pipeline")


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('add_extra_components', network='elec',
                                  simpl='', clusters=5)
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760.
    costs = load_costs(Nyears, tech_costs=snakemake.input.tech_costs,
                       config=snakemake.config['costs'],
                       elec_config=snakemake.config['electricity'])

    attach_storageunits(n, costs)
    attach_stores(n, costs)
    attach_hydrogen_pipelines(n, costs)

    add_nice_carrier_names(n, config=snakemake.config)

    n.export_to_netcdf(snakemake.output[0])

## base network
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Creates the network topology from a `ENTSO-E map extract <https://github.com/PyPSA/GridKit/tree/master/entsoe>`_ (January 2020) as a PyPSA network.

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    countries:

    electricity:
        voltages:

    lines:
        types:
        s_max_pu:
        under_construction:

    links:
        p_max_pu:
        under_construction:
        include_tyndp:

    transformers:
        x:
        s_nom:
        type:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`snapshots_cf`, :ref:`toplevel_cf`, :ref:`electricity_cf`, :ref:`load_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`transformers_cf`

Inputs
------

- ``data/entsoegridkit``:  Extract from the geographical vector data of the online `ENTSO-E Interactive Map <https://www.entsoe.eu/data/map/>`_ by the `GridKit <https://github.com/pypsa/gridkit>`_ toolkit dating back to January 2020.
- ``data/parameter_corrections.yaml``: Corrections for ``data/entsoegridkit``
- ``data/links_p_nom.csv``: confer :ref:`links`
- ``data/links_tyndp.csv``: List of projects in the `TYNDP 2018 <https://tyndp.entsoe.eu/tyndp2018/>`_ that are at least *in permitting* with fields for start- and endpoint (names and coordinates), length, capacity, construction status, and project reference ID.
- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/europe_shape.geojson``: confer :ref:`shapes`

Outputs
-------

- ``networks/base.nc``

    .. image:: ../img/base.png
        :scale: 33 %

Description
-----------

"""

import logging
from _helpers import configure_logging

import pypsa
import yaml
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx

from scipy import spatial
from scipy.sparse import csgraph
from itertools import product

from shapely.geometry import Point, LineString
import shapely, shapely.prepared, shapely.wkt

logger = logging.getLogger(__name__)


def _get_oid(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"oid"=>"(\d+)"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _get_country(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"country"=>"([A-Z]{2})"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _find_closest_links(links, new_links, distance_upper_bound=1.5):
    treecoords = np.asarray([np.asarray(shapely.wkt.loads(s).coords)[[0, -1]].flatten()
                              for s in links.geometry])
    querycoords = np.vstack([new_links[['x1', 'y1', 'x2', 'y2']],
                            new_links[['x2', 'y2', 'x1', 'y1']]])
    tree = spatial.KDTree(treecoords)
    dist, ind = tree.query(querycoords, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(links)
    found_i = np.arange(len(new_links)*2)[found_b] % len(new_links)
    return pd.DataFrame(dict(D=dist[found_b],
                             i=links.index[ind[found_b] % len(links)]),
                        index=new_links.index[found_i]).sort_values(by='D')\
                        [lambda ds: ~ds.index.duplicated(keep='first')]\
                         .sort_index()['i']


def _load_buses_from_eg():
    buses = (pd.read_csv(snakemake.input.eg_buses, quotechar="'",
                         true_values=['t'], false_values=['f'],
                         dtype=dict(bus_id="str"))
            .set_index("bus_id")
            .drop(['station_id'], axis=1)
            .rename(columns=dict(voltage='v_nom')))

    buses['carrier'] = buses.pop('dc').map({True: 'DC', False: 'AC'})
    buses['under_construction'] = buses['under_construction'].fillna(False).astype(bool)

    # remove all buses outside of all countries including exclusive economic zones (offshore)
    europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, 'geometry']
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    buses_in_europe_b = buses[['x', 'y']].apply(lambda p: europe_shape_prepped.contains(Point(p)), axis=1)

    buses_with_v_nom_to_keep_b = buses.v_nom.isin(snakemake.config['electricity']['voltages']) | buses.v_nom.isnull()
    logger.info("Removing buses with voltages {}".format(pd.Index(buses.v_nom.unique()).dropna().difference(snakemake.config['electricity']['voltages'])))

    return pd.DataFrame(buses.loc[buses_in_europe_b & buses_with_v_nom_to_keep_b])


def _load_transformers_from_eg(buses):
    transformers = (pd.read_csv(snakemake.input.eg_transformers, quotechar="'",
                                true_values=['t'], false_values=['f'],
                                dtype=dict(transformer_id='str', bus0='str', bus1='str'))
                    .set_index('transformer_id'))

    transformers = _remove_dangling_branches(transformers, buses)

    return transformers


def _load_converters_from_eg(buses):
    converters = (pd.read_csv(snakemake.input.eg_converters, quotechar="'",
                              true_values=['t'], false_values=['f'],
                              dtype=dict(converter_id='str', bus0='str', bus1='str'))
                  .set_index('converter_id'))

    converters = _remove_dangling_branches(converters, buses)

    converters['carrier'] = 'B2B'

    return converters


def _load_links_from_eg(buses):
    links = (pd.read_csv(snakemake.input.eg_links, quotechar="'", true_values=['t'], false_values=['f'],
                         dtype=dict(link_id='str', bus0='str', bus1='str', under_construction="bool"))
             .set_index('link_id'))

    links['length'] /= 1e3

    # hotfix
    links.loc[links.bus1=='6271', 'bus1'] = '6273'

    links = _remove_dangling_branches(links, buses)

    # Add DC line parameters
    links['carrier'] = 'DC'

    return links


def _add_links_from_tyndp(buses, links):
    links_tyndp = pd.read_csv(snakemake.input.links_tyndp)

    # remove all links from list which lie outside all of the desired countries
    europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, 'geometry']
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    x1y1_in_europe_b = links_tyndp[['x1', 'y1']].apply(lambda p: europe_shape_prepped.contains(Point(p)), axis=1)
    x2y2_in_europe_b = links_tyndp[['x2', 'y2']].apply(lambda p: europe_shape_prepped.contains(Point(p)), axis=1)
    is_within_covered_countries_b = x1y1_in_europe_b & x2y2_in_europe_b

    if not is_within_covered_countries_b.all():
        logger.info("TYNDP links outside of the covered area (skipping): " +
                    ", ".join(links_tyndp.loc[~ is_within_covered_countries_b, "Name"]))

        links_tyndp = links_tyndp.loc[is_within_covered_countries_b]
        if links_tyndp.empty:
            return buses, links

    has_replaces_b = links_tyndp.replaces.notnull()
    oids = dict(Bus=_get_oid(buses), Link=_get_oid(links))
    keep_b = dict(Bus=pd.Series(True, index=buses.index),
                  Link=pd.Series(True, index=links.index))
    for reps in links_tyndp.loc[has_replaces_b, 'replaces']:
        for comps in reps.split(':'):
            oids_to_remove = comps.split('.')
            c = oids_to_remove.pop(0)
            keep_b[c] &= ~oids[c].isin(oids_to_remove)
    buses = buses.loc[keep_b['Bus']]
    links = links.loc[keep_b['Link']]

    links_tyndp["j"] = _find_closest_links(links, links_tyndp, distance_upper_bound=0.20)
    # Corresponds approximately to 20km tolerances

    if links_tyndp["j"].notnull().any():
        logger.info("TYNDP links already in the dataset (skipping): " + ", ".join(links_tyndp.loc[links_tyndp["j"].notnull(), "Name"]))
        links_tyndp = links_tyndp.loc[links_tyndp["j"].isnull()]
        if links_tyndp.empty: return buses, links

    tree = spatial.KDTree(buses[['x', 'y']])
    _, ind0 = tree.query(links_tyndp[["x1", "y1"]])
    ind0_b = ind0 < len(buses)
    links_tyndp.loc[ind0_b, "bus0"] = buses.index[ind0[ind0_b]]

    _, ind1 = tree.query(links_tyndp[["x2", "y2"]])
    ind1_b = ind1 < len(buses)
    links_tyndp.loc[ind1_b, "bus1"] = buses.index[ind1[ind1_b]]

    links_tyndp_located_b = links_tyndp["bus0"].notnull() & links_tyndp["bus1"].notnull()
    if not links_tyndp_located_b.all():
        logger.warning("Did not find connected buses for TYNDP links (skipping): " + ", ".join(links_tyndp.loc[~links_tyndp_located_b, "Name"]))
        links_tyndp = links_tyndp.loc[links_tyndp_located_b]

    logger.info("Adding the following TYNDP links: " + ", ".join(links_tyndp["Name"]))

    links_tyndp = links_tyndp[["bus0", "bus1"]].assign(
        carrier='DC',
        p_nom=links_tyndp["Power (MW)"],
        length=links_tyndp["Length (given) (km)"].fillna(links_tyndp["Length (distance*1.2) (km)"]),
        under_construction=True,
        underground=False,
        geometry=(links_tyndp[["x1", "y1", "x2", "y2"]]
                  .apply(lambda s: str(LineString([[s.x1, s.y1], [s.x2, s.y2]])), axis=1)),
        tags=('"name"=>"' + links_tyndp["Name"] + '", ' +
              '"ref"=>"' + links_tyndp["Ref"] + '", ' +
              '"status"=>"' + links_tyndp["status"] + '"')
    )

    links_tyndp.index = "T" + links_tyndp.index.astype(str)

    return buses, links.append(links_tyndp, sort=True)


def _load_lines_from_eg(buses):
    lines = (pd.read_csv(snakemake.input.eg_lines, quotechar="'", true_values=['t'], false_values=['f'],
                         dtype=dict(line_id='str', bus0='str', bus1='str',
                                    underground="bool", under_construction="bool"))
             .set_index('line_id')
             .rename(columns=dict(voltage='v_nom', circuits='num_parallel')))

    lines['length'] /= 1e3

    lines = _remove_dangling_branches(lines, buses)

    return lines


def _apply_parameter_corrections(n):
    with open(snakemake.input.parameter_corrections) as f:
        corrections = yaml.safe_load(f)

    if corrections is None: return

    for component, attrs in corrections.items():
        df = n.df(component)
        oid = _get_oid(df)
        if attrs is None: continue

        for attr, repls in attrs.items():
            for i, r in repls.items():
                if i == 'oid':
                    r = oid.map(repls["oid"]).dropna()
                elif i == 'index':
                    r = pd.Series(repls["index"])
                else:
                    raise NotImplementedError()
                inds = r.index.intersection(df.index)
                df.loc[inds, attr] = r[inds].astype(df[attr].dtype)


def _set_electrical_parameters_lines(lines):
    v_noms = snakemake.config['electricity']['voltages']
    linetypes = snakemake.config['lines']['types']

    for v_nom in v_noms:
        lines.loc[lines["v_nom"] == v_nom, 'type'] = linetypes[v_nom]

    lines['s_max_pu'] = snakemake.config['lines']['s_max_pu']

    return lines


def _set_lines_s_nom_from_linetypes(n):
    n.lines['s_nom'] = (
        np.sqrt(3) * n.lines['type'].map(n.line_types.i_nom) *
        n.lines['v_nom'] * n.lines.num_parallel
    )


def _set_electrical_parameters_links(links):
    if links.empty: return links

    p_max_pu = snakemake.config['links'].get('p_max_pu', 1.)
    links['p_max_pu'] = p_max_pu
    links['p_min_pu'] = -p_max_pu

    links_p_nom = pd.read_csv(snakemake.input.links_p_nom)

    # filter links that are not in operation anymore
    removed_b = links_p_nom.Remarks.str.contains('Shut down|Replaced', na=False)
    links_p_nom = links_p_nom[~removed_b]

    # find closest link for all links in links_p_nom
    links_p_nom['j'] = _find_closest_links(links, links_p_nom)

    links_p_nom = links_p_nom.groupby(['j'],as_index=False).agg({'Power (MW)': 'sum'})

    p_nom = links_p_nom.dropna(subset=["j"]).set_index("j")["Power (MW)"]

    # Don't update p_nom if it's already set
    p_nom_unset = p_nom.drop(links.index[links.p_nom.notnull()], errors='ignore') if "p_nom" in links else p_nom
    links.loc[p_nom_unset.index, "p_nom"] = p_nom_unset

    return links


def _set_electrical_parameters_converters(converters):
    p_max_pu = snakemake.config['links'].get('p_max_pu', 1.)
    converters['p_max_pu'] = p_max_pu
    converters['p_min_pu'] = -p_max_pu

    converters['p_nom'] = 2000

    # Converters are combined with links
    converters['under_construction'] = False
    converters['underground'] = False

    return converters


def _set_electrical_parameters_transformers(transformers):
    config = snakemake.config['transformers']

    ## Add transformer parameters
    transformers["x"] = config.get('x', 0.1)
    transformers["s_nom"] = config.get('s_nom', 2000)
    transformers['type'] = config.get('type', '')

    return transformers


def _remove_dangling_branches(branches, buses):
    return pd.DataFrame(branches.loc[branches.bus0.isin(buses.index) & branches.bus1.isin(buses.index)])


def _remove_unconnected_components(network):
    _, labels = csgraph.connected_components(network.adjacency_matrix(), directed=False)
    component = pd.Series(labels, index=network.buses.index)

    component_sizes = component.value_counts()
    components_to_remove = component_sizes.iloc[1:]

    logger.info("Removing {} unconnected network components with less than {} buses. In total {} buses."
                .format(len(components_to_remove), components_to_remove.max(), components_to_remove.sum()))

    return network[component == component_sizes.index[0]]


def _set_countries_and_substations(n):

    buses = n.buses

    def buses_in_shape(shape):
        shape = shapely.prepared.prep(shape)
        return pd.Series(
            np.fromiter((shape.contains(Point(x, y))
                        for x, y in buses.loc[:,["x", "y"]].values),
                        dtype=bool, count=len(buses)),
            index=buses.index
        )

    countries = snakemake.config['countries']
    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index('name')['geometry']
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index('name')['geometry']
    substation_b = buses['symbol'].str.contains('substation|converter station', case=False)

    def prefer_voltage(x, which):
        index = x.index
        if len(index) == 1:
            return pd.Series(index, index)
        key = (x.index[0]
               if x['v_nom'].isnull().all()
               else getattr(x['v_nom'], 'idx' + which)())
        return pd.Series(key, index)

    gb = buses.loc[substation_b].groupby(['x', 'y'], as_index=False,
                                         group_keys=False, sort=False)
    bus_map_low = gb.apply(prefer_voltage, 'min')
    lv_b = (bus_map_low == bus_map_low.index).reindex(buses.index, fill_value=False)
    bus_map_high = gb.apply(prefer_voltage, 'max')
    hv_b = (bus_map_high == bus_map_high.index).reindex(buses.index, fill_value=False)

    onshore_b = pd.Series(False, buses.index)
    offshore_b = pd.Series(False, buses.index)

    for country in countries:
        onshore_shape = country_shapes[country]
        onshore_country_b = buses_in_shape(onshore_shape)
        onshore_b |= onshore_country_b

        buses.loc[onshore_country_b, 'country'] = country

        if country not in offshore_shapes.index: continue
        offshore_country_b = buses_in_shape(offshore_shapes[country])
        offshore_b |= offshore_country_b

        buses.loc[offshore_country_b, 'country'] = country

    # Only accept buses as low-voltage substations (where load is attached), if
    # they have at least one connection which is not under_construction
    has_connections_b = pd.Series(False, index=buses.index)
    for b, df in product(('bus0', 'bus1'), (n.lines, n.links)):
        has_connections_b |= ~ df.groupby(b).under_construction.min()

    buses['substation_lv'] = lv_b & onshore_b & (~ buses['under_construction']) & has_connections_b
    buses['substation_off'] = (offshore_b | (hv_b & onshore_b)) & (~ buses['under_construction'])

    c_nan_b = buses.country.isnull()
    if c_nan_b.sum() > 0:
        c_tag = _get_country(buses.loc[c_nan_b])
        c_tag.loc[~c_tag.isin(countries)] = np.nan
        n.buses.loc[c_nan_b, 'country'] = c_tag

        c_tag_nan_b = n.buses.country.isnull()

        # Nearest country in path length defines country of still homeless buses
        # Work-around until commit 705119 lands in pypsa release
        n.transformers['length'] = 0.
        graph = n.graph(weight='length')
        n.transformers.drop('length', axis=1, inplace=True)

        for b in n.buses.index[c_tag_nan_b]:
            df = (pd.DataFrame(dict(pathlength=nx.single_source_dijkstra_path_length(graph, b, cutoff=200)))
                  .join(n.buses.country).dropna())
            assert not df.empty, "No buses with defined country within 200km of bus `{}`".format(b)
            n.buses.at[b, 'country'] = df.loc[df.pathlength.idxmin(), 'country']

        logger.warning("{} buses are not in any country or offshore shape,"
                       " {} have been assigned from the tag of the entsoe map,"
                       " the rest from the next bus in terms of pathlength."
                       .format(c_nan_b.sum(), c_nan_b.sum() - c_tag_nan_b.sum()))

    return buses


def _replace_b2b_converter_at_country_border_by_link(n):
    # Affects only the B2B converter in Lithuania at the Polish border at the moment
    buscntry = n.buses.country
    linkcntry = n.links.bus0.map(buscntry)
    converters_i = n.links.index[(n.links.carrier == 'B2B') & (linkcntry == n.links.bus1.map(buscntry))]

    def findforeignbus(G, i):
        cntry = linkcntry.at[i]
        for busattr in ('bus0', 'bus1'):
            b0 = n.links.at[i, busattr]
            for b1 in G[b0]:
                if buscntry[b1] != cntry:
                    return busattr, b0, b1
        return None, None, None

    for i in converters_i:
        G = n.graph()
        busattr, b0, b1 = findforeignbus(G, i)
        if busattr is not None:
            comp, line = next(iter(G[b0][b1]))
            if comp != "Line":
                logger.warning("Unable to replace B2B `{}` expected a Line, but found a {}"
                            .format(i, comp))
                continue

            n.links.at[i, busattr] = b1
            n.links.at[i, 'p_nom'] = min(n.links.at[i, 'p_nom'], n.lines.at[line, 's_nom'])
            n.links.at[i, 'carrier'] = 'DC'
            n.links.at[i, 'underwater_fraction'] = 0.
            n.links.at[i, 'length'] = n.lines.at[line, 'length']

            n.remove("Line", line)
            n.remove("Bus", b0)

            logger.info("Replacing B2B converter `{}` together with bus `{}` and line `{}` by an HVDC tie-line {}-{}"
                        .format(i, b0, line, linkcntry.at[i], buscntry.at[b1]))


def _set_links_underwater_fraction(n):
    if n.links.empty: return

    if not hasattr(n.links, 'geometry'):
        n.links['underwater_fraction'] = 0.
    else:
        offshore_shape = gpd.read_file(snakemake.input.offshore_shapes).unary_union
        links = gpd.GeoSeries(n.links.geometry.dropna().map(shapely.wkt.loads))
        n.links['underwater_fraction'] = links.intersection(offshore_shape).length / links.length


def _adjust_capacities_of_under_construction_branches(n):
    lines_mode = snakemake.config['lines'].get('under_construction', 'undef')
    if lines_mode == 'zero':
        n.lines.loc[n.lines.under_construction, 'num_parallel'] = 0.
        n.lines.loc[n.lines.under_construction, 's_nom'] = 0.
    elif lines_mode == 'remove':
        n.mremove("Line", n.lines.index[n.lines.under_construction])
    elif lines_mode != 'keep':
        logger.warning("Unrecognized configuration for `lines: under_construction` = `{}`. Keeping under construction lines.")

    links_mode = snakemake.config['links'].get('under_construction', 'undef')
    if links_mode == 'zero':
        n.links.loc[n.links.under_construction, "p_nom"] = 0.
    elif links_mode == 'remove':
        n.mremove("Link", n.links.index[n.links.under_construction])
    elif links_mode != 'keep':
        logger.warning("Unrecognized configuration for `links: under_construction` = `{}`. Keeping under construction links.")

    if lines_mode == 'remove' or links_mode == 'remove':
        # We might need to remove further unconnected components
        n = _remove_unconnected_components(n)

    return n


def base_network():
    buses = _load_buses_from_eg()

    links = _load_links_from_eg(buses)
    if snakemake.config['links'].get('include_tyndp'):
        buses, links = _add_links_from_tyndp(buses, links)

    converters = _load_converters_from_eg(buses)

    lines = _load_lines_from_eg(buses)
    transformers = _load_transformers_from_eg(buses)

    lines = _set_electrical_parameters_lines(lines)
    transformers = _set_electrical_parameters_transformers(transformers)
    links = _set_electrical_parameters_links(links)
    converters = _set_electrical_parameters_converters(converters)

    n = pypsa.Network()
    n.name = 'PyPSA-Eur'

    n.set_snapshots(pd.date_range(freq='h', **snakemake.config['snapshots']))
    n.snapshot_weightings[:] *= 8760. / n.snapshot_weightings.sum()

    n.import_components_from_dataframe(buses, "Bus")
    n.import_components_from_dataframe(lines, "Line")
    n.import_components_from_dataframe(transformers, "Transformer")
    n.import_components_from_dataframe(links, "Link")
    n.import_components_from_dataframe(converters, "Link")

    _set_lines_s_nom_from_linetypes(n)

    _apply_parameter_corrections(n)

    n = _remove_unconnected_components(n)

    _set_countries_and_substations(n)

    _set_links_underwater_fraction(n)

    _replace_b2b_converter_at_country_border_by_link(n)

    n = _adjust_capacities_of_under_construction_branches(n)

    return n

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('base_network')
    configure_logging(snakemake)

    n = base_network()

    n.export_to_netcdf(snakemake.output[0])

## build_bus_regions
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates Voronoi shapes for each bus representing both onshore and offshore regions.

Relevant Settings
-----------------

.. code:: yaml

    countries:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

Inputs
------

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``resources/regions_onshore.geojson``:

    .. image:: ../img/regions_onshore.png
        :scale: 33 %

- ``resources/regions_offshore.geojson``:

    .. image:: ../img/regions_offshore.png
        :scale: 33 %

Description
-----------

"""

import logging
from _helpers import configure_logging

import pypsa
import os
import pandas as pd
import geopandas as gpd

from vresutils.graph import voronoi_partition_pts

logger = logging.getLogger(__name__)


def save_to_geojson(s, fn):
    if os.path.exists(fn):
        os.unlink(fn)
    schema = {**gpd.io.file.infer_schema(s), 'geometry': 'Unknown'}
    s.to_file(fn, driver='GeoJSON', schema=schema)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_bus_regions')
    configure_logging(snakemake)

    countries = snakemake.config['countries']

    n = pypsa.Network(snakemake.input.base_network)

    country_shapes = gpd.read_file(snakemake.input.country_shapes).set_index('name')['geometry']
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index('name')['geometry']

    onshore_regions = []
    offshore_regions = []

    for country in countries:
        c_b = n.buses.country == country

        onshore_shape = country_shapes[country]
        onshore_locs = n.buses.loc[c_b & n.buses.substation_lv, ["x", "y"]]
        onshore_regions.append(gpd.GeoDataFrame({
                'name': onshore_locs.index,
                'x': onshore_locs['x'],
                'y': onshore_locs['y'],
                'geometry': voronoi_partition_pts(onshore_locs.values, onshore_shape),
                'country': country
            }))

        if country not in offshore_shapes.index: continue
        offshore_shape = offshore_shapes[country]
        offshore_locs = n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]]
        offshore_regions_c = gpd.GeoDataFrame({
                'name': offshore_locs.index,
                'x': offshore_locs['x'],
                'y': offshore_locs['y'],
                'geometry': voronoi_partition_pts(offshore_locs.values, offshore_shape),
                'country': country
            })
        offshore_regions_c = offshore_regions_c.loc[offshore_regions_c.area > 1e-2]
        offshore_regions.append(offshore_regions_c)

    save_to_geojson(pd.concat(onshore_regions, ignore_index=True), snakemake.output.regions_onshore)

    save_to_geojson(pd.concat(offshore_regions, ignore_index=True), snakemake.output.regions_offshore)


##build_cutout
# SPDX-FileCopyrightText: : 2017-2021 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Create cutouts with `atlite <https://atlite.readthedocs.io/en/latest/>`_.

For this rule to work you must have

- installed the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_ ``cdsapi`` package  (`install with `pip``) and
- registered and setup your CDS API key as described `on their website <https://cds.climate.copernicus.eu/api-how-to>`_.

.. seealso::
    For details on the weather data read the `atlite documentation <https://atlite.readthedocs.io/en/latest/>`_.
    If you need help specifically for creating cutouts `the corresponding section in the atlite documentation <https://atlite.readthedocs.io/en/latest/examples/create_cutout.html>`_ should be helpful.

Relevant Settings
-----------------

.. code:: yaml

    atlite:
        nprocesses:
        cutouts:
            {cutout}:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`atlite_cf`

Inputs
------

*None*

Outputs
-------

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_
  reanalysis weather dataset or `SARAH-2 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`_
  satellite-based historic weather data with the following structure:

**ERA5 cutout:**

    ===================  ==========  ==========  =========================================================
    Field                Dimensions  Unit        Description
    ===================  ==========  ==========  =========================================================
    pressure             time, y, x  Pa          Surface pressure
    -------------------  ----------  ----------  ---------------------------------------------------------
    temperature          time, y, x  K           Air temperature 2 meters above the surface.
    -------------------  ----------  ----------  ---------------------------------------------------------
    soil temperature     time, y, x  K           Soil temperature between 1 meters and 3 meters
                                                 depth (layer 4).
    -------------------  ----------  ----------  ---------------------------------------------------------
    influx_toa           time, y, x  Wm**-2      Top of Earth's atmosphere TOA incident solar radiation
    -------------------  ----------  ----------  ---------------------------------------------------------
    influx_direct        time, y, x  Wm**-2      Total sky direct solar radiation at surface
    -------------------  ----------  ----------  ---------------------------------------------------------
    runoff               time, y, x  m           `Runoff <https://en.wikipedia.org/wiki/Surface_runoff>`_
                                                 (volume per area)
    -------------------  ----------  ----------  ---------------------------------------------------------
    roughness            y, x        m           Forecast surface roughness
                                                 (`roughness length <https://en.wikipedia.org/wiki/Roughness_length>`_)
    -------------------  ----------  ----------  ---------------------------------------------------------
    height               y, x        m           Surface elevation above sea level
    -------------------  ----------  ----------  ---------------------------------------------------------
    albedo               time, y, x  --          `Albedo <https://en.wikipedia.org/wiki/Albedo>`_
                                                 measure of diffuse reflection of solar radiation.
                                                 Calculated from relation between surface solar radiation
                                                 downwards (Jm**-2) and surface net solar radiation
                                                 (Jm**-2). Takes values between 0 and 1.
    -------------------  ----------  ----------  ---------------------------------------------------------
    influx_diffuse       time, y, x  Wm**-2      Diffuse solar radiation at surface.
                                                 Surface solar radiation downwards minus
                                                 direct solar radiation.
    -------------------  ----------  ----------  ---------------------------------------------------------
    wnd100m              time, y, x  ms**-1      Wind speeds at 100 meters (regardless of direction)
    ===================  ==========  ==========  =========================================================

    .. image:: ../img/era5.png
        :scale: 40 %

A **SARAH-2 cutout** can be used to amend the fields ``temperature``, ``influx_toa``, ``influx_direct``, ``albedo``,
``influx_diffuse`` of ERA5 using satellite-based radiation observations.

    .. image:: ../img/sarah.png
        :scale: 40 %

Description
-----------

"""

import logging
import atlite
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_cutout', cutout='europe-2013-era5')
    configure_logging(snakemake)

    cutout_params = snakemake.config['atlite']['cutouts'][snakemake.wildcards.cutout]

    snapshots = pd.date_range(freq='h', **snakemake.config['snapshots'])
    time = [snapshots[0], snapshots[-1]]
    cutout_params['time'] = slice(*cutout_params.get('time', time))

    if {'x', 'y', 'bounds'}.isdisjoint(cutout_params):
        # Determine the bounds from bus regions with a buffer of two grid cells
        onshore = gpd.read_file(snakemake.input.regions_onshore)
        offshore = gpd.read_file(snakemake.input.regions_offshore)
        regions =  onshore.append(offshore)
        d = max(cutout_params.get('dx', 0.25), cutout_params.get('dy', 0.25))*2
        cutout_params['bounds'] = regions.total_bounds + [-d, -d, d, d]
    elif {'x', 'y'}.issubset(cutout_params):
        cutout_params['x'] = slice(*cutout_params['x'])
        cutout_params['y'] = slice(*cutout_params['y'])


    logging.info(f"Preparing cutout with parameters {cutout_params}.")
    features = cutout_params.pop('features', None)
    cutout = atlite.Cutout(snakemake.output[0], **cutout_params)
    cutout.prepare(features=features)

#build_hydro_profile
#!/usr/bin/env python

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Build hydroelectric inflow time-series for each country.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    renewable:
        hydro:
            cutout:
            clip_min_inflow:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/bundle/EIA_hydro_generation_2000_2014.csv``: Hydroelectricity net generation per country and year (`EIA <https://www.eia.gov/beta/international/data/browser/#/?pa=000000000000000000000000000000g&c=1028i008006gg6168g80a4k000e0ag00gg0004g800ho00g8&ct=0&ug=8&tl_id=2-A&vs=INTL.33-12-ALB-BKWH.A&cy=2014&vo=0&v=H&start=2000&end=2016>`_)

    .. image:: ../img/hydrogeneration.png
        :scale: 33 %

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``"cutouts/" + config["renewable"]['hydro']['cutout']``: confer :ref:`cutout`

Outputs
-------

- ``resources/profile_hydro.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               countries, time   Inflow to the state of charge (in MW),
                                           e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================

    .. image:: ../img/inflow-ts.png
        :scale: 33 %

    .. image:: ../img/inflow-box.png
        :scale: 33 %

Description
-----------

.. seealso::
    :mod:`build_renewable_profiles`
"""

import logging
from _helpers import configure_logging

import atlite
import geopandas as gpd
from vresutils import hydro as vhydro

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_hydro_profile')
    configure_logging(snakemake)

    config = snakemake.config['renewable']['hydro']
    cutout = atlite.Cutout(snakemake.input.cutout)

    countries = snakemake.config['countries']
    country_shapes = (gpd.read_file(snakemake.input.country_shapes)
                      .set_index('name')['geometry'].reindex(countries))
    country_shapes.index.name = 'countries'

    eia_stats = vhydro.get_eia_annual_hydro_generation(
        snakemake.input.eia_hydro_generation).reindex(columns=countries)
    inflow = cutout.runoff(shapes=country_shapes,
                           smooth=True,
                           lower_threshold_quantile=True,
                           normalize_using_yearly=eia_stats)

    if 'clip_min_inflow' in config:
        inflow = inflow.where(inflow > config['clip_min_inflow'], 0)

    inflow.to_netcdf(snakemake.output[0])

#build_load_data
# SPDX-FileCopyrightText: : 2020 @JanFrederickUnnewehr, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""

This rule downloads the load data from `Open Power System Data Time series <https://data.open-power-system-data.org/time_series/>`_. For all countries in the network, the per country load timeseries with suffix ``_load_actual_entsoe_transparency`` are extracted from the dataset. After filling small gaps linearly and large gaps by copying time-slice of a given period, the load data is exported to a ``.csv`` file.

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    load:
        url:
        interpolate_limit:
        time_shift_for_large_gaps:
        manual_adjustments:


.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`load_cf`

Inputs
------


Outputs
-------

- ``resource/time_series_60min_singleindex_filtered.csv``:


"""

import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging

import pandas as pd
import numpy as np
import dateutil
from pandas import Timedelta as Delta


def load_timeseries(fn, years, countries, powerstatistics=True):
    """
    Read load data from OPSD time-series package version 2020-10-06.

    Parameters
    ----------
    years : None or slice()
        Years for which to read load data (defaults to
        slice("2018","2019"))
    fn : str
        File name or url location (file format .csv)
    countries : listlike
        Countries for which to read load data.
    powerstatistics: bool
        Whether the electricity consumption data of the ENTSOE power
        statistics (if true) or of the ENTSOE transparency map (if false)
        should be parsed.

    Returns
    -------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries
    """
    logger.info(f"Retrieving load data from '{fn}'.")

    pattern = 'power_statistics' if powerstatistics else '_transparency'
    pattern = f'_load_actual_entsoe_{pattern}'
    rename = lambda s: s[:-len(pattern)]
    date_parser = lambda x: dateutil.parser.parse(x, ignoretz=True)
    return (pd.read_csv(fn, index_col=0, parse_dates=[0], date_parser=date_parser)
            .filter(like=pattern)
            .rename(columns=rename)
            .dropna(how="all", axis=0)
            .rename(columns={'GB_UKM' : 'GB'})
            .filter(items=countries)
            .loc[years])


def consecutive_nans(ds):
    return (ds.isnull().astype(int)
            .groupby(ds.notnull().astype(int).cumsum()[ds.isnull()])
            .transform('sum').fillna(0))


def fill_large_gaps(ds, shift):
    """
    Fill up large gaps with load data from the previous week.

    This function fills gaps ragning from 3 to 168 hours (one week).
    """
    shift = Delta(shift)
    nhours = shift / np.timedelta64(1, 'h')
    if (consecutive_nans(ds) > nhours).any():
        logger.warning('There exist gaps larger then the time shift used for '
                       'copying time slices.')
    time_shift = pd.Series(ds.values, ds.index + shift)
    return ds.where(ds.notnull(), time_shift.reindex_like(ds))


def nan_statistics(df):
    def max_consecutive_nans(ds):
        return (ds.isnull().astype(int)
                  .groupby(ds.notnull().astype(int).cumsum())
                  .sum().max())
    consecutive = df.apply(max_consecutive_nans)
    total = df.isnull().sum()
    max_total_per_month = df.isnull().resample('m').sum().max()
    return pd.concat([total, consecutive, max_total_per_month],
                 keys=['total', 'consecutive', 'max_total_per_month'], axis=1)


def copy_timeslice(load, cntry, start, stop, delta):
    start = pd.Timestamp(start)
    stop = pd.Timestamp(stop)
    if start-delta in load.index and stop in load.index and cntry in load:
        load.loc[start:stop, cntry] = load.loc[start-delta:stop-delta, cntry].values


def manual_adjustment(load, powerstatistics):
    """
    Adjust gaps manual for load data from OPSD time-series package.

    1. For the ENTSOE power statistics load data (if powerstatistics is True)

    Kosovo (KV) and Albania (AL) do not exist in the data set. Kosovo gets the
    same load curve as Serbia and Albania the same as Macdedonia, both scaled
    by the corresponding ratio of total energy consumptions reported by
    IEA Data browser [0] for the year 2013.

    2. For the ENTSOE transparency load data (if powerstatistics is False)

    Albania (AL) and Macedonia (MK) do not exist in the data set. Both get the
    same load curve as Montenegro,  scaled by the corresponding ratio of total energy
    consumptions reported by  IEA Data browser [0] for the year 2016.

    [0] https://www.iea.org/data-and-statistics?country=WORLD&fuel=Electricity%20and%20heat&indicator=TotElecCons


    Parameters
    ----------
    load : pd.DataFrame
        Load time-series with UTC timestamps x ISO-2 countries
    powerstatistics: bool
        Whether argument load comprises the electricity consumption data of
        the ENTSOE power statistics or of the ENTSOE transparency map

    Returns
    -------
    load : pd.DataFrame
        Manual adjusted and interpolated load time-series with UTC
        timestamps x ISO-2 countries
    """

    if powerstatistics:
        if 'MK' in load.columns:
            if 'AL' not in load.columns or load.AL.isnull().values.all():
                load['AL'] = load['MK'] * (4.1 / 7.4)
        if 'RS' in load.columns:
            if 'KV' not in load.columns or load.KV.isnull().values.all():
                load['KV'] = load['RS'] * (4.8 / 27.)

        copy_timeslice(load, 'GR', '2015-08-11 21:00', '2015-08-15 20:00', Delta(weeks=1))
        copy_timeslice(load, 'AT', '2018-12-31 22:00', '2019-01-01 22:00', Delta(days=2))
        copy_timeslice(load, 'CH', '2010-01-19 07:00', '2010-01-19 22:00', Delta(days=1))
        copy_timeslice(load, 'CH', '2010-03-28 00:00', '2010-03-28 21:00', Delta(days=1))
        # is a WE, so take WE before
        copy_timeslice(load, 'CH', '2010-10-08 13:00', '2010-10-10 21:00', Delta(weeks=1))
        copy_timeslice(load, 'CH', '2010-11-04 04:00', '2010-11-04 22:00', Delta(days=1))
        copy_timeslice(load, 'NO', '2010-12-09 11:00', '2010-12-09 18:00', Delta(days=1))
        # whole january missing
        copy_timeslice(load, 'GB', '2009-12-31 23:00', '2010-01-31 23:00', Delta(days=-364))

    else:
        if 'ME' in load:
            if 'AL' not in load and 'AL' in countries:
                load['AL'] = load.ME * (5.7/2.9)
            if 'MK' not in load and 'MK' in countries:
                load['MK'] = load.ME * (6.7/2.9)
        copy_timeslice(load, 'BG', '2018-10-27 21:00', '2018-10-28 22:00', Delta(weeks=1))

    return load


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_load_data')

    configure_logging(snakemake)

    config = snakemake.config
    powerstatistics = config['load']['power_statistics']
    url = config['load']['url']
    interpolate_limit = config['load']['interpolate_limit']
    countries = config['countries']
    snapshots = pd.date_range(freq='h', **config['snapshots'])
    years = slice(snapshots[0], snapshots[-1])
    time_shift = config['load']['time_shift_for_large_gaps']

    load = load_timeseries(url, years, countries, powerstatistics)

    if config['load']['manual_adjustments']:
        load = manual_adjustment(load, powerstatistics)

    logger.info(f"Linearly interpolate gaps of size {interpolate_limit} and less.")
    load = load.interpolate(method='linear', limit=interpolate_limit)

    logger.info("Filling larger gaps by copying time-slices of period "
                f"'{time_shift}'.")
    load = load.apply(fill_large_gaps, shift=time_shift)

    assert not load.isna().any().any(), (
        'Load data contains nans. Adjust the parameters '
        '`time_shift_for_large_gaps` or modify the `manual_adjustment` function '
        'for implementing the needed load data modifications.')

    load.to_csv(snakemake.output[0])

#build_natura_raster
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Rasters the vector data of the `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas onto all cutout regions.

Relevant Settings
-----------------

.. code:: yaml

    renewable:
        {technology}:
            cutout:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`renewable_cf`

Inputs
------

- ``data/bundle/natura/Natura2000_end2015.shp``: `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas.

    .. image:: ../img/natura.png
        :scale: 33 %

Outputs
-------

- ``resources/natura.tiff``: Rasterized version of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times.

    .. image:: ../img/natura.png
        :scale: 33 %

Description
-----------

"""

import logging
from _helpers import configure_logging

import atlite
import geopandas as gpd
import rasterio as rio
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)


def determine_cutout_xXyY(cutout_name):
    cutout = atlite.Cutout(cutout_name)
    assert cutout.crs.to_epsg() == 4326
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    return [x - dx/2., X + dx/2., y - dy/2., Y + dy/2.]


def get_transform_and_shape(bounds, res):
    left, bottom = [(b // res)* res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    shape = int((top - bottom) // res), int((right - left) / res)
    transform = rio.Affine(res, 0, left, 0, -res, top)
    return transform, shape


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_natura_raster')
    configure_logging(snakemake)


    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))
    bounds = transform_bounds(4326, 3035, min(xs), min(ys), max(Xs), max(Ys))
    transform, out_shape = get_transform_and_shape(bounds, res=100)

    # adjusted boundaries
    shapes = gpd.read_file(snakemake.input.natura).to_crs(3035)
    raster = ~geometry_mask(shapes.geometry, out_shape[::-1], transform)
    raster = raster.astype(rio.uint8)

    with rio.open(snakemake.output[0], 'w', driver='GTiff', dtype=rio.uint8,
                  count=1, transform=transform, crs=3035, compress='lzw',
                  width=raster.shape[1], height=raster.shape[0]) as dst:
        dst.write(raster, indexes=1)

#build_powerplants
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Retrieves conventional powerplant capacities and locations from `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_, assigns these to buses and creates a ``.csv`` file. It is possible to amend the powerplant database with custom entries provided in ``data/custom_powerplants.csv``.

Relevant Settings
-----------------

.. code:: yaml

    electricity:
      powerplants_filter:
      custom_powerplants:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity`

Inputs
------

- ``networks/base.nc``: confer :ref:`base`.
- ``data/custom_powerplants.csv``: custom powerplants in the same format as `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`_ provides

Outputs
-------

- ``resource/powerplants.csv``: A list of conventional power plants (i.e. neither wind nor solar) with fields for name, fuel type, technology, country, capacity in MW, duration, commissioning year, retrofit year, latitude, longitude, and dam information as documented in the `powerplantmatching README <https://github.com/FRESNA/powerplantmatching/blob/master/README.md>`_; additionally it includes information on the closest substation/bus in ``networks/base.nc``.

    .. image:: ../img/powerplantmatching.png
        :scale: 30 %

    **Source:** `powerplantmatching on GitHub <https://github.com/FRESNA/powerplantmatching>`_

Description
-----------

The configuration options ``electricity: powerplants_filter`` and ``electricity: custom_powerplants`` can be used to control whether data should be retrieved from the original powerplants database or from custom amendmends. These specify `pandas.query <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html>`_ commands.

1. Adding all powerplants from custom:

    .. code:: yaml

        powerplants_filter: false
        custom_powerplants: true

2. Replacing powerplants in e.g. Germany by custom data:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: true

    or

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: Country in ['Germany']


3. Adding additional built year constraints:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany'] and YearCommissioned <= 2015
        custom_powerplants: YearCommissioned <= 2015

"""

import logging
from _helpers import configure_logging

import pypsa
import powerplantmatching as pm
import pandas as pd
import numpy as np

from scipy.spatial import cKDTree as KDTree

logger = logging.getLogger(__name__)


def add_custom_powerplants(ppl):
    custom_ppl_query = snakemake.config['electricity']['custom_powerplants']
    if not custom_ppl_query:
        return ppl
    add_ppls = pd.read_csv(snakemake.input.custom_powerplants, index_col=0,
                           dtype={'bus': 'str'})
    if isinstance(custom_ppl_query, str):
        add_ppls.query(custom_ppl_query, inplace=True)
    return ppl.append(add_ppls, sort=False, ignore_index=True, verify_integrity=True)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_powerplants')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    countries = n.buses.country.unique()

    ppl = (pm.powerplants(from_url=True)
           .powerplant.fill_missing_decommyears()
           .powerplant.convert_country_to_alpha2()
           .query('Fueltype not in ["Solar", "Wind"] and Country in @countries')
           .replace({'Technology': {'Steam Turbine': 'OCGT'}})
            .assign(Fueltype=lambda df: (
                    df.Fueltype
                      .where(df.Fueltype != 'Natural Gas',
                             df.Technology.replace('Steam Turbine',
                                                   'OCGT').fillna('OCGT')))))

    ppl_query = snakemake.config['electricity']['powerplants_filter']
    if isinstance(ppl_query, str):
        ppl.query(ppl_query, inplace=True)

    ppl = add_custom_powerplants(ppl) # add carriers from own powerplant files

    cntries_without_ppl = [c for c in countries if c not in ppl.Country.unique()]

    for c in countries:
        substation_i = n.buses.query('substation_lv and country == @c').index
        kdtree = KDTree(n.buses.loc[substation_i, ['x','y']].values)
        ppl_i = ppl.query('Country == @c').index

        tree_i = kdtree.query(ppl.loc[ppl_i, ['lon','lat']].values)[1]
        ppl.loc[ppl_i, 'bus'] = substation_i.append(pd.Index([np.nan]))[tree_i]

    if cntries_without_ppl:
        logging.warning(f"No powerplants known in: {', '.join(cntries_without_ppl)}")

    bus_null_b = ppl["bus"].isnull()
    if bus_null_b.any():
        logging.warning(f"Couldn't find close bus for {bus_null_b.sum()} powerplants")

    ppl.to_csv(snakemake.output[0])

## build_renewable_profiles
#!/usr/bin/env python

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""Calculates for each network node the
(i) installable capacity (based on land-use), (ii) the available generation time
series (based on weather data), and (iii) the average distance from the node for
onshore wind, AC-connected offshore wind, DC-connected offshore wind and solar
PV generators. In addition for offshore wind it calculates the fraction of the
grid connection which is under water.

.. note:: Hydroelectric profiles are built in script :mod:`build_hydro_profiles`.

Relevant settings
-----------------

.. code:: yaml

    snapshots:

    atlite:
        nprocesses:

    renewable:
        {technology}:
            cutout:
            corine:
            grid_codes:
            distance:
            natura:
            max_depth:
            max_shore_distance:
            min_shore_distance:
            capacity_per_sqkm:
            correction_factor:
            potential:
            min_p_max_pu:
            clip_p_max_pu:
            resource:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`snapshots_cf`, :ref:`atlite_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/bundle/corine/g250_clc06_V18_5.tif``: `CORINE Land Cover (CLC) <https://land.copernicus.eu/pan-european/corine-land-cover>`_ inventory on `44 classes <https://wiki.openstreetmap.org/wiki/Corine_Land_Cover#Tagging>`_ of land use (e.g. forests, arable land, industrial, urban areas).

    .. image:: ../img/corine.png
        :scale: 33 %

- ``data/bundle/GEBCO_2014_2D.nc``: A `bathymetric <https://en.wikipedia.org/wiki/Bathymetry>`_ data set with a global terrain model for ocean and land at 15 arc-second intervals by the `General Bathymetric Chart of the Oceans (GEBCO) <https://www.gebco.net/data_and_products/gridded_bathymetry_data/>`_.

    .. image:: ../img/gebco_2019_grid_image.jpg
        :scale: 50 %

    **Source:** `GEBCO <https://www.gebco.net/data_and_products/images/gebco_2019_grid_image.jpg>`_

- ``resources/natura.tiff``: confer :ref:`natura`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/regions_onshore.geojson``: (if not offshore wind), confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: (if offshore wind), :ref:`busregions`
- ``"cutouts/" + config["renewable"][{technology}]['cutout']``: :ref:`cutout`
- ``networks/base.nc``: :ref:`base`

Outputs
-------

- ``resources/profile_{technology}.nc`` with the following structure

    ===================  ==========  =========================================================
    Field                Dimensions  Description
    ===================  ==========  =========================================================
    profile              bus, time   the per unit hourly availability factors for each node
    -------------------  ----------  ---------------------------------------------------------
    weight               bus         sum of the layout weighting for each node
    -------------------  ----------  ---------------------------------------------------------
    p_nom_max            bus         maximal installable capacity at the node (in MW)
    -------------------  ----------  ---------------------------------------------------------
    potential            y, x        layout of generator units at cutout grid cells inside the
                                     Voronoi cell (maximal installable capacity at each grid
                                     cell multiplied by capacity factor)
    -------------------  ----------  ---------------------------------------------------------
    average_distance     bus         average distance of units in the Voronoi cell to the
                                     grid node (in km)
    -------------------  ----------  ---------------------------------------------------------
    underwater_fraction  bus         fraction of the average connection distance which is
                                     under water (only for offshore)
    ===================  ==========  =========================================================

    - **profile**

    .. image:: ../img/profile_ts.png
        :scale: 33 %
        :align: center

    - **p_nom_max**

    .. image:: ../img/p_nom_max_hist.png
        :scale: 33 %
        :align: center

    - **potential**

    .. image:: ../img/potential_heatmap.png
        :scale: 33 %
        :align: center

    - **average_distance**

    .. image:: ../img/distance_hist.png
        :scale: 33 %
        :align: center

    - **underwater_fraction**

    .. image:: ../img/underwater_hist.png
        :scale: 33 %
        :align: center

Description
-----------

This script functions at two main spatial resolutions: the resolution of the
network nodes and their `Voronoi cells
<https://en.wikipedia.org/wiki/Voronoi_diagram>`_, and the resolution of the
cutout grid cells for the weather data. Typically the weather data grid is
finer than the network nodes, so we have to work out the distribution of
generators across the grid cells within each Voronoi cell. This is done by
taking account of a combination of the available land at each grid cell and the
capacity factor there.

First the script computes how much of the technology can be installed at each
cutout grid cell and each node using the `GLAES
<https://github.com/FZJ-IEK3-VSA/glaes>`_ library. This uses the CORINE land use data,
Natura2000 nature reserves and GEBCO bathymetry data.

.. image:: ../img/eligibility.png
    :scale: 50 %
    :align: center

To compute the layout of generators in each node's Voronoi cell, the
installable potential in each grid cell is multiplied with the capacity factor
at each grid cell. This is done since we assume more generators are installed
at cells with a higher capacity factor.

.. image:: ../img/offwinddc-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/offwindac-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/onwind-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/solar-gridcell.png
    :scale: 50 %
    :align: center

This layout is then used to compute the generation availability time series
from the weather data cutout from ``atlite``.

Two methods are available to compute the maximal installable potential for the
node (`p_nom_max`): ``simple`` and ``conservative``:

- ``simple`` adds up the installable potentials of the individual grid cells.
  If the model comes close to this limit, then the time series may slightly
  overestimate production since it is assumed the geographical distribution is
  proportional to capacity factor.

- ``conservative`` assertains the nodal limit by increasing capacities
  proportional to the layout until the limit of an individual grid cell is
  reached.

"""
import progressbar as pgb
import geopandas as gpd
import xarray as xr
import numpy as np
import functools
import atlite
import logging
from pypsa.geo import haversine
from shapely.geometry import LineString
import time

from _helpers import configure_logging

logger = logging.getLogger(__name__)


if __name__ == '__main__':
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_renewable_profiles', technology='solar')
    configure_logging(snakemake)
    pgb.streams.wrap_stderr()
    paths = snakemake.input
    nprocesses = snakemake.config['atlite'].get('nprocesses')
    noprogress = not snakemake.config['atlite'].get('show_progress', True)
    config = snakemake.config['renewable'][snakemake.wildcards.technology]
    resource = config['resource'] # pv panel config / wind turbine config
    correction_factor = config.get('correction_factor', 1.)
    capacity_per_sqkm = config['capacity_per_sqkm']
    p_nom_max_meth = config.get('potential', 'conservative')

    if isinstance(config.get("corine", {}), list):
        config['corine'] = {'grid_codes': config['corine']}

    if correction_factor != 1.:
        logger.info(f'correction_factor is set as {correction_factor}')


    cutout = atlite.Cutout(paths['cutout'])
    regions = gpd.read_file(paths.regions).set_index('name').rename_axis('bus')
    buses = regions.index

    excluder = atlite.ExclusionContainer(crs=3035, res=100)

    if config['natura']:
        excluder.add_raster(paths.natura, nodata=0, allow_no_overlap=True)

    corine = config.get("corine", {})
    if "grid_codes" in corine:
        codes = corine["grid_codes"]
        excluder.add_raster(paths.corine, codes=codes, invert=True, crs=3035)
    if corine.get("distance", 0.) > 0.:
        codes = corine["distance_grid_codes"]
        buffer = corine["distance"]
        excluder.add_raster(paths.corine, codes=codes, buffer=buffer, crs=3035)

    if "max_depth" in config:
        # lambda not supported for atlite + multiprocessing
        # use named function np.greater with partially frozen argument instead
        # and exclude areas where: -max_depth > grid cell depth
        func = functools.partial(np.greater,-config['max_depth'])
        excluder.add_raster(paths.gebco, codes=func, crs=4236, nodata=-1000)

    if 'min_shore_distance' in config:
        buffer = config['min_shore_distance']
        excluder.add_geometry(paths.country_shapes, buffer=buffer)

    if 'max_shore_distance' in config:
        buffer = config['max_shore_distance']
        excluder.add_geometry(paths.country_shapes, buffer=buffer, invert=True)

    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)
    if noprogress:
        logger.info('Calculate landuse availabilities...')
        start = time.time()
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)
        duration = time.time() - start
        logger.info(f'Completed availability calculation ({duration:2.2f}s)')
    else:
        availability = cutout.availabilitymatrix(regions, excluder, **kwargs)

    area = cutout.grid.to_crs(3035).area / 1e6
    area = xr.DataArray(area.values.reshape(cutout.shape),
                        [cutout.coords['y'], cutout.coords['x']])

    potential = capacity_per_sqkm * availability.sum('bus') * area
    func = getattr(cutout, resource.pop('method'))
    resource['dask_kwargs'] = {'num_workers': nprocesses}
    capacity_factor = correction_factor * func(capacity_factor=True, **resource)
    layout = capacity_factor * area * capacity_per_sqkm
    profile, capacities = func(matrix=availability.stack(spatial=['y','x']),
                                layout=layout, index=buses,
                                per_unit=True, return_capacity=True, **resource)

    logger.info(f"Calculating maximal capacity per bus (method '{p_nom_max_meth}')")
    if p_nom_max_meth == 'simple':
        p_nom_max = capacity_per_sqkm * availability @ area
    elif p_nom_max_meth == 'conservative':
        max_cap_factor = capacity_factor.where(availability!=0).max(['x', 'y'])
        p_nom_max = capacities / max_cap_factor
    else:
        raise AssertionError('Config key `potential` should be one of "simple" '
                        f'(default) or "conservative", not "{p_nom_max_meth}"')



    logger.info('Calculate average distances.')
    layoutmatrix = (layout * availability).stack(spatial=['y','x'])

    coords = cutout.grid[['x', 'y']]
    bus_coords = regions[['x', 'y']]

    average_distance = []
    centre_of_mass = []
    for bus in buses:
        row = layoutmatrix.sel(bus=bus).data
        nz_b = row != 0
        row = row[nz_b]
        co = coords[nz_b]
        distances = haversine(bus_coords.loc[bus],  co)
        average_distance.append((distances * (row / row.sum())).sum())
        centre_of_mass.append(co.values.T @ (row / row.sum()))

    average_distance = xr.DataArray(average_distance, [buses])
    centre_of_mass = xr.DataArray(centre_of_mass, [buses, ('spatial', ['x', 'y'])])


    ds = xr.merge([(correction_factor * profile).rename('profile'),
                    capacities.rename('weight'),
                    p_nom_max.rename('p_nom_max'),
                    potential.rename('potential'),
                    average_distance.rename('average_distance')])


    if snakemake.wildcards.technology.startswith("offwind"):
        logger.info('Calculate underwater fraction of connections.')
        offshore_shape = gpd.read_file(paths['offshore_shapes']).unary_union
        underwater_fraction = []
        for bus in buses:
            p = centre_of_mass.sel(bus=bus).data
            line = LineString([p, regions.loc[bus, ['x', 'y']]])
            frac = line.intersection(offshore_shape).length/line.length
            underwater_fraction.append(frac)

        ds['underwater_fraction'] = xr.DataArray(underwater_fraction, [buses])

    # select only buses with some capacity and minimal capacity factor
    ds = ds.sel(bus=((ds['profile'].mean('time') > config.get('min_p_max_pu', 0.)) &
                      (ds['p_nom_max'] > config.get('min_p_nom_max', 0.))))

    if 'clip_p_max_pu' in config:
        min_p_max_pu = config['clip_p_max_pu']
        ds['profile'] = ds['profile'].where(ds['profile'] >= min_p_max_pu, 0)

    ds.to_netcdf(snakemake.output.profile)

## build_shapes

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates GIS shape files of the countries, exclusive economic zones and `NUTS3 <https://en.wikipedia.org/wiki/Nomenclature_of_Territorial_Units_for_Statistics>`_ areas.

Relevant Settings
-----------------

.. code:: yaml

    countries:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

Inputs
------

- ``data/bundle/naturalearth/ne_10m_admin_0_countries.shp``: World country shapes

    .. image:: ../img/countries.png
        :scale: 33 %

- ``data/bundle/eez/World_EEZ_v8_2014.shp``: World `exclusive economic zones <https://en.wikipedia.org/wiki/Exclusive_economic_zone>`_ (EEZ)

    .. image:: ../img/eez.png
        :scale: 33 %

- ``data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp``: Europe NUTS3 regions

    .. image:: ../img/nuts3.png
        :scale: 33 %

- ``data/bundle/nama_10r_3popgdp.tsv.gz``: Average annual population by NUTS3 region (`eurostat <http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nama_10r_3popgdp&lang=en>`__)
- ``data/bundle/nama_10r_3gdp.tsv.gz``: Gross domestic product (GDP) by NUTS 3 regions (`eurostat <http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nama_10r_3gdp&lang=en>`__)
- ``data/bundle/ch_cantons.csv``: Mapping between Swiss Cantons and NUTS3 regions
- ``data/bundle/je-e-21.03.02.xls``: Population and GDP data per Canton (`BFS - Swiss Federal Statistical Office <https://www.bfs.admin.ch/bfs/en/home/news/whats-new.assetdetail.7786557.html>`_ )

Outputs
-------

- ``resources/country_shapes.geojson``: country shapes out of country selection

    .. image:: ../img/country_shapes.png
        :scale: 33 %

- ``resources/offshore_shapes.geojson``: EEZ shapes out of country selection

    .. image:: ../img/offshore_shapes.png
        :scale: 33 %

- ``resources/europe_shape.geojson``: Shape of Europe including countries and EEZ

    .. image:: ../img/europe_shape.png
        :scale: 33 %

- ``resources/nuts3_shapes.geojson``: NUTS3 shapes out of country selection including population and GDP data.

    .. image:: ../img/nuts3_shapes.png
        :scale: 33 %

Description
-----------

"""

import logging
from _helpers import configure_logging

import os
import numpy as np
from operator import attrgetter
from functools import reduce
from itertools import takewhile

import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import cascaded_union
import pycountry as pyc

logger = logging.getLogger(__name__)


def _get_country(target, **keys):
    assert len(keys) == 1
    try:
        return getattr(pyc.countries.get(**keys), target)
    except (KeyError, AttributeError):
        return np.nan


def _simplify_polys(polys, minarea=0.1, tolerance=0.01, filterremote=True):
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys, key=attrgetter('area'), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area/(2.*np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon([p
                                  for p in takewhile(lambda p: p.area > minarea, polys)
                                  if not filterremote or (mainpoly.distance(p) < mainlength)])
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)


def countries():
    cntries = snakemake.config['countries']
    if 'RS' in cntries: cntries.append('KV')

    df = gpd.read_file(snakemake.input.naturalearth)

    # Names are a hassle in naturalearth, try several fields
    fieldnames = (df[x].where(lambda s: s!='-99') for x in ('ISO_A2', 'WB_A2', 'ADM0_A3'))
    df['name'] = reduce(lambda x,y: x.fillna(y), fieldnames, next(fieldnames)).str[0:2]

    df = df.loc[df.name.isin(cntries) & ((df['scalerank'] == 0) | (df['scalerank'] == 5))]
    s = df.set_index('name')['geometry'].map(_simplify_polys)
    if 'RS' in cntries: s['RS'] = s['RS'].union(s.pop('KV'))

    return s


def eez(country_shapes):
    df = gpd.read_file(snakemake.input.eez)
    df = df.loc[df['ISO_3digit'].isin([_get_country('alpha_3', alpha_2=c) for c in snakemake.config['countries']])]
    df['name'] = df['ISO_3digit'].map(lambda c: _get_country('alpha_2', alpha_3=c))
    s = df.set_index('name').geometry.map(lambda s: _simplify_polys(s, filterremote=False))
    s = gpd.GeoSeries({k:v for k,v in s.iteritems() if v.distance(country_shapes[k]) < 1e-3})
    s.index.name = "name"
    return s


def country_cover(country_shapes, eez_shapes=None):
    shapes = list(country_shapes)
    if eez_shapes is not None:
        shapes += list(eez_shapes)

    europe_shape = cascaded_union(shapes)
    if isinstance(europe_shape, MultiPolygon):
        europe_shape = max(europe_shape, key=attrgetter('area'))
    return Polygon(shell=europe_shape.exterior)


def nuts3(country_shapes):
    df = gpd.read_file(snakemake.input.nuts3)
    df = df.loc[df['STAT_LEVL_'] == 3]
    df['geometry'] = df['geometry'].map(_simplify_polys)
    df = df.rename(columns={'NUTS_ID': 'id'})[['id', 'geometry']].set_index('id')

    pop = pd.read_table(snakemake.input.nuts3pop, na_values=[':'], delimiter=' ?\t', engine='python')
    pop = (pop
           .set_index(pd.MultiIndex.from_tuples(pop.pop('unit,geo\\time').str.split(','))).loc['THS']
           .applymap(lambda x: pd.to_numeric(x, errors='coerce'))
           .fillna(method='bfill', axis=1))['2014']

    gdp = pd.read_table(snakemake.input.nuts3gdp, na_values=[':'], delimiter=' ?\t', engine='python')
    gdp = (gdp
           .set_index(pd.MultiIndex.from_tuples(gdp.pop('unit,geo\\time').str.split(','))).loc['EUR_HAB']
           .applymap(lambda x: pd.to_numeric(x, errors='coerce'))
           .fillna(method='bfill', axis=1))['2014']

    cantons = pd.read_csv(snakemake.input.ch_cantons)
    cantons = cantons.set_index(cantons['HASC'].str[3:])['NUTS']
    cantons = cantons.str.pad(5, side='right', fillchar='0')

    swiss = pd.read_excel(snakemake.input.ch_popgdp, skiprows=3, index_col=0)
    swiss.columns = swiss.columns.to_series().map(cantons)

    pop = pop.append(pd.to_numeric(swiss.loc['Residents in 1000', 'CH040':]))
    gdp = gdp.append(pd.to_numeric(swiss.loc['Gross domestic product per capita in Swiss francs', 'CH040':]))

    df = df.join(pd.DataFrame(dict(pop=pop, gdp=gdp)))

    df['country'] = df.index.to_series().str[:2].replace(dict(UK='GB', EL='GR'))

    excludenuts = pd.Index(('FRA10', 'FRA20', 'FRA30', 'FRA40', 'FRA50',
                            'PT200', 'PT300',
                            'ES707', 'ES703', 'ES704','ES705', 'ES706', 'ES708', 'ES709',
                            'FI2', 'FR9'))
    excludecountry = pd.Index(('MT', 'TR', 'LI', 'IS', 'CY', 'KV'))

    df = df.loc[df.index.difference(excludenuts)]
    df = df.loc[~df.country.isin(excludecountry)]

    manual = gpd.GeoDataFrame(
        [['BA1', 'BA', 3871.],
         ['RS1', 'RS', 7210.],
         ['AL1', 'AL', 2893.]],
        columns=['NUTS_ID', 'country', 'pop']
    ).set_index('NUTS_ID')
    manual['geometry'] = manual['country'].map(country_shapes)
    manual = manual.dropna()

    df = df.append(manual, sort=False)

    df.loc['ME000', 'pop'] = 650.

    return df


def save_to_geojson(df, fn):
    if os.path.exists(fn):
        os.unlink(fn)
    if not isinstance(df, gpd.GeoDataFrame):
        df = gpd.GeoDataFrame(dict(geometry=df))
    df = df.reset_index()
    schema = {**gpd.io.file.infer_schema(df), 'geometry': 'Unknown'}
    df.to_file(fn, driver='GeoJSON', schema=schema)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_shapes')
    configure_logging(snakemake)

    out = snakemake.output

    country_shapes = countries()
    save_to_geojson(country_shapes, out.country_shapes)

    offshore_shapes = eez(country_shapes)
    save_to_geojson(offshore_shapes, out.offshore_shapes)

    europe_shape = country_cover(country_shapes, offshore_shapes)
    save_to_geojson(gpd.GeoSeries(europe_shape), out.europe_shape)

    nuts3_shapes = nuts3(country_shapes)
    save_to_geojson(nuts3_shapes, out.nuts3_shapes)

## cluster_network
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Creates networks clustered to ``{cluster}`` number of zones with aggregated buses, generators and transmission corridors.

Relevant Settings
-----------------

.. code:: yaml

    focus_weights:

    renewable: (keys)
        {technology}:
            potential:

    solving:
        solver:
            name:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`, :ref:`solving_cf`, :ref:`lines_cf`

Inputs
------

- ``resources/regions_onshore_elec_s{simpl}.geojson``: confer :ref:`simplify`
- ``resources/regions_offshore_elec_s{simpl}.geojson``: confer :ref:`simplify`
- ``resources/busmap_elec_s{simpl}.csv``: confer :ref:`simplify`
- ``networks/elec_s{simpl}.nc``: confer :ref:`simplify`
- ``data/custom_busmap_elec_s{simpl}_{clusters}.csv``: optional input

Outputs
-------

- ``resources/regions_onshore_elec_s{simpl}_{clusters}.geojson``:

    .. image:: ../img/regions_onshore_elec_s_X.png
        :scale: 33 %

- ``resources/regions_offshore_elec_s{simpl}_{clusters}.geojson``:

    .. image:: ../img/regions_offshore_elec_s_X.png
        :scale: 33 %

- ``resources/busmap_elec_s{simpl}_{clusters}.csv``: Mapping of buses from ``networks/elec_s{simpl}.nc`` to ``networks/elec_s{simpl}_{clusters}.nc``;
- ``resources/linemap_elec_s{simpl}_{clusters}.csv``: Mapping of lines from ``networks/elec_s{simpl}.nc`` to ``networks/elec_s{simpl}_{clusters}.nc``;
- ``networks/elec_s{simpl}_{clusters}.nc``:

    .. image:: ../img/elec_s_X.png
        :scale: 40  %

Description
-----------

.. note::

    **Why is clustering used both in** ``simplify_network`` **and** ``cluster_network`` **?**

        Consider for example a network ``networks/elec_s100_50.nc`` in which
        ``simplify_network`` clusters the network to 100 buses and in a second
        step ``cluster_network``` reduces it down to 50 buses.

        In preliminary tests, it turns out, that the principal effect of
        changing spatial resolution is actually only partially due to the
        transmission network. It is more important to differentiate between
        wind generators with higher capacity factors from those with lower
        capacity factors, i.e. to have a higher spatial resolution in the
        renewable generation than in the number of buses.

        The two-step clustering allows to study this effect by looking at
        networks like ``networks/elec_s100_50m.nc``. Note the additional
        ``m`` in the ``{cluster}`` wildcard. So in the example network
        there are still up to 100 different wind generators.

        In combination these two features allow you to study the spatial
        resolution of the transmission network separately from the
        spatial resolution of renewable generators.

    **Is it possible to run the model without the** ``simplify_network`` **rule?**

        No, the network clustering methods in the PyPSA module
        `pypsa.networkclustering <https://github.com/PyPSA/PyPSA/blob/master/pypsa/networkclustering.py>`_
        do not work reliably with multiple voltage levels and transformers.

.. tip::
    The rule :mod:`cluster_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`cluster_network`.

Exemplary unsolved network clustered to 512 nodes:

.. image:: ../img/elec_s_512.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 256 nodes:

.. image:: ../img/elec_s_256.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 128 nodes:

.. image:: ../img/elec_s_128.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 37 nodes:

.. image:: ../img/elec_s_37.png
    :scale: 40  %
    :align: center

"""

import logging
from _helpers import configure_logging, update_p_nom_max

import pypsa
import os
import shapely

import pandas as pd
import numpy as np
import geopandas as gpd
import pyomo.environ as po
import matplotlib.pyplot as plt
import seaborn as sns

from functools import reduce

from pypsa.networkclustering import (busmap_by_kmeans, busmap_by_spectral_clustering,
                                     _make_consense, get_clustering_from_busmap)

from add_electricity import load_costs

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(x): return (x/x.sum()).fillna(0.)


def weighting_for_country(n, x):
    conv_carriers = {'OCGT','CCGT','PHS', 'hydro'}
    gen = (n
           .generators.loc[n.generators.carrier.isin(conv_carriers)]
           .groupby('bus').p_nom.sum()
           .reindex(n.buses.index, fill_value=0.) +
           n
           .storage_units.loc[n.storage_units.carrier.isin(conv_carriers)]
           .groupby('bus').p_nom.sum()
           .reindex(n.buses.index, fill_value=0.))
    load = n.loads_t.p_set.mean().groupby(n.loads.bus).sum()

    b_i = x.index
    g = normed(gen.reindex(b_i, fill_value=0))
    l = normed(load.reindex(b_i, fill_value=0))

    w = g + l
    return (w * (100. / w.max())).clip(lower=1.).astype(int)


def distribute_clusters(n, n_clusters, focus_weights=None, solver_name=None):
    """Determine the number of clusters per country"""

    if solver_name is None:
        solver_name = snakemake.config['solving']['solver']['name']

    L = (n.loads_t.p_set.mean()
         .groupby(n.loads.bus).sum()
         .groupby([n.buses.country, n.buses.sub_network]).sum()
         .pipe(normed))

    N = n.buses.groupby(['country', 'sub_network']).size()

    assert n_clusters >= len(N) and n_clusters <= N.sum(), \
        f"Number of clusters must be {len(N)} <= n_clusters <= {N.sum()} for this selection of countries."

    if focus_weights is not None:

        total_focus = sum(list(focus_weights.values()))

        assert total_focus <= 1.0, "The sum of focus weights must be less than or equal to 1."

        for country, weight in focus_weights.items():
            L[country] = weight / len(L[country])

        remainder = [c not in focus_weights.keys() for c in L.index.get_level_values('country')]
        L[remainder] = L.loc[remainder].pipe(normed) * (1 - total_focus)

        logger.warning('Using custom focus weights for determining number of clusters.')

    assert np.isclose(L.sum(), 1.0, rtol=1e-3), f"Country weights L must sum up to 1.0 when distributing clusters. Is {L.sum()}."

    m = po.ConcreteModel()
    def n_bounds(model, *n_id):
        return (1, N[n_id])
    m.n = po.Var(list(L.index), bounds=n_bounds, domain=po.Integers)
    m.tot = po.Constraint(expr=(po.summation(m.n) == n_clusters))
    m.objective = po.Objective(expr=sum((m.n[i] - L.loc[i]*n_clusters)**2 for i in L.index),
                               sense=po.minimize)

    opt = po.SolverFactory(solver_name)
    if not opt.has_capability('quadratic_objective'):
        logger.warning(f'The configured solver `{solver_name}` does not support quadratic objectives. Falling back to `ipopt`.')
        opt = po.SolverFactory('ipopt')

    results = opt.solve(m)
    assert results['Solver'][0]['Status'] == 'ok', f"Solver returned non-optimally: {results}"

    return pd.Series(m.n.get_values(), index=L.index).astype(int)


def busmap_for_n_clusters(n, n_clusters, solver_name, focus_weights=None, algorithm="kmeans", **algorithm_kwds):
    if algorithm == "kmeans":
        algorithm_kwds.setdefault('n_init', 1000)
        algorithm_kwds.setdefault('max_iter', 30000)
        algorithm_kwds.setdefault('tol', 1e-6)

    n.determine_network_topology()

    n_clusters = distribute_clusters(n, n_clusters, focus_weights=focus_weights, solver_name=solver_name)

    def reduce_network(n, buses):
        nr = pypsa.Network()
        nr.import_components_from_dataframe(buses, "Bus")
        nr.import_components_from_dataframe(n.lines.loc[n.lines.bus0.isin(buses.index) & n.lines.bus1.isin(buses.index)], "Line")
        return nr

    def busmap_for_country(x):
        prefix = x.name[0] + x.name[1] + ' '
        logger.debug(f"Determining busmap for country {prefix[:-1]}")
        if len(x) == 1:
            return pd.Series(prefix + '0', index=x.index)
        weight = weighting_for_country(n, x)

        if algorithm == "kmeans":
            return prefix + busmap_by_kmeans(n, weight, n_clusters[x.name], buses_i=x.index, **algorithm_kwds)
        elif algorithm == "spectral":
            return prefix + busmap_by_spectral_clustering(reduce_network(n, x), n_clusters[x.name], **algorithm_kwds)
        elif algorithm == "louvain":
            return prefix + busmap_by_louvain(reduce_network(n, x), n_clusters[x.name], **algorithm_kwds)
        else:
            raise ValueError(f"`algorithm` must be one of 'kmeans', 'spectral' or 'louvain'. Is {algorithm}.")

    return (n.buses.groupby(['country', 'sub_network'], group_keys=False)
            .apply(busmap_for_country).squeeze().rename('busmap'))


def clustering_for_n_clusters(n, n_clusters, custom_busmap=False, aggregate_carriers=None,
                              line_length_factor=1.25, potential_mode='simple', solver_name="cbc",
                              algorithm="kmeans", extended_link_costs=0, focus_weights=None):

    if potential_mode == 'simple':
        p_nom_max_strategy = np.sum
    elif potential_mode == 'conservative':
        p_nom_max_strategy = np.min
    else:
        raise AttributeError(f"potential_mode should be one of 'simple' or 'conservative' but is '{potential_mode}'")

    if custom_busmap:
        busmap = pd.read_csv(snakemake.input.custom_busmap, index_col=0, squeeze=True)
        busmap.index = busmap.index.astype(str)
        logger.info(f"Imported custom busmap from {snakemake.input.custom_busmap}")
    else:
        busmap = busmap_for_n_clusters(n, n_clusters, solver_name, focus_weights, algorithm)

    clustering = get_clustering_from_busmap(
        n, busmap,
        bus_strategies=dict(country=_make_consense("Bus", "country")),
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=aggregate_carriers,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=line_length_factor,
        generator_strategies={'p_nom_max': p_nom_max_strategy, 'p_nom_min': np.sum},
        scale_link_capital_costs=False)

    if not n.links.empty:
        nc = clustering.network
        nc.links['underwater_fraction'] = (n.links.eval('underwater_fraction * length')
                                        .div(nc.links.length).dropna())
        nc.links['capital_cost'] = (nc.links['capital_cost']
                                    .add((nc.links.length - n.links.length)
                                        .clip(lower=0).mul(extended_link_costs),
                                        fill_value=0))

    return clustering


def save_to_geojson(s, fn):
    if os.path.exists(fn):
        os.unlink(fn)
    df = s.reset_index()
    schema = {**gpd.io.file.infer_schema(df), 'geometry': 'Unknown'}
    df.to_file(fn, driver='GeoJSON', schema=schema)


def cluster_regions(busmaps, input=None, output=None):
    if input is None: input = snakemake.input
    if output is None: output = snakemake.output

    busmap = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])

    for which in ('regions_onshore', 'regions_offshore'):
        regions = gpd.read_file(getattr(input, which)).set_index('name')
        geom_c = regions.geometry.groupby(busmap).apply(shapely.ops.cascaded_union)
        regions_c = gpd.GeoDataFrame(dict(geometry=geom_c))
        regions_c.index.name = 'name'
        save_to_geojson(regions_c, getattr(output, which))


def plot_busmap_for_n_clusters(n, n_clusters, fn=None):
    busmap = busmap_for_n_clusters(n, n_clusters)
    cs = busmap.unique()
    cr = sns.color_palette("hls", len(cs))
    n.plot(bus_colors=busmap.map(dict(zip(cs, cr))))
    if fn is not None:
        plt.savefig(fn, bbox_inches='tight')
    del cs, cr


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('cluster_network', network='elec', simpl='', clusters='5')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    focus_weights = snakemake.config.get('focus_weights', None)

    renewable_carriers = pd.Index([tech
                                   for tech in n.generators.carrier.unique()
                                   if tech in snakemake.config['renewable']])

    if snakemake.wildcards.clusters.endswith('m'):
        n_clusters = int(snakemake.wildcards.clusters[:-1])
        aggregate_carriers = pd.Index(n.generators.carrier.unique()).difference(renewable_carriers)
    else:
        n_clusters = int(snakemake.wildcards.clusters)
        aggregate_carriers = None # All

    if n_clusters == len(n.buses):
        # Fast-path if no clustering is necessary
        busmap = n.buses.index.to_series()
        linemap = n.lines.index.to_series()
        clustering = pypsa.networkclustering.Clustering(n, busmap, linemap, linemap, pd.Series(dtype='O'))
    else:
        line_length_factor = snakemake.config['lines']['length_factor']
        Nyears = n.snapshot_weightings.objective.sum()/8760
        hvac_overhead_cost = (load_costs(Nyears,
                                   tech_costs=snakemake.input.tech_costs,
                                   config=snakemake.config['costs'],
                                   elec_config=snakemake.config['electricity'])
                              .at['HVAC overhead', 'capital_cost'])

        def consense(x):
            v = x.iat[0]
            assert ((x == v).all() or x.isnull().all()), (
                "The `potential` configuration option must agree for all renewable carriers, for now!"
            )
            return v
        potential_mode = consense(pd.Series([snakemake.config['renewable'][tech]['potential']
                                             for tech in renewable_carriers]))
        custom_busmap = snakemake.config["enable"].get("custom_busmap", False)
        clustering = clustering_for_n_clusters(n, n_clusters, custom_busmap, aggregate_carriers,
                                               line_length_factor=line_length_factor,
                                               potential_mode=potential_mode,
                                               solver_name=snakemake.config['solving']['solver']['name'],
                                               extended_link_costs=hvac_overhead_cost,
                                               focus_weights=focus_weights)

    update_p_nom_max(n)

    clustering.network.export_to_netcdf(snakemake.output.network)
    for attr in ('busmap', 'linemap'): #also available: linemap_positive, linemap_negative
        getattr(clustering, attr).to_csv(snakemake.output[attr])

    cluster_regions((clustering.busmap,))


#make_summary
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Creates summaries of aggregated energy and costs as ``.csv`` files.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        USD2013_to_EUR2013:
        discountrate:
        marginal_cost:
        capital_cost:

    electricity:
        max_hours:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`

Inputs
------

Outputs
-------

Description
-----------

The following rule can be used to summarize the results in seperate .csv files:

.. code::

    snakemake results/summaries/elec_s_all_lall_Co2L-3H_all
                                         clusters
                                             line volume or cost cap
                                                - options
                                                        - all countries

the line volume/cost cap field can be set to one of the following:
* ``lv1.25`` for a particular line volume extension by 25%
* ``lc1.25`` for a line cost extension by 25 %
* ``lall`` for all evalutated caps
* ``lvall`` for all line volume caps
* ``lcall`` for all line cost caps

Replacing '/summaries/' with '/plots/' creates nice colored maps of the results.

"""

import logging
from _helpers import configure_logging

import os
import pypsa
import pandas as pd

from add_electricity import load_costs, update_transmission_costs

idx = pd.IndexSlice

logger = logging.getLogger(__name__)

opt_name = {"Store": "e", "Line" : "s", "Transformer" : "s"}


def _add_indexed_rows(df, raw_index):
    new_index = df.index.union(pd.MultiIndex.from_product(raw_index))
    if isinstance(new_index, pd.Index):
        new_index = pd.MultiIndex.from_tuples(new_index)

    return df.reindex(new_index)


def assign_carriers(n):

    if "carrier" not in n.loads:
        n.loads["carrier"] = "electricity"
        for carrier in ["transport","heat","urban heat"]:
            n.loads.loc[n.loads.index.str.contains(carrier),"carrier"] = carrier

    n.storage_units['carrier'].replace({'hydro': 'hydro+PHS', 'PHS': 'hydro+PHS'}, inplace=True)

    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"

    n.lines["carrier"].replace({"AC": "lines"}, inplace=True)

    if n.links.empty: n.links["carrier"] = pd.Series(dtype=str)
    n.links["carrier"].replace({"DC": "lines"}, inplace=True)

    if "EU gas store" in n.stores.index and n.stores.loc["EU gas Store","carrier"] == "":
        n.stores.loc["EU gas Store","carrier"] = "gas Store"


def calculate_costs(n, label, costs):

    for c in n.iterate_components(n.branch_components|n.controllable_one_port_components^{"Load"}):
        capital_costs = c.df.capital_cost*c.df[opt_name.get(c.name,"p") + "_nom_opt"]
        capital_costs_grouped = capital_costs.groupby(c.df.carrier).sum()

        # Index tuple(s) indicating the newly to-be-added row(s)
        raw_index = tuple([[c.list_name],["capital"],list(capital_costs_grouped.index)])
        costs = _add_indexed_rows(costs, raw_index)

        costs.loc[idx[raw_index],label] = capital_costs_grouped.values

        if c.name == "Link":
            p = c.pnl.p0.multiply(n.snapshot_weightings.generators,axis=0).sum()
        elif c.name == "Line":
            continue
        elif c.name == "StorageUnit":
            p_all = c.pnl.p.multiply(n.snapshot_weightings.generators,axis=0)
            p_all[p_all < 0.] = 0.
            p = p_all.sum()
        else:
            p = c.pnl.p.multiply(n.snapshot_weightings.generators,axis=0).sum()

        marginal_costs = p*c.df.marginal_cost

        marginal_costs_grouped = marginal_costs.groupby(c.df.carrier).sum()

        costs = costs.reindex(costs.index.union(pd.MultiIndex.from_product([[c.list_name],["marginal"],marginal_costs_grouped.index])))

        costs.loc[idx[c.list_name,"marginal",list(marginal_costs_grouped.index)],label] = marginal_costs_grouped.values

    return costs

def calculate_curtailment(n, label, curtailment):

    avail = n.generators_t.p_max_pu.multiply(n.generators.p_nom_opt).sum().groupby(n.generators.carrier).sum()
    used = n.generators_t.p.sum().groupby(n.generators.carrier).sum()

    curtailment[label] = (((avail - used)/avail)*100).round(3)

    return curtailment

def calculate_energy(n, label, energy):

    for c in n.iterate_components(n.one_port_components|n.branch_components):

        if c.name in {'Generator', 'Load', 'ShuntImpedance'}:
            c_energies = c.pnl.p.multiply(n.snapshot_weightings.generators,axis=0).sum().multiply(c.df.sign).groupby(c.df.carrier).sum()
        elif c.name in {'StorageUnit', 'Store'}:
            c_energies = c.pnl.p.multiply(n.snapshot_weightings.stores,axis=0).sum().multiply(c.df.sign).groupby(c.df.carrier).sum()
        else:
            c_energies = (-c.pnl.p1.multiply(n.snapshot_weightings.generators,axis=0).sum() - c.pnl.p0.multiply(n.snapshot_weightings.generators,axis=0).sum()).groupby(c.df.carrier).sum()

        energy = include_in_summary(energy, [c.list_name], label, c_energies)

    return energy

def include_in_summary(summary, multiindexprefix, label, item):

    # Index tuple(s) indicating the newly to-be-added row(s)
    raw_index = tuple([multiindexprefix,list(item.index)])
    summary = _add_indexed_rows(summary, raw_index)

    summary.loc[idx[raw_index], label] = item.values

    return summary

def calculate_capacity(n,label,capacity):

    for c in n.iterate_components(n.one_port_components):
        if 'p_nom_opt' in c.df.columns:
            c_capacities = abs(c.df.p_nom_opt.multiply(c.df.sign)).groupby(c.df.carrier).sum()
            capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    for c in n.iterate_components(n.passive_branch_components):
        c_capacities = c.df['s_nom_opt'].groupby(c.df.carrier).sum()
        capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    for c in n.iterate_components(n.controllable_branch_components):
        c_capacities = c.df.p_nom_opt.groupby(c.df.carrier).sum()
        capacity = include_in_summary(capacity, [c.list_name], label, c_capacities)

    return capacity

def calculate_supply(n, label, supply):
    """calculate the max dispatch of each component at the buses where the loads are attached"""

    load_types = n.loads.carrier.value_counts().index

    for i in load_types:

        buses = n.loads.bus[n.loads.carrier == i].values

        bus_map = pd.Series(False,index=n.buses.index)

        bus_map.loc[buses] = True

        for c in n.iterate_components(n.one_port_components):

            items = c.df.index[c.df.bus.map(bus_map)]

            if len(items) == 0 or c.pnl.p.empty:
                continue

            s = c.pnl.p[items].max().multiply(c.df.loc[items,'sign']).groupby(c.df.loc[items,'carrier']).sum()

            # Index tuple(s) indicating the newly to-be-added row(s)
            raw_index = tuple([[i],[c.list_name],list(s.index)])
            supply = _add_indexed_rows(supply, raw_index)

            supply.loc[idx[raw_index],label] = s.values


        for c in n.iterate_components(n.branch_components):

            for end in ["0","1"]:

                items = c.df.index[c.df["bus" + end].map(bus_map)]

                if len(items) == 0 or c.pnl["p"+end].empty:
                    continue

                #lots of sign compensation for direction and to do maximums
                s = (-1)**(1-int(end))*((-1)**int(end)*c.pnl["p"+end][items]).max().groupby(c.df.loc[items,'carrier']).sum()

                supply = supply.reindex(supply.index.union(pd.MultiIndex.from_product([[i],[c.list_name],s.index])))
                supply.loc[idx[i,c.list_name,list(s.index)],label] = s.values

    return supply


def calculate_supply_energy(n, label, supply_energy):
    """calculate the total dispatch of each component at the buses where the loads are attached"""

    load_types = n.loads.carrier.value_counts().index

    for i in load_types:

        buses = n.loads.bus[n.loads.carrier == i].values

        bus_map = pd.Series(False,index=n.buses.index)

        bus_map.loc[buses] = True

        for c in n.iterate_components(n.one_port_components):

            items = c.df.index[c.df.bus.map(bus_map)]

            if len(items) == 0 or c.pnl.p.empty:
                continue

            s = c.pnl.p[items].sum().multiply(c.df.loc[items,'sign']).groupby(c.df.loc[items,'carrier']).sum()

            # Index tuple(s) indicating the newly to-be-added row(s)
            raw_index = tuple([[i],[c.list_name],list(s.index)])
            supply_energy = _add_indexed_rows(supply_energy, raw_index)

            supply_energy.loc[idx[raw_index],label] = s.values


        for c in n.iterate_components(n.branch_components):

            for end in ["0","1"]:

                items = c.df.index[c.df["bus" + end].map(bus_map)]

                if len(items) == 0  or c.pnl['p' + end].empty:
                    continue

                s = (-1)*c.pnl["p"+end][items].sum().groupby(c.df.loc[items,'carrier']).sum()

                supply_energy = supply_energy.reindex(supply_energy.index.union(pd.MultiIndex.from_product([[i],[c.list_name],s.index])))
                supply_energy.loc[idx[i,c.list_name,list(s.index)],label] = s.values

    return supply_energy


def calculate_metrics(n,label,metrics):

    metrics = metrics.reindex(metrics.index.union(pd.Index(["line_volume","line_volume_limit","line_volume_AC","line_volume_DC","line_volume_shadow","co2_shadow"])))

    metrics.at["line_volume_DC",label] = (n.links.length*n.links.p_nom_opt)[n.links.carrier == "DC"].sum()
    metrics.at["line_volume_AC",label] = (n.lines.length*n.lines.s_nom_opt).sum()
    metrics.at["line_volume",label] = metrics.loc[["line_volume_AC","line_volume_DC"],label].sum()

    if hasattr(n,"line_volume_limit"):
        metrics.at["line_volume_limit",label] = n.line_volume_limit

    if hasattr(n,"line_volume_limit_dual"):
        metrics.at["line_volume_shadow",label] = n.line_volume_limit_dual

    if "CO2Limit" in n.global_constraints.index:
        metrics.at["co2_shadow",label] = n.global_constraints.at["CO2Limit","mu"]

    return metrics


def calculate_prices(n,label,prices):

    bus_type = pd.Series(n.buses.index.str[3:],n.buses.index).replace("","electricity")

    prices = prices.reindex(prices.index.union(bus_type.value_counts().index))

    logger.warning("Prices are time-averaged, not load-weighted")
    prices[label] = n.buses_t.marginal_price.mean().groupby(bus_type).mean()

    return prices


def calculate_weighted_prices(n,label,weighted_prices):

    logger.warning("Weighted prices don't include storage units as loads")

    weighted_prices = weighted_prices.reindex(pd.Index(["electricity","heat","space heat","urban heat","space urban heat","gas","H2"]))

    link_loads = {"electricity" :  ["heat pump", "resistive heater", "battery charger", "H2 Electrolysis"],
                  "heat" : ["water tanks charger"],
                  "urban heat" : ["water tanks charger"],
                  "space heat" : [],
                  "space urban heat" : [],
                  "gas" : ["OCGT","gas boiler","CHP electric","CHP heat"],
                  "H2" : ["Sabatier", "H2 Fuel Cell"]}

    for carrier in link_loads:

        if carrier == "electricity":
            suffix = ""
        elif carrier[:5] == "space":
            suffix = carrier[5:]
        else:
            suffix =  " " + carrier

        buses = n.buses.index[n.buses.index.str[2:] == suffix]

        if buses.empty:
            continue

        if carrier in ["H2","gas"]:
            load = pd.DataFrame(index=n.snapshots,columns=buses,data=0.)
        elif carrier[:5] == "space":
            load = heat_demand_df[buses.str[:2]].rename(columns=lambda i: str(i)+suffix)
        else:
            load = n.loads_t.p_set[buses]


        for tech in link_loads[carrier]:

            names = n.links.index[n.links.index.to_series().str[-len(tech):] == tech]

            if names.empty:
                continue

            load += n.links_t.p0[names].groupby(n.links.loc[names,"bus0"],axis=1).sum(axis=1)

        # Add H2 Store when charging
        if carrier == "H2":
            stores = n.stores_t.p[buses+ " Store"].groupby(n.stores.loc[buses+ " Store","bus"],axis=1).sum(axis=1)
            stores[stores > 0.] = 0.
            load += -stores

        weighted_prices.loc[carrier,label] = (load*n.buses_t.marginal_price[buses]).sum().sum()/load.sum().sum()

        if carrier[:5] == "space":
            print(load*n.buses_t.marginal_price[buses])

    return weighted_prices


outputs = ["costs",
           "curtailment",
           "energy",
           "capacity",
           "supply",
           "supply_energy",
           "prices",
           "weighted_prices",
           "metrics",
           ]


def make_summaries(networks_dict, country='all'):

    columns = pd.MultiIndex.from_tuples(networks_dict.keys(),names=["simpl","clusters","ll","opts"])

    dfs = {}

    for output in outputs:
        dfs[output] = pd.DataFrame(columns=columns,dtype=float)

    for label, filename in networks_dict.items():
        print(label, filename)
        if not os.path.exists(filename):
            print("does not exist!!")
            continue

        try:
            n = pypsa.Network(filename)
        except OSError:
            logger.warning("Skipping {filename}".format(filename=filename))
            continue

        if country != 'all':
            n = n[n.buses.country == country]

        Nyears = n.snapshot_weightings.objective.sum() / 8760.
        costs = load_costs(Nyears, snakemake.input[0],
                           snakemake.config['costs'], snakemake.config['electricity'])
        update_transmission_costs(n, costs, simple_hvdc_costs=False)

        assign_carriers(n)

        for output in outputs:
            dfs[output] = globals()["calculate_" + output](n, label, dfs[output])

    return dfs


def to_csv(dfs):
    dir = snakemake.output[0]
    os.makedirs(dir, exist_ok=True)
    for key, df in dfs.items():
        df.to_csv(os.path.join(dir, f"{key}.csv"))


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('make_summary', network='elec', simpl='',
                           clusters='5', ll='copt', opts='Co2L-24H', country='all')
        network_dir = os.path.join('..', 'results', 'networks')
    else:
        network_dir = os.path.join('results', 'networks')
    configure_logging(snakemake)

    def expand_from_wildcard(key):
        w = getattr(snakemake.wildcards, key)
        return snakemake.config["scenario"][key] if w == "all" else [w]

    if snakemake.wildcards.ll.endswith("all"):
        ll = snakemake.config["scenario"]["ll"]
        if len(snakemake.wildcards.ll) == 4:
            ll = [l for l in ll if l[0] == snakemake.wildcards.ll[0]]
    else:
        ll = [snakemake.wildcards.ll]

    networks_dict = {(simpl,clusters,l,opts) :
        os.path.join(network_dir, f'elec_s{simpl}_'
                                  f'{clusters}_ec_l{l}_{opts}.nc')
                     for simpl in expand_from_wildcard("simpl")
                     for clusters in expand_from_wildcard("clusters")
                     for l in ll
                     for opts in expand_from_wildcard("opts")}

    dfs = make_summaries(networks_dict, country=snakemake.wildcards.country)

    to_csv(dfs)

## override_components

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 13:08:35 2022

@author: mdomenech

script to override components
(links with more than one input/output)
data needed
"""
import os
import pandas as pd
from pypsa.descriptors import Dict
from pypsa.components import components, component_attrs

import logging
logger = logging.getLogger(__name__)

def override_component_attrs(directory):
    """Tell PyPSA that links can have multiple outputs by
    overriding the component_attrs. This can be done for
    as many buses as you need with format busi for i = 2,3,4,5,....
    See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs

    Parameters
    ----------
    directory : string
        Folder where component attributes to override are stored
        analogous to ``pypsa/component_attrs``, e.g. `links.csv`.

    Returns
    -------
    Dictionary of overriden component attributes.
    """

    attrs = Dict({k : v.copy() for k,v in component_attrs.items()})

    for component, list_name in components.list_name.items():
        fn = f"{directory}/{list_name}.csv"
        if os.path.isfile(fn):
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs

## plot_network
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Plots map with pie charts and cost box bar charts.

Relevant Settings
-----------------

Inputs
------

Outputs
-------

Description
-----------

"""

import logging
from _helpers import (load_network_for_plots, aggregate_p, aggregate_costs,
                      configure_logging)

import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, Ellipse
from matplotlib.legend_handler import HandlerPatch
to_rgba = mpl.colors.colorConverter.to_rgba

logger = logging.getLogger(__name__)


def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    fig = ax.get_figure()
    def axes2pt():
        return np.diff(ax.transData.transform([(0,0), (1,1)]), axis=0)[0] * (72./fig.dpi)

    ellipses = []
    if not dont_resize_actively:
        def update_width_height(event):
            dist = axes2pt()
            for e, radius in ellipses: e.width, e.height = 2. * radius * dist
        fig.canvas.mpl_connect('resize_event', update_width_height)
        ax.callbacks.connect('xlim_changed', update_width_height)
        ax.callbacks.connect('ylim_changed', update_width_height)

    def legend_circle_handler(legend, orig_handle, xdescent, ydescent,
                              width, height, fontsize):
        w, h = 2. * orig_handle.get_radius() * axes2pt()
        e = Ellipse(xy=(0.5*width-0.5*xdescent, 0.5*height-0.5*ydescent), width=w, height=w)
        ellipses.append((e, orig_handle.get_radius()))
        return e
    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}


def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0,0), radius=(s/scale)**0.5, **kw) for s in sizes]


def set_plot_style():
    plt.style.use(['classic', 'seaborn-white',
                {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 12,
                    'legend.fontsize': 'medium',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                }])


def plot_map(n, ax=None, attribute='p_nom', opts={}):
    if ax is None:
        ax = plt.gca()

    ## DATA
    line_colors = {'cur': "purple",
                   'exp': mpl.colors.rgb2hex(to_rgba("red", 0.7), True)}
    tech_colors = opts['tech_colors']

    if attribute == 'p_nom':
        # bus_sizes = n.generators_t.p.sum().loc[n.generators.carrier == "load"].groupby(n.generators.bus).sum()
        bus_sizes = pd.concat((n.generators.query('carrier != "load"').groupby(['bus', 'carrier']).p_nom_opt.sum(),
                               n.storage_units.groupby(['bus', 'carrier']).p_nom_opt.sum()))
        line_widths_exp = n.lines.s_nom_opt
        line_widths_cur = n.lines.s_nom_min
        link_widths_exp = n.links.p_nom_opt
        link_widths_cur = n.links.p_nom_min
    else:
        raise 'plotting of {} has not been implemented yet'.format(attribute)


    line_colors_with_alpha = \
        ((line_widths_cur / n.lines.s_nom > 1e-3)
         .map({True: line_colors['cur'], False: to_rgba(line_colors['cur'], 0.)}))
    link_colors_with_alpha = \
        ((link_widths_cur / n.links.p_nom > 1e-3)
        .map({True: line_colors['cur'], False: to_rgba(line_colors['cur'], 0.)}))


    ## FORMAT
    linewidth_factor = opts['map'][attribute]['linewidth_factor']
    bus_size_factor  = opts['map'][attribute]['bus_size_factor']

    ## PLOT
    n.plot(line_widths=line_widths_exp/linewidth_factor,
           link_widths=link_widths_exp/linewidth_factor,
           line_colors=line_colors['exp'],
           link_colors=line_colors['exp'],
           bus_sizes=bus_sizes/bus_size_factor,
           bus_colors=tech_colors,
           boundaries=map_boundaries,
           color_geomap=True, geomap=True,
           ax=ax)
    n.plot(line_widths=line_widths_cur/linewidth_factor,
           link_widths=link_widths_cur/linewidth_factor,
           line_colors=line_colors_with_alpha,
           link_colors=link_colors_with_alpha,
           bus_sizes=0,
           boundaries=map_boundaries,
           color_geomap=True, geomap=False,
           ax=ax)
    ax.set_aspect('equal')
    ax.axis('off')

    # Rasterize basemap
    # TODO : Check if this also works with cartopy
    for c in ax.collections[:2]: c.set_rasterized(True)

    # LEGEND
    handles = []
    labels = []

    for s in (10, 1):
        handles.append(plt.Line2D([0],[0],color=line_colors['exp'],
                                linewidth=s*1e3/linewidth_factor))
        labels.append("{} GW".format(s))
    l1_1 = ax.legend(handles, labels,
                     loc="upper left", bbox_to_anchor=(0.24, 1.01),
                     frameon=False,
                     labelspacing=0.8, handletextpad=1.5,
                     title='Transmission Exp./Exist.             ')
    ax.add_artist(l1_1)

    handles = []
    labels = []
    for s in (10, 5):
        handles.append(plt.Line2D([0],[0],color=line_colors['cur'],
                                linewidth=s*1e3/linewidth_factor))
        labels.append("/")
    l1_2 = ax.legend(handles, labels,
                loc="upper left", bbox_to_anchor=(0.26, 1.01),
                frameon=False,
                labelspacing=0.8, handletextpad=0.5,
                title=' ')
    ax.add_artist(l1_2)

    handles = make_legend_circles_for([10e3, 5e3, 1e3], scale=bus_size_factor, facecolor="w")
    labels = ["{} GW".format(s) for s in (10, 5, 3)]
    l2 = ax.legend(handles, labels,
                loc="upper left", bbox_to_anchor=(0.01, 1.01),
                frameon=False, labelspacing=1.0,
                title='Generation',
                handler_map=make_handler_map_to_scale_circles_as_in(ax))
    ax.add_artist(l2)

    techs =  (bus_sizes.index.levels[1]).intersection(pd.Index(opts['vre_techs'] + opts['conv_techs'] + opts['storage_techs']))
    handles = []
    labels = []
    for t in techs:
        handles.append(plt.Line2D([0], [0], color=tech_colors[t], marker='o', markersize=8, linewidth=0))
        labels.append(opts['nice_names'].get(t, t))
    l3 = ax.legend(handles, labels, loc="upper center",  bbox_to_anchor=(0.5, -0.), # bbox_to_anchor=(0.72, -0.05),
                handletextpad=0., columnspacing=0.5, ncol=4, title='Technology')

    return fig


def plot_total_energy_pie(n, ax=None):
    if ax is None: ax = plt.gca()

    ax.set_title('Energy per technology', fontdict=dict(fontsize="medium"))

    e_primary = aggregate_p(n).drop('load', errors='ignore').loc[lambda s: s>0]

    patches, texts, autotexts = ax.pie(e_primary,
        startangle=90,
        labels = e_primary.rename(opts['nice_names']).index,
        autopct='%.0f%%',
        shadow=False,
        colors = [opts['tech_colors'][tech] for tech in e_primary.index])
    for t1, t2, i in zip(texts, autotexts, e_primary.index):
        if e_primary.at[i] < 0.04 * e_primary.sum():
            t1.remove()
            t2.remove()

def plot_total_cost_bar(n, ax=None):
    if ax is None: ax = plt.gca()

    total_load = (n.snapshot_weightings.generators * n.loads_t.p.sum(axis=1)).sum()
    tech_colors = opts['tech_colors']

    def split_costs(n):
        costs = aggregate_costs(n).reset_index(level=0, drop=True)
        costs_ex = aggregate_costs(n, existing_only=True).reset_index(level=0, drop=True)
        return (costs['capital'].add(costs['marginal'], fill_value=0.),
                costs_ex['capital'], costs['capital'] - costs_ex['capital'], costs['marginal'])

    costs, costs_cap_ex, costs_cap_new, costs_marg = split_costs(n)

    costs_graph = pd.DataFrame(dict(a=costs.drop('load', errors='ignore')),
                            index=['AC-AC', 'AC line', 'onwind', 'offwind-ac',
                                   'offwind-dc', 'solar', 'OCGT','CCGT', 'battery', 'H2']).dropna()
    bottom = np.array([0., 0.])
    texts = []

    for i,ind in enumerate(costs_graph.index):
        data = np.asarray(costs_graph.loc[ind])/total_load
        ax.bar([0.5], data, bottom=bottom, color=tech_colors[ind],
               width=0.7, zorder=-1)
        bottom_sub = bottom
        bottom = bottom+data

        if ind in opts['conv_techs'] + ['AC line']:
            for c in [costs_cap_ex, costs_marg]:
                if ind in c:
                    data_sub = np.asarray([c.loc[ind]])/total_load
                    ax.bar([0.5], data_sub, linewidth=0,
                        bottom=bottom_sub, color=tech_colors[ind],
                        width=0.7, zorder=-1, alpha=0.8)
                    bottom_sub += data_sub

        if abs(data[-1]) < 5:
            continue

        text = ax.text(1.1,(bottom-0.5*data)[-1]-3,opts['nice_names'].get(ind,ind))
        texts.append(text)

    ax.set_ylabel("Average system cost [Eur/MWh]")
    ax.set_ylim([0, opts.get('costs_max', 80)])
    ax.set_xlim([0, 1])
    ax.set_xticklabels([])
    ax.grid(True, axis="y", color='k', linestyle='dotted')


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_network', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-24H',
                                  attr='p_nom', ext="pdf")
    configure_logging(snakemake)

    set_plot_style()

    opts = snakemake.config['plotting']
    map_figsize = opts['map']['figsize']
    map_boundaries = opts['map']['boundaries']

    n = load_network_for_plots(snakemake.input.network, snakemake.input.tech_costs, snakemake.config)

    scenario_opts = snakemake.wildcards.opts.split('-')

    fig, ax = plt.subplots(figsize=map_figsize, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_map(n, ax, snakemake.wildcards.attr, opts)

    fig.savefig(snakemake.output.only_map, dpi=150, bbox_inches='tight')

    ax1 = fig.add_axes([-0.115, 0.625, 0.2, 0.2])
    plot_total_energy_pie(n, ax1)

    ax2 = fig.add_axes([-0.075, 0.1, 0.1, 0.45])
    plot_total_cost_bar(n, ax2)

    ll = snakemake.wildcards.ll
    ll_type = ll[0]
    ll_factor = ll[1:]
    lbl = dict(c='line cost', v='line volume')[ll_type]
    amnt = '{ll} x today\'s'.format(ll=ll_factor) if ll_factor != 'opt' else 'optimal'
    fig.suptitle('Expansion to {amount} {label} at {clusters} clusters'
                .format(amount=amnt, label=lbl, clusters=snakemake.wildcards.clusters))

    fig.savefig(snakemake.output.ext, transparent=True, bbox_inches='tight')

## plot_p_nom_max
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Plots renewable installation potentials per capacity factor.

Relevant Settings
-----------------

Inputs
------

Outputs
-------

Description
-----------

"""
import logging
from _helpers import configure_logging

import pypsa
import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def cum_p_nom_max(net, tech, country=None):
    carrier_b = net.generators.carrier == tech

    generators = pd.DataFrame(dict(
        p_nom_max=net.generators.loc[carrier_b, 'p_nom_max'],
        p_max_pu=net.generators_t.p_max_pu.loc[:,carrier_b].mean(),
        country=net.generators.loc[carrier_b, 'bus'].map(net.buses.country)
    )).sort_values("p_max_pu", ascending=False)

    if country is not None:
        generators = generators.loc[generators.country == country]

    generators["cum_p_nom_max"] = generators["p_nom_max"].cumsum() / 1e6

    return generators


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_p_nom_max', network='elec', simpl='',
                                  techs='solar,onwind,offwind-dc', ext='png',
                                  clusts= '5,full', country= 'all')
    configure_logging(snakemake)

    plot_kwds = dict(drawstyle="steps-post")

    clusters = snakemake.wildcards.clusts.split(',')
    techs = snakemake.wildcards.techs.split(',')
    country = snakemake.wildcards.country
    if country == 'all':
        country = None
    else:
        plot_kwds['marker'] = 'x'

    fig, axes = plt.subplots(1, len(techs))

    for j, cluster in enumerate(clusters):
        net = pypsa.Network(snakemake.input[j])

        for i, tech in enumerate(techs):
            cum_p_nom_max(net, tech, country).plot(x="p_max_pu", y="cum_p_nom_max",
                         label=cluster, ax=axes[i], **plot_kwds)

    for i, tech in enumerate(techs):
        ax = axes[i]
        ax.set_xlabel(f"Capacity factor of {tech}")
        ax.set_ylabel("Cumulative installable capacity / TW")

    plt.legend(title="Cluster level")

    fig.savefig(snakemake.output[0], transparent=True, bbox_inches='tight')

##plot_summary

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Plots energy and cost summaries for solved networks.

Relevant Settings
-----------------

Inputs
------

Outputs
-------

Description
-----------

"""

import os
import logging
from _helpers import configure_logging

import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def rename_techs(label):
    if "H2" in label:
        label = "hydrogen storage"
    elif label == "solar":
        label = "solar PV"
    elif label == "offwind-ac":
        label = "offshore wind ac"
    elif label == "offwind-dc":
        label = "offshore wind dc"
    elif label == "onwind":
        label = "onshore wind"
    elif label == "ror":
        label = "hydroelectricity"
    elif label == "hydro":
        label = "hydroelectricity"
    elif label == "PHS":
        label = "hydroelectricity"
    elif "battery" in label:
        label = "battery storage"

    return label


preferred_order = pd.Index(["transmission lines","hydroelectricity","hydro reservoir","run of river","pumped hydro storage","onshore wind","offshore wind ac", "offshore wind dc","solar PV","solar thermal","OCGT","hydrogen storage","battery storage"])


def plot_costs(infn, fn=None):

    ## For now ignore the simpl header
    cost_df = pd.read_csv(infn,index_col=list(range(3)),header=[1,2,3])

    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()

    #convert to billions
    df = df/1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.config['plotting']['costs_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    new_index = (preferred_order&df.index).append(df.index.difference(preferred_order))

    new_columns = df.sum().sort_values().index

    fig, ax = plt.subplots()
    fig.set_size_inches((12,8))

    df.loc[new_index,new_columns].T.plot(kind="bar",ax=ax,stacked=True,color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index])


    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([0,snakemake.config['plotting']['costs_max']])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles,labels,ncol=4,loc="upper left")


    fig.tight_layout()

    if fn is not None:
        fig.savefig(fn, transparent=True)


def plot_energy(infn, fn=None):

    energy_df = pd.read_csv(infn, index_col=list(range(2)),header=[1,2,3])

    df = energy_df.groupby(energy_df.index.get_level_values(1)).sum()

    #convert MWh to TWh
    df = df/1e6

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.abs().max(axis=1) < snakemake.config['plotting']['energy_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    new_index = (preferred_order&df.index).append(df.index.difference(preferred_order))

    new_columns = df.columns.sort_values()

    fig, ax = plt.subplots()
    fig.set_size_inches((12,8))

    df.loc[new_index,new_columns].T.plot(kind="bar",ax=ax,stacked=True,color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index])


    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([snakemake.config['plotting']['energy_min'],snakemake.config['plotting']['energy_max']])

    ax.set_ylabel("Energy [TWh/a]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles,labels,ncol=4,loc="upper left")


    fig.tight_layout()

    if fn is not None:
        fig.savefig(fn, transparent=True)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_summary', summary='energy', network='elec',
                                  simpl='', clusters=5, ll='copt', opts='Co2L-24H',
                                  attr='', ext='png', country='all')
    configure_logging(snakemake)

    summary = snakemake.wildcards.summary
    try:
        func = globals()[f"plot_{summary}"]
    except KeyError:
        raise RuntimeError(f"plotting function for {summary} has not been defined")

    func(os.path.join(snakemake.input[0], f"{summary}.csv"), snakemake.output[0])

## prepare_links_p_nom

#!/usr/bin/env python

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Extracts capacities of HVDC links from `Wikipedia <https://en.wikipedia.org/wiki/List_of_HVDC_projects>`_.

Relevant Settings
-----------------

.. code:: yaml

    enable:
        prepare_links_p_nom:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

Inputs
------

*None*

Outputs
-------

- ``data/links_p_nom.csv``: A plain download of https://en.wikipedia.org/wiki/List_of_HVDC_projects#Europe plus extracted coordinates.

Description
-----------

*None*

"""

import logging
from _helpers import configure_logging

import pandas as pd

logger = logging.getLogger(__name__)


def multiply(s):
    return s.str[0].astype(float) * s.str[1].astype(float)


def extract_coordinates(s):
    regex = (r"(\d{1,2})°(\d{1,2})′(\d{1,2})″(N|S) "
             r"(\d{1,2})°(\d{1,2})′(\d{1,2})″(E|W)")
    e = s.str.extract(regex, expand=True)
    lat = (e[0].astype(float) + (e[1].astype(float) + e[2].astype(float)/60.)/60.)*e[3].map({'N': +1., 'S': -1.})
    lon = (e[4].astype(float) + (e[5].astype(float) + e[6].astype(float)/60.)/60.)*e[7].map({'E': +1., 'W': -1.})
    return lon, lat


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake #rule must be enabled in config
        snakemake = mock_snakemake('prepare_links_p_nom', simpl='', network='elec')
    configure_logging(snakemake)

    links_p_nom = pd.read_html('https://en.wikipedia.org/wiki/List_of_HVDC_projects', header=0, match="SwePol")[0]

    mw = "Power (MW)"
    m_b = links_p_nom[mw].str.contains('x').fillna(False)

    links_p_nom.loc[m_b, mw] = links_p_nom.loc[m_b, mw].str.split('x').pipe(multiply)
    links_p_nom[mw] = links_p_nom[mw].str.extract("[-/]?([\d.]+)", expand=False).astype(float)

    links_p_nom['x1'], links_p_nom['y1'] = extract_coordinates(links_p_nom['Converterstation 1'])
    links_p_nom['x2'], links_p_nom['y2'] = extract_coordinates(links_p_nom['Converterstation 2'])

    links_p_nom.dropna(subset=['x1', 'y1', 'x2', 'y2']).to_csv(snakemake.output[0], index=False)


## prepare_network
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Prepare PyPSA network for solving according to :ref:`opts` and :ref:`ll`, such as

- adding an annual **limit** of carbon-dioxide emissions,
- adding an exogenous **price** per tonne emissions of carbon-dioxide (or other kinds),
- setting an **N-1 security margin** factor for transmission line capacities,
- specifying an expansion limit on the **cost** of transmission expansion,
- specifying an expansion limit on the **volume** of transmission expansion, and
- reducing the **temporal** resolution by averaging over multiple hours
  or segmenting time series into chunks of varying lengths using ``tsam``.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        emission_prices:
        USD2013_to_EUR2013:
        discountrate:
        marginal_cost:
        capital_cost:

    electricity:
        co2limit:
        max_hours:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``networks/elec_s{simpl}_{clusters}.nc``: confer :ref:`cluster`

Outputs
-------

- ``networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Complete PyPSA network that will be handed to the ``solve_network`` rule.

Description
-----------

.. tip::
    The rule :mod:`prepare_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`prepare_network`.

"""

import logging
from _helpers import configure_logging

import re
import pypsa
import numpy as np
import pandas as pd

from add_electricity import load_costs, update_transmission_costs

idx = pd.IndexSlice

logger = logging.getLogger(__name__)
from override_components import override_component_attrs #new


def add_co2limit(n, Nyears=1., factor=None):

    if factor is not None:
        annual_emissions = factor*snakemake.config['electricity']['co2base']
    else:
        annual_emissions = snakemake.config['electricity']['co2limit']

    n.add("GlobalConstraint", "CO2Limit",
          carrier_attribute="co2_emissions", sense="<=",
          constant=annual_emissions * Nyears)


def add_emission_prices(n, emission_prices=None, exclude_co2=False):
    if emission_prices is None:
        emission_prices = snakemake.config['costs']['emission_prices']
    if exclude_co2: emission_prices.pop('co2')
    ep = (pd.Series(emission_prices).rename(lambda x: x+'_emissions') *n.carriers.filter(like='_emissions')).sum(axis=1)
    gen_ep = n.generators.carrier.map(ep) / n.generators.efficiency
    n.generators['marginal_cost'] += gen_ep
    su_ep = n.storage_units.carrier.map(ep) / n.storage_units.efficiency_dispatch
    n.storage_units['marginal_cost'] += su_ep


def set_line_s_max_pu(n):
    s_max_pu = snakemake.config['lines']['s_max_pu']
    n.lines['s_max_pu'] = s_max_pu
    logger.info(f"N-1 security margin of lines set to {s_max_pu}")


def set_transmission_limit(n, ll_type, factor, Nyears=1):
    links_dc_b = n.links.carrier == 'DC' if not n.links.empty else pd.Series()

    _lines_s_nom = (np.sqrt(3) * n.lines.type.map(n.line_types.i_nom) *
                   n.lines.num_parallel *  n.lines.bus0.map(n.buses.v_nom))
    lines_s_nom = n.lines.s_nom.where(n.lines.type == '', _lines_s_nom)


    col = 'capital_cost' if ll_type == 'c' else 'length'
    ref = (lines_s_nom @ n.lines[col] +
           n.links.loc[links_dc_b, "p_nom"] @ n.links.loc[links_dc_b, col])

    costs = load_costs(Nyears, snakemake.input.tech_costs,
                       snakemake.config['costs'],
                       snakemake.config['electricity'])
    update_transmission_costs(n, costs, simple_hvdc_costs=False)

    if factor == 'opt' or float(factor) > 1.0:
        n.lines['s_nom_min'] = lines_s_nom
        n.lines['s_nom_extendable'] = True

        n.links.loc[links_dc_b, 'p_nom_min'] = n.links.loc[links_dc_b, 'p_nom']
        n.links.loc[links_dc_b, 'p_nom_extendable'] = True

    if factor != 'opt':
        con_type = 'expansion_cost' if ll_type == 'c' else 'volume_expansion'
        rhs = float(factor) * ref
        n.add('GlobalConstraint', f'l{ll_type}_limit',
              type=f'transmission_{con_type}_limit',
              sense='<=', constant=rhs, carrier_attribute='AC, DC')

    return n


def average_every_nhours(n, offset):
    logger.info(f"Resampling the network to {offset}")
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name+"_t")
        for k, df in c.pnl.items():
            if not df.empty:
                pnl[k] = df.resample(offset).mean()

    return m


def apply_time_segmentation(n, segments):
    logger.info(f"Aggregating time series to {segments} segments.")
    try:
        import tsam.timeseriesaggregation as tsam
    except:
        raise ModuleNotFoundError("Optional dependency 'tsam' not found."
                                  "Install via 'pip install tsam'")

    p_max_pu_norm = n.generators_t.p_max_pu.max()
    p_max_pu = n.generators_t.p_max_pu / p_max_pu_norm

    load_norm = n.loads_t.p_set.max()
    load = n.loads_t.p_set / load_norm

    inflow_norm = n.storage_units_t.inflow.max()
    inflow = n.storage_units_t.inflow / inflow_norm

    raw = pd.concat([p_max_pu, load, inflow], axis=1, sort=False)

    solver_name = snakemake.config["solving"]["solver"]["name"]

    agg = tsam.TimeSeriesAggregation(raw, hoursPerPeriod=len(raw),
                                     noTypicalPeriods=1, noSegments=int(segments),
                                     segmentation=True, solver=solver_name)

    segmented = agg.createTypicalPeriods()

    weightings = segmented.index.get_level_values("Segment Duration")
    offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
    snapshots = [n.snapshots[0] + pd.Timedelta(f"{offset}h") for offset in offsets]

    n.set_snapshots(pd.DatetimeIndex(snapshots, name='name'))
    n.snapshot_weightings = pd.Series(weightings, index=snapshots, name="weightings", dtype="float64")

    segmented.index = snapshots
    n.generators_t.p_max_pu = segmented[n.generators_t.p_max_pu.columns] * p_max_pu_norm
    n.loads_t.p_set = segmented[n.loads_t.p_set.columns] * load_norm
    n.storage_units_t.inflow = segmented[n.storage_units_t.inflow.columns] * inflow_norm

    return n

def enforce_autarky(n, only_crossborder=False):
    if only_crossborder:
        lines_rm = n.lines.loc[
                        n.lines.bus0.map(n.buses.country) !=
                        n.lines.bus1.map(n.buses.country)
                    ].index
        links_rm = n.links.loc[
                        n.links.bus0.map(n.buses.country) !=
                        n.links.bus1.map(n.buses.country)
                    ].index
    else:
        lines_rm = n.lines.index
        links_rm = n.links.loc[n.links.carrier=="DC"].index
    n.mremove("Line", lines_rm)
    n.mremove("Link", links_rm)

def set_line_nom_max(n):
    s_nom_max_set = snakemake.config["lines"].get("s_nom_max,", np.inf)
    p_nom_max_set = snakemake.config["links"].get("p_nom_max", np.inf)
    n.lines.s_nom_max.clip(upper=s_nom_max_set, inplace=True)
    n.links.p_nom_max.clip(upper=p_nom_max_set, inplace=True)

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('prepare_network', network='elec', simpl='',
                                  clusters='40', ll='v0.3', opts='Co2L-24H')
    configure_logging(snakemake)

    opts = snakemake.wildcards.opts.split('-')

    #n = pypsa.Network(snakemake.input[0])
    overrides = override_component_attrs(snakemake.input.overrides)   #new
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides) #new

    Nyears = n.snapshot_weightings.objective.sum() / 8760.

    set_line_s_max_pu(n)

    for o in opts:
        m = re.match(r'^\d+h$', o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break

    for o in opts:
        m = re.match(r'^\d+seg$', o, re.IGNORECASE)
        if m is not None:
            n = apply_time_segmentation(n, m.group(0)[:-3])
            break

    for o in opts:
        if "Co2L" in o:
            m = re.findall("[0-9]*\.?[0-9]+$", o)
            if len(m) > 0:
                add_co2limit(n, Nyears, float(m[0]))
            else:
                add_co2limit(n, Nyears)
            break

    for o in opts:
        oo = o.split("+")
        suptechs = map(lambda c: c.split("-", 2)[0], n.carriers.index)
        if oo[0].startswith(tuple(suptechs)):
            carrier = oo[0]
            # handles only p_nom_max as stores and lines have no potentials
            attr_lookup = {"p": "p_nom_max", "c": "capital_cost"}
            attr = attr_lookup[oo[1][0]]
            factor = float(oo[1][1:])
            if carrier == "AC":  # lines do not have carrier
                n.lines[attr] *= factor
            else:
                comps = {"Generator", "Link", "StorageUnit", "Store"}
                for c in n.iterate_components(comps):
                    sel = c.df.carrier.str.contains(carrier)
                    c.df.loc[sel,attr] *= factor

    if 'Ep' in opts:
        add_emission_prices(n)

    ll_type, factor = snakemake.wildcards.ll[0], snakemake.wildcards.ll[1:]
    set_transmission_limit(n, ll_type, factor, Nyears)

    set_line_nom_max(n)

    if "ATK" in opts:
        enforce_autarky(n)
    elif "ATKc" in opts:
        enforce_autarky(n, only_crossborder=True)

    n.export_to_netcdf(snakemake.output[0])

## retrieve_databundle
# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517935.svg
   :target: https://doi.org/10.5281/zenodo.3517935

The data bundle (1.4 GB) contains common GIS datasets like NUTS3 shapes, EEZ shapes, CORINE Landcover, Natura 2000 and also electricity specific summary statistics like historic per country yearly totals of hydro generation, GDP and POP on NUTS3 levels and per-country load time-series.

This rule downloads the data bundle from `zenodo <https://doi.org/10.5281/zenodo.3517935>`_ and extracts it in the ``data`` sub-directory, such that all files of the bundle are stored in the ``data/bundle`` subdirectory.

The :ref:`tutorial` uses a smaller `data bundle <https://zenodo.org/record/3517921/files/pypsa-eur-tutorial-data-bundle.tar.xz>`_ than required for the full model (19 MB)

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517921.svg
    :target: https://doi.org/10.5281/zenodo.3517921

**Relevant Settings**

.. code:: yaml

    tutorial:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``cutouts/bundle``: input data collected from various sources

"""

import logging
from _helpers import progress_retrieve, configure_logging

import tarfile
from pathlib import Path

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('retrieve_databundle')
        rootpath = '..'
    else:
        rootpath = '.'
    configure_logging(snakemake) # TODO Make logging compatible with progressbar (see PR #102)

    if snakemake.config['tutorial']:
        url = "https://zenodo.org/record/3517921/files/pypsa-eur-tutorial-data-bundle.tar.xz"
    else:
        url = "https://zenodo.org/record/3517935/files/pypsa-eur-data-bundle.tar.xz"

    # Save locations
    tarball_fn = Path(f"{rootpath}/bundle.tar.xz")
    to_fn = Path(f"{rootpath}/data")

    logger.info(f"Downloading databundle from '{url}'.")
    progress_retrieve(url, tarball_fn)

    logger.info(f"Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")


## simplify_network
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Lifts electrical transmission network to a single 380 kV voltage layer,
removes dead-ends of the network,
and reduces multi-hop HVDC connections to a single link.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        USD2013_to_EUR2013:
        discountrate:
        marginal_cost:
        capital_cost:

    electricity:
        max_hours:

    renewables: (keys)
        {technology}:
            potential:

    lines:
        length_factor:

    links:
        p_max_pu:

    solving:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`, :ref:`renewable_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`solving_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: confer :ref:`busregions`
- ``networks/elec.nc``: confer :ref:`electricity`

Outputs
-------

- ``resources/regions_onshore_elec_s{simpl}.geojson``:

    .. image:: ../img/regions_onshore_elec_s.png
            :scale: 33 %

- ``resources/regions_offshore_elec_s{simpl}.geojson``:

    .. image:: ../img/regions_offshore_elec_s  .png
            :scale: 33 %

- ``resources/busmap_elec_s{simpl}.csv``: Mapping of buses from ``networks/elec.nc`` to ``networks/elec_s{simpl}.nc``;
- ``networks/elec_s{simpl}.nc``:

    .. image:: ../img/elec_s.png
        :scale: 33 %

Description
-----------

The rule :mod:`simplify_network` does up to four things:

1. Create an equivalent transmission network in which all voltage levels are mapped to the 380 kV level by the function ``simplify_network(...)``.

2. DC only sub-networks that are connected at only two buses to the AC network are reduced to a single representative link in the function ``simplify_links(...)``. The components attached to buses in between are moved to the nearest endpoint. The grid connection cost of offshore wind generators are added to the captial costs of the generator.

3. Stub lines and links, i.e. dead-ends of the network, are sequentially removed from the network in the function ``remove_stubs(...)``. Components are moved along.

4. Optionally, if an integer were provided for the wildcard ``{simpl}`` (e.g. ``networks/elec_s500.nc``), the network is clustered to this number of clusters with the routines from the ``cluster_network`` rule with the function ``cluster_network.cluster(...)``. This step is usually skipped!
"""

import logging
from _helpers import configure_logging, update_p_nom_max

from cluster_network import clustering_for_n_clusters, cluster_regions
from add_electricity import load_costs

import pandas as pd
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import connected_components, dijkstra

from functools import reduce

import pypsa
from pypsa.io import import_components_from_dataframe, import_series_from_dataframe
from pypsa.networkclustering import busmap_by_stubs, aggregategenerators, aggregateoneport, get_clustering_from_busmap, _make_consense

logger = logging.getLogger(__name__)


def simplify_network_to_380(n):
    ## All goes to v_nom == 380
    logger.info("Mapping all network lines onto a single 380kV layer")

    n.buses['v_nom'] = 380.

    linetype_380, = n.lines.loc[n.lines.v_nom == 380., 'type'].unique()
    lines_v_nom_b = n.lines.v_nom != 380.
    n.lines.loc[lines_v_nom_b, 'num_parallel'] *= (n.lines.loc[lines_v_nom_b, 'v_nom'] / 380.)**2
    n.lines.loc[lines_v_nom_b, 'v_nom'] = 380.
    n.lines.loc[lines_v_nom_b, 'type'] = linetype_380
    n.lines.loc[lines_v_nom_b, 's_nom'] = (
        np.sqrt(3) * n.lines['type'].map(n.line_types.i_nom) *
        n.lines.bus0.map(n.buses.v_nom) * n.lines.num_parallel
    )

    # Replace transformers by lines
    trafo_map = pd.Series(n.transformers.bus1.values, index=n.transformers.bus0.values)
    trafo_map = trafo_map[~trafo_map.index.duplicated(keep='first')]
    several_trafo_b = trafo_map.isin(trafo_map.index)
    trafo_map.loc[several_trafo_b] = trafo_map.loc[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    trafo_map = trafo_map.append(pd.Series(missing_buses_i, missing_buses_i))

    for c in n.one_port_components|n.branch_components:
        df = n.df(c)
        for col in df.columns:
            if col.startswith('bus'):
                df[col] = df[col].map(trafo_map)

    n.mremove("Transformer", n.transformers.index)
    n.mremove("Bus", n.buses.index.difference(trafo_map))

    return n, trafo_map


def _prepare_connection_costs_per_link(n):
    if n.links.empty: return {}

    Nyears = n.snapshot_weightings.objective.sum() / 8760
    costs = load_costs(Nyears, snakemake.input.tech_costs,
                       snakemake.config['costs'], snakemake.config['electricity'])

    connection_costs_per_link = {}

    for tech in snakemake.config['renewable']:
        if tech.startswith('offwind'):
            connection_costs_per_link[tech] = (
                n.links.length * snakemake.config['lines']['length_factor'] *
                (n.links.underwater_fraction * costs.at[tech + '-connection-submarine', 'capital_cost'] +
                 (1. - n.links.underwater_fraction) * costs.at[tech + '-connection-underground', 'capital_cost'])
            )

    return connection_costs_per_link


def _compute_connection_costs_to_bus(n, busmap, connection_costs_per_link=None, buses=None):
    if connection_costs_per_link is None:
        connection_costs_per_link = _prepare_connection_costs_per_link(n)

    if buses is None:
        buses = busmap.index[busmap.index != busmap.values]

    connection_costs_to_bus = pd.DataFrame(index=buses)

    for tech in connection_costs_per_link:
        adj = n.adjacency_matrix(weights=pd.concat(dict(Link=connection_costs_per_link[tech].reindex(n.links.index),
                                                        Line=pd.Series(0., n.lines.index))))

        costs_between_buses = dijkstra(adj, directed=False, indices=n.buses.index.get_indexer(buses))
        connection_costs_to_bus[tech] = costs_between_buses[np.arange(len(buses)),
                                                            n.buses.index.get_indexer(busmap.loc[buses])]

    return connection_costs_to_bus


def _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus):
    connection_costs = {}
    for tech in connection_costs_to_bus:
        tech_b = n.generators.carrier == tech
        costs = n.generators.loc[tech_b, "bus"].map(connection_costs_to_bus[tech]).loc[lambda s: s>0]
        if not costs.empty:
            n.generators.loc[costs.index, "capital_cost"] += costs
            logger.info("Displacing {} generator(s) and adding connection costs to capital_costs: {} "
                        .format(tech, ", ".join("{:.0f} Eur/MW/a for `{}`".format(d, b) for b, d in costs.iteritems())))
            connection_costs[tech] = costs
    pd.DataFrame(connection_costs).to_csv(snakemake.output.connection_costs)



def _aggregate_and_move_components(n, busmap, connection_costs_to_bus, aggregate_one_ports={"Load", "StorageUnit"}):
    def replace_components(n, c, df, pnl):
        n.mremove(c, n.df(c).index)

        import_components_from_dataframe(n, df, c)
        for attr, df in pnl.items():
            if not df.empty:
                import_series_from_dataframe(n, df, c, attr)

    _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus)

    generators, generators_pnl = aggregategenerators(n, busmap, custom_strategies={'p_nom_min': np.sum})
    replace_components(n, "Generator", generators, generators_pnl)

    for one_port in aggregate_one_ports:
        df, pnl = aggregateoneport(n, busmap, component=one_port)
        replace_components(n, one_port, df, pnl)

    buses_to_del = n.buses.index.difference(busmap)
    n.mremove("Bus", buses_to_del)
    for c in n.branch_components:
        df = n.df(c)
        n.mremove(c, df.index[df.bus0.isin(buses_to_del) | df.bus1.isin(buses_to_del)])


def simplify_links(n):
    ## Complex multi-node links are folded into end-points
    logger.info("Simplifying connected link components")

    if n.links.empty:
        return n, n.buses.index.to_series()

    # Determine connected link components, ignore all links but DC
    adjacency_matrix = n.adjacency_matrix(branch_components=['Link'],
                                          weights=dict(Link=(n.links.carrier == 'DC').astype(float)))

    _, labels = connected_components(adjacency_matrix, directed=False)
    labels = pd.Series(labels, n.buses.index)

    G = n.graph()

    def split_links(nodes):
        nodes = frozenset(nodes)

        seen = set()
        supernodes = {m for m in nodes
                      if len(G.adj[m]) > 2 or (set(G.adj[m]) - nodes)}

        for u in supernodes:
            for m, ls in G.adj[u].items():
                if m not in nodes or m in seen: continue

                buses = [u, m]
                links = [list(ls)] #[name for name in ls]]

                while m not in (supernodes | seen):
                    seen.add(m)
                    for m2, ls in G.adj[m].items():
                        if m2 in seen or m2 == u: continue
                        buses.append(m2)
                        links.append(list(ls)) # [name for name in ls])
                        break
                    else:
                        # stub
                        break
                    m = m2
                if m != u:
                    yield pd.Index((u, m)), buses, links
            seen.add(u)

    busmap = n.buses.index.to_series()

    connection_costs_per_link = _prepare_connection_costs_per_link(n)
    connection_costs_to_bus = pd.DataFrame(0., index=n.buses.index, columns=list(connection_costs_per_link))

    for lbl in labels.value_counts().loc[lambda s: s > 2].index:

        for b, buses, links in split_links(labels.index[labels == lbl]):
            if len(buses) <= 2: continue

            logger.debug('nodes = {}'.format(labels.index[labels == lbl]))
            logger.debug('b = {}\nbuses = {}\nlinks = {}'.format(b, buses, links))

            m = sp.spatial.distance_matrix(n.buses.loc[b, ['x', 'y']],
                                           n.buses.loc[buses[1:-1], ['x', 'y']])
            busmap.loc[buses] = b[np.r_[0, m.argmin(axis=0), 1]]
            connection_costs_to_bus.loc[buses] += _compute_connection_costs_to_bus(n, busmap, connection_costs_per_link, buses)

            all_links = [i for _, i in sum(links, [])]

            p_max_pu = snakemake.config['links'].get('p_max_pu', 1.)
            lengths = n.links.loc[all_links, 'length']
            name = lengths.idxmax() + '+{}'.format(len(links) - 1)
            params = dict(
                carrier='DC',
                bus0=b[0], bus1=b[1],
                length=sum(n.links.loc[[i for _, i in l], 'length'].mean() for l in links),
                p_nom=min(n.links.loc[[i for _, i in l], 'p_nom'].sum() for l in links),
                underwater_fraction=sum(lengths/lengths.sum() * n.links.loc[all_links, 'underwater_fraction']),
                p_max_pu=p_max_pu,
                p_min_pu=-p_max_pu,
                underground=False,
                under_construction=False
            )

            logger.info("Joining the links {} connecting the buses {} to simple link {}".format(", ".join(all_links), ", ".join(buses), name))

            n.mremove("Link", all_links)

            static_attrs = n.components["Link"]["attrs"].loc[lambda df: df.static]
            for attr, default in static_attrs.default.iteritems(): params.setdefault(attr, default)
            n.links.loc[name] = pd.Series(params)

            # n.add("Link", **params)

    logger.debug("Collecting all components using the busmap")

    _aggregate_and_move_components(n, busmap, connection_costs_to_bus)
    return n, busmap

def remove_stubs(n):
    logger.info("Removing stubs")

    busmap = busmap_by_stubs(n) #  ['country'])

    connection_costs_to_bus = _compute_connection_costs_to_bus(n, busmap)

    _aggregate_and_move_components(n, busmap, connection_costs_to_bus)

    return n, busmap

def aggregate_to_substations(n, buses_i=None):
    # can be used to aggregate a selection of buses to electrically closest neighbors
    # if no buses are given, nodes that are no substations or without offshore connection are aggregated

    if buses_i is None:
        logger.info("Aggregating buses that are no substations or have no valid offshore connection")
        buses_i = list(set(n.buses.index)-set(n.generators.bus)-set(n.loads.bus))

    weight = pd.concat({'Line': n.lines.length/n.lines.s_nom.clip(1e-3),
                        'Link': n.links.length/n.links.p_nom.clip(1e-3)})

    adj = n.adjacency_matrix(branch_components=['Line', 'Link'], weights=weight)

    bus_indexer = n.buses.index.get_indexer(buses_i)
    dist = pd.DataFrame(dijkstra(adj, directed=False, indices=bus_indexer), buses_i, n.buses.index)

    dist[buses_i] = np.inf # bus in buses_i should not be assigned to different bus in buses_i

    for c in n.buses.country.unique():
        incountry_b = n.buses.country == c
        dist.loc[incountry_b, ~incountry_b] = np.inf

    busmap = n.buses.index.to_series()
    busmap.loc[buses_i] = dist.idxmin(1)

    clustering = get_clustering_from_busmap(n, busmap,
                                            bus_strategies=dict(country=_make_consense("Bus", "country")),
                                            aggregate_generators_weighted=True,
                                            aggregate_generators_carriers=None,
                                            aggregate_one_ports=["Load", "StorageUnit"],
                                            line_length_factor=1.0,
                                            generator_strategies={'p_nom_max': 'sum'},
                                            scale_link_capital_costs=False)

    return clustering.network, busmap


def cluster(n, n_clusters):
    logger.info(f"Clustering to {n_clusters} buses")

    focus_weights = snakemake.config.get('focus_weights', None)

    renewable_carriers = pd.Index([tech
                                    for tech in n.generators.carrier.unique()
                                    if tech.split('-', 2)[0] in snakemake.config['renewable']])
    def consense(x):
        v = x.iat[0]
        assert ((x == v).all() or x.isnull().all()), (
            "The `potential` configuration option must agree for all renewable carriers, for now!"
        )
        return v
    potential_mode = (consense(pd.Series([snakemake.config['renewable'][tech]['potential']
                                            for tech in renewable_carriers]))
                        if len(renewable_carriers) > 0 else 'conservative')
    clustering = clustering_for_n_clusters(n, n_clusters, custom_busmap=False, potential_mode=potential_mode,
                                           solver_name=snakemake.config['solving']['solver']['name'],
                                           focus_weights=focus_weights)

    return clustering.network, clustering.busmap


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('simplify_network', simpl='', network='elec')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    n, trafo_map = simplify_network_to_380(n)

    n, simplify_links_map = simplify_links(n)

    n, stub_map = remove_stubs(n)

    busmaps = [trafo_map, simplify_links_map, stub_map]

    if snakemake.config.get('clustering', {}).get('simplify', {}).get('to_substations', False):
        n, substation_map = aggregate_to_substations(n)
        busmaps.append(substation_map)

    if snakemake.wildcards.simpl:
        n, cluster_map = cluster(n, int(snakemake.wildcards.simpl))
        busmaps.append(cluster_map)

    # some entries in n.buses are not updated in previous functions, therefore can be wrong. as they are not needed
    # and are lost when clustering (for example with the simpl wildcard), we remove them for consistency:
    buses_c = {'symbol', 'tags', 'under_construction', 'substation_lv', 'substation_off'}.intersection(n.buses.columns)
    n.buses = n.buses.drop(buses_c, axis=1)

    update_p_nom_max(n)

    n.export_to_netcdf(snakemake.output.network)

    busmap_s = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    busmap_s.to_csv(snakemake.output.busmap)

    cluster_regions(busmaps, snakemake.input, snakemake.output)


## _helpers
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import pandas as pd
from pathlib import Path


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    import logging

    kwargs = snakemake.config.get('logging', dict())
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath('..', 'logs', f"{snakemake.rule}.log")
        logfile = snakemake.log.get('python', snakemake.log[0] if snakemake.log
                                    else fallback_path)
        kwargs.update(
            {'handlers': [
                # Prefer the 'python' log, otherwise take the first log for each
                # Snakemake rule
                logging.FileHandler(logfile),
                logging.StreamHandler()
                ]
            })
    logging.basicConfig(**kwargs)


def load_network(import_name=None, custom_components=None):
    """
    Helper for importing a pypsa.Network with additional custom components.

    Parameters
    ----------
    import_name : str
        As in pypsa.Network(import_name)
    custom_components : dict
        Dictionary listing custom components.
        For using ``snakemake.config['override_components']``
        in ``config.yaml`` define:

        .. code:: yaml

            override_components:
                ShadowPrice:
                    component: ["shadow_prices","Shadow price for a global constraint.",np.nan]
                    attributes:
                    name: ["string","n/a","n/a","Unique name","Input (required)"]
                    value: ["float","n/a",0.,"shadow value","Output"]

    Returns
    -------
    pypsa.Network
    """
    import pypsa
    from pypsa.descriptors import Dict

    override_components = None
    override_component_attrs = None

    if custom_components is not None:
        override_components = pypsa.components.components.copy()
        override_component_attrs = Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})
        for k, v in custom_components.items():
            override_components.loc[k] = v['component']
            override_component_attrs[k] = pd.DataFrame(columns = ["type","unit","default","description","status"])
            for attr, val in v['attributes'].items():
                override_component_attrs[k].loc[attr] = val

    return pypsa.Network(import_name=import_name,
                         override_components=override_components,
                         override_component_attrs=override_component_attrs)


def pdbcast(v, h):
    return pd.DataFrame(v.values.reshape((-1, 1)) * h.values,
                        index=v.index, columns=h.index)


def load_network_for_plots(fn, tech_costs, config, combine_hydro_ps=True):
    import pypsa
    from add_electricity import update_transmission_costs, load_costs

    n = pypsa.Network(fn)

    n.loads["carrier"] = n.loads.bus.map(n.buses.carrier) + " load"
    n.stores["carrier"] = n.stores.bus.map(n.buses.carrier)

    n.links["carrier"] = (n.links.bus0.map(n.buses.carrier) + "-" + n.links.bus1.map(n.buses.carrier))
    n.lines["carrier"] = "AC line"
    n.transformers["carrier"] = "AC transformer"

    n.lines['s_nom'] = n.lines['s_nom_min']
    n.links['p_nom'] = n.links['p_nom_min']

    if combine_hydro_ps:
        n.storage_units.loc[n.storage_units.carrier.isin({'PHS', 'hydro'}), 'carrier'] = 'hydro+PHS'

    # if the carrier was not set on the heat storage units
    # bus_carrier = n.storage_units.bus.map(n.buses.carrier)
    # n.storage_units.loc[bus_carrier == "heat","carrier"] = "water tanks"

    Nyears = n.snapshot_weightings.objective.sum() / 8760.
    costs = load_costs(Nyears, tech_costs, config['costs'], config['electricity'])
    update_transmission_costs(n, costs)

    return n

def update_p_nom_max(n):
    # if extendable carriers (solar/onwind/...) have capacity >= 0,
    # e.g. existing assets from the OPSD project are included to the network,
    # the installed capacity might exceed the expansion limit.
    # Hence, we update the assumptions.
    
    n.generators.p_nom_max = n.generators[['p_nom_min', 'p_nom_max']].max(1)

def aggregate_p_nom(n):
    return pd.concat([
        n.generators.groupby("carrier").p_nom_opt.sum(),
        n.storage_units.groupby("carrier").p_nom_opt.sum(),
        n.links.groupby("carrier").p_nom_opt.sum(),
        n.loads_t.p.groupby(n.loads.carrier,axis=1).sum().mean()
    ])

def aggregate_p(n):
    return pd.concat([
        n.generators_t.p.sum().groupby(n.generators.carrier).sum(),
        n.storage_units_t.p.sum().groupby(n.storage_units.carrier).sum(),
        n.stores_t.p.sum().groupby(n.stores.carrier).sum(),
        -n.loads_t.p.sum().groupby(n.loads.carrier).sum()
    ])

def aggregate_e_nom(n):
    return pd.concat([
        (n.storage_units["p_nom_opt"]*n.storage_units["max_hours"]).groupby(n.storage_units["carrier"]).sum(),
        n.stores["e_nom_opt"].groupby(n.stores.carrier).sum()
    ])

def aggregate_p_curtailed(n):
    return pd.concat([
        ((n.generators_t.p_max_pu.sum().multiply(n.generators.p_nom_opt) - n.generators_t.p.sum())
         .groupby(n.generators.carrier).sum()),
        ((n.storage_units_t.inflow.sum() - n.storage_units_t.p.sum())
         .groupby(n.storage_units.carrier).sum())
    ])

def aggregate_costs(n, flatten=False, opts=None, existing_only=False):

    components = dict(Link=("p_nom", "p0"),
                      Generator=("p_nom", "p"),
                      StorageUnit=("p_nom", "p"),
                      Store=("e_nom", "p"),
                      Line=("s_nom", None),
                      Transformer=("s_nom", None))

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components.keys(), skip_empty=False),
        components.values()
    ):
        if c.df.empty: continue
        if not existing_only: p_nom += "_opt"
        costs[(c.list_name, 'capital')] = (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == 'StorageUnit':
                p = p.loc[p > 0]
            costs[(c.list_name, 'marginal')] = (p*c.df.marginal_cost).groupby(c.df.carrier).sum()
    costs = pd.concat(costs)

    if flatten:
        assert opts is not None
        conv_techs = opts['conv_techs']

        costs = costs.reset_index(level=0, drop=True)
        costs = costs['capital'].add(
            costs['marginal'].rename({t: t + ' marginal' for t in conv_techs}),
            fill_value=0.
        )

    return costs

def progress_retrieve(url, file):
    import urllib
    from progressbar import ProgressBar

    pbar = ProgressBar(0, 100)

    def dlProgress(count, blockSize, totalSize):
        pbar.update( int(count * blockSize * 100 / totalSize) )

    urllib.request.urlretrieve(url, file, reporthook=dlProgress)


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import snakemake as sm
    import os
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert Path.cwd().resolve() == script_dir, \
      f'mock_snakemake has to be run from the repository scripts directory {script_dir}'
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    workflow = sm.Workflow(snakefile)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(job.input, job.output, job.params, job.wildcards,
                          job.threads, job.resources, job.log,
                          job.dag.workflow.config, job.rule.name, None,)
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake

## solve_network

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Solves linear optimal power flow for a network iteratively while updating reactances.

Relevant Settings
-----------------

.. code:: yaml

    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
            skip_iterations:
            track_iterations:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity_cf`, :ref:`solving_cf`, :ref:`plotting_cf`

Inputs
------

- ``networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: confer :ref:`prepare`

Outputs
-------

- ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Solved PyPSA network including optimisation results

    .. image:: ../img/results.png
        :scale: 40 %

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
The optimization is based on the ``pyomo=False`` setting in the :func:`network.lopf` and  :func:`pypsa.linopf.ilopf` function.
Additionally, some extra constraints specified in :mod:`prepare_network` are added.

Solving the network in multiple iterations is motivated through the dependence of transmission line capacities and impedances.
As lines are expanded their electrical parameters change, which renders the optimisation bilinear even if the power flow
equations are linearized.
To retain the computational advantage of continuous linear programming, a sequential linear programming technique
is used, where in between iterations the line impedances are updated.
Details (and errors made through this heuristic) are discussed in the paper

- Fabian Neumann and Tom Brown. `Heuristics for Transmission Expansion Planning in Low-Carbon Energy System Models <https://arxiv.org/abs/1907.10548>`_), *16th International Conference on the European Energy Market*, 2019. `arXiv:1907.10548 <https://arxiv.org/abs/1907.10548>`_.

.. warning::
    Capital costs of existing network components are not included in the objective function,
    since for the optimisation problem they are just a constant term (no influence on optimal result).

    Therefore, these capital costs are not included in ``network.objective``!

    If you want to calculate the full total annual system costs add these to the objective value.

.. tip::
    The rule :mod:`solve_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`solve_network`.

"""

import logging
from _helpers import configure_logging

import numpy as np
import pandas as pd
import re

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr, join_exprs,
                          network_lopf, ilopf)

from pathlib import Path
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)
from override_components import override_component_attrs #new


def prepare_network(n, solve_opts):

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        buses_i = n.buses.query("carrier == 'AC'").index
        n.madd("Generator", buses_i, " load",
               bus=buses_i,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
               )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components(n.one_port_components):
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                t.df['marginal_cost'] += (1e-2 + 2e-3 *
                                          (np.random.random(len(t.df)) - 0.5))

        for t in n.iterate_components(['Line', 'Link']):
            t.df['capital_cost'] += (1e-1 +
                2e-2*(np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760. / nhours

    return n


def add_CCL_constraints(n, config):
    agg_p_nom_limits = config['electricity'].get('agg_p_nom_limits')

    try:
        agg_p_nom_minmax = pd.read_csv(agg_p_nom_limits,
                                       index_col=list(range(2)))
    except IOError:
        logger.exception("Need to specify the path to a .csv file containing "
                          "aggregate capacity limits per country in "
                          "config['electricity']['agg_p_nom_limit'].")
    logger.info("Adding per carrier generation capacity constraints for "
                "individual countries")

    gen_country = n.generators.bus.map(n.buses.country)
    # cc means country and carrier
    p_nom_per_cc = (pd.DataFrame(
                    {'p_nom': linexpr((1, get_var(n, 'Generator', 'p_nom'))),
                    'country': gen_country, 'carrier': n.generators.carrier})
                    .dropna(subset=['p_nom'])
                    .groupby(['country', 'carrier']).p_nom
                    .apply(join_exprs))
    minimum = agg_p_nom_minmax['min'].dropna()
    if not minimum.empty:
        minconstraint = define_constraints(n, p_nom_per_cc[minimum.index],
                                           '>=', minimum, 'agg_p_nom', 'min')
    maximum = agg_p_nom_minmax['max'].dropna()
    if not maximum.empty:
        maxconstraint = define_constraints(n, p_nom_per_cc[maximum.index],
                                           '<=', maximum, 'agg_p_nom', 'max')


def add_EQ_constraints(n, o, scaling=1e-1):
    float_regex = "[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == 'c':
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = n.snapshot_weightings.generators @ \
           n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    inflow = n.snapshot_weightings.stores @ \
             n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    inflow = inflow.reindex(load.index).fillna(0.)
    rhs = scaling * ( level * load - inflow )
    lhs_gen = linexpr((n.snapshot_weightings.generators * scaling,
                       get_var(n, "Generator", "p").T)
              ).T.groupby(ggrouper, axis=1).apply(join_exprs)
    lhs_spill = linexpr((-n.snapshot_weightings.stores * scaling,
                         get_var(n, "StorageUnit", "spill").T)
                ).T.groupby(sgrouper, axis=1).apply(join_exprs)
    lhs_spill = lhs_spill.reindex(lhs_gen.index).fillna("")
    lhs = lhs_gen + lhs_spill
    define_constraints(n, lhs, ">=", rhs, "equity", "min")


def add_BAU_constraints(n, config):
    mincaps = pd.Series(config['electricity']['BAU_mincapacities'])
    lhs = (linexpr((1, get_var(n, 'Generator', 'p_nom')))
           .groupby(n.generators.carrier).apply(join_exprs))
    define_constraints(n, lhs, '>=', mincaps[lhs.index], 'Carrier', 'bau_mincaps')


def add_SAFE_constraints(n, config):
    peakdemand = (1. + config['electricity']['SAFE_reservemargin']) *\
                  n.loads_t.p_set.sum(axis=1).max()
    conv_techs = config['plotting']['conv_techs']
    exist_conv_caps = n.generators.query('~p_nom_extendable & carrier in @conv_techs')\
                       .p_nom.sum()
    ext_gens_i = n.generators.query('carrier in @conv_techs & p_nom_extendable').index
    lhs = linexpr((1, get_var(n, 'Generator', 'p_nom')[ext_gens_i])).sum()
    rhs = peakdemand - exist_conv_caps
    define_constraints(n, lhs, '>=', rhs, 'Safe', 'mintotalcap')


def add_battery_constraints(n):
    nodes = n.buses.index[n.buses.carrier == "battery"]
    if nodes.empty or ('Link', 'p_nom') not in n.variables.index:
        return
    link_p_nom = get_var(n, "Link", "p_nom")
    lhs = linexpr((1,link_p_nom[nodes + " charger"]),
                  (-n.links.loc[nodes + " discharger", "efficiency"].values,
                   link_p_nom[nodes + " discharger"].values))
    define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to ``pypsa.linopf.network_lopf``.
    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    if 'BAU' in opts and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if 'SAFE' in opts and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if 'CCL' in opts and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config)
    for o in opts:
        if "EQ" in o:
            add_EQ_constraints(n, o)
    add_battery_constraints(n)


def solve_network(n, config, opts='', **kwargs):
    solver_options = config['solving']['solver'].copy()
    solver_name = solver_options.pop('name')
    cf_solving = config['solving']['options']
    track_iterations = cf_solving.get('track_iterations', False)
    min_iterations = cf_solving.get('min_iterations', 4)
    max_iterations = cf_solving.get('max_iterations', 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get('skip_iterations', False):
        network_lopf(n, solver_name=solver_name, solver_options=solver_options,
                     extra_functionality=extra_functionality, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=extra_functionality, **kwargs)
    return n


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-BAU-CCL-24H')
    configure_logging(snakemake)

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split('-')
    solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:
        #n = pypsa.Network(snakemake.input[0])
        overrides = override_component_attrs(snakemake.input.overrides) #new
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides) #new
        n = prepare_network(n, solve_opts)
        n = solve_network(n, config=snakemake.config, opts=opts,
                          solver_dir=tmpdir,
                          solver_logfile=snakemake.log.solver)
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
