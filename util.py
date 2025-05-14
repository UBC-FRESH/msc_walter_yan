##################################################################################
# This module contain local utility function defintions that we can reuse 
# in example notebooks to help reduce clutter.
# #################################################################################

import ws3
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

##########################################################
# Implement a priority queue heuristic harvest scheduler
# #########################################################

def schedule_harvest_areacontrol(fm, period=None, acode='harvest', util_rate=0.85, 
                                 target_masks=None, target_areas=None,
                                 target_scalefactors=None,
                                 mask_area_thresh=0.,
                                 verbose=0):
    if not target_areas:
        if not target_masks: # default to AU-wise THLB 
            au_vals = []
            au_agg = []
            for au in fm.theme_basecodes(2):
                mask = '? 1 %s ? ?' % au
                masked_area = fm.inventory(0, mask=mask)
                if masked_area > mask_area_thresh:
                    au_vals.append(au)
                else:
                    au_agg.append(au)
                    if verbose > 0:
                        print('adding to au_agg', mask, masked_area)
            if au_agg:
                fm._themes[2]['areacontrol_au_agg'] = au_agg 
                if fm.inventory(0, mask='? ? areacontrol_au_agg ? ?') > mask_area_thresh:
                    au_vals.append('areacontrol_au_agg')
            target_masks = ['? 1 %s ? ?' % au for au in au_vals]
        target_areas = []
        for i, mask in enumerate(target_masks): # compute area-weighted mean CMAI age for each masked DT set
            masked_area = fm.inventory(0, mask=mask, verbose=verbose)
            if not masked_area: continue
            r = sum((fm.dtypes[dtk].ycomp('totvol').mai().ytp().lookup(0) * fm.dtypes[dtk].area(0)) for dtk in fm.unmask(mask))
            r /= masked_area
            asf = 1. if not target_scalefactors else target_scalefactors[i]  
            ta = (1/r) * fm.period_length * masked_area * asf
            target_areas.append(ta)
    periods = fm.periods if not period else [period]
    for period in periods:
        for mask, target_area in zip(target_masks, target_areas):
            if verbose > 0:
                print('calling areaselector', period, acode, target_area, mask)
            fm.areaselector.operate(period, acode, target_area, mask=mask, verbose=verbose)
    sch = fm.compile_schedule()
    return sch



##############################################################
# Implement an LP optimization harvest scheduler
# #############################################################

def cmp_c_z(fm, path, expr):
    """
    Compile objective function coefficient (given ForestModel instance, 
    leaf-to-root-node path, and expression to evaluate).
    """
    result = 0.
    for t, n in enumerate(path, start=1):
        d = n.data()
        if fm.is_harvest(d['acode']):
            result += fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result

def cmp_c_cflw(fm, path, expr, mask=None): # product, all harvest actions
    """
    Compile flow constraint coefficient for product indicator (given ForestModel 
    instance, leaf-to-root-node path, expression to evaluate, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if fm.is_harvest(d['acode']):
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result


def cmp_c_caa(fm, path, expr, acodes, mask=None): # product, named actions
    """
    Compile constraint coefficient for product indicator (given ForestModel 
    instance, leaf-to-root-node path, expression to evaluate, list of action codes, 
    and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if d['acode'] in acodes:
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result


def cmp_c_ci(fm, path, yname, mask=None): # product, named actions
    """
    Compile constraint coefficient for inventory indicator (given ForestModel instance, 
    leaf-to-root-node path, expression to evaluate, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['_dtk']): continue
        result[t] = fm.inventory(t, yname=yname, age=d['_age'], dtype_keys=[d['_dtk']]) 
        #result[t] = fm.inventory(t, yname=yname, age=d['age'], dtype_keys=[d['dtk']]) 
    return result


def compile_scenario(fm):
    oha = [fm.compile_product(period, '1.', acode='harvest') for period in fm.periods]
    ohv = [fm.compile_product(period, 'totvol * 0.85', acode='harvest') for period in fm.periods]
    ogs = [fm.inventory(period, 'totvol') for period in fm.periods]
    ocp = [fm.inventory(period, 'ecosystem') for period in fm.periods]
    ocf = [fm.inventory(period, 'total_emissions') for period in fm.periods]
    data = {'period':fm.periods, 
            'oha':oha, 
            'ohv':ohv, 
            'ogs':ogs,
            'ocp':ocp,
            'ocf':ocf}
    df = pd.DataFrame(data)
    return df


def plot_scenario(df):
    fig, ax = plt.subplots(1, 4, figsize=(20, 5))
    # Plot and label the first subplot for harvested area
    ax[0].bar(df.period, df.oha)
    ax[0].set_ylim(0, None)
    ax[0].set_title('Harvested area')
    ax[0].set_xlabel('Period')
    ax[0].set_ylabel('Area (ha)')
    
    # Plot and label the second subplot for harvested volume
    ax[1].bar(df.period, df.ohv)
    ax[1].set_ylim(0, None)
    ax[1].set_title('Harvested volume')
    ax[1].set_xlabel('Period')
    ax[1].set_ylabel('Volume (m3)')

    # Plot and label the third subplot for growing stock
    ax[2].bar(df.period, df.ogs)
    ax[2].set_ylim(0, None)
    ax[2].set_xlabel('Period')
    ax[2].set_title('Growing Stock')
    ax[2].set_ylabel('Volume (m3)')

    # Plot and label the fourth subplot for ecosystem carbon stock
    ax[3].bar(df.period, df.ocp)
    ax[3].set_ylim(0, None)
    ax[3].set_title('Ecosystem C stock')
    ax[3].set_xlabel('Period')
    ax[3].set_ylabel('Stock (ton)')

    return fig, ax


##############################################################
# Implement simple functions to run CBM from ws3 export data and output resutls
# #############################################################
def run_cbm(sit_config, sit_tables, n_steps):
    from libcbm.input.sit import sit_reader
    from libcbm.input.sit import sit_cbm_factory 
    from libcbm.model.cbm.cbm_output import CBMOutput
    from libcbm.storage.backends import BackendType
    from libcbm.model.cbm import cbm_simulator
    sit_data = sit_reader.parse(sit_classifiers=sit_tables['sit_classifiers'],
                                sit_disturbance_types=sit_tables['sit_disturbance_types'],
                                sit_age_classes=sit_tables['sit_age_classes'],
                                sit_inventory=sit_tables['sit_inventory'],
                                sit_yield=sit_tables['sit_yield'],
                                sit_events=sit_tables['sit_events'],
                                sit_transitions=sit_tables['sit_transitions'],
                                sit_eligibilities=None)
    sit = sit_cbm_factory.initialize_sit(sit_data=sit_data, config=sit_config)
    classifiers, inventory = sit_cbm_factory.initialize_inventory(sit)
    cbm_output = CBMOutput(classifier_map=sit.classifier_value_names,
                           disturbance_type_map=sit.disturbance_name_map)
    with sit_cbm_factory.initialize_cbm(sit) as cbm:
        # Create a function to apply rule based disturbance events and transition rules based on the SIT input
        rule_based_processor = sit_cbm_factory.create_sit_rule_based_processor(sit, cbm)
        # The following line of code spins up the CBM inventory and runs it through 200 timesteps.
        cbm_simulator.simulate(cbm,
                               n_steps=n_steps,
                               classifiers=classifiers,
                               inventory=inventory,
                               pre_dynamics_func=rule_based_processor.pre_dynamics_func,
                               reporting_func=cbm_output.append_simulation_result,
                               backend_type=BackendType.numpy)      
    return cbm_output


def cbm_report(fm, cbm_output, biomass_pools, dom_pools, total_emission, gross_growth, production, n_steps):
    
    # Add carbon pools indicators 
    pi = cbm_output.classifiers.to_pandas().merge(cbm_output.pools.to_pandas(), 
                                                  left_on=["identifier", "timestep"], 
                                                  right_on=["identifier", "timestep"])

    annual_carbon_stock = pd.DataFrame({'Year': pi['timestep'],
                                        'Biomass': pi[biomass_pools].sum(axis=1),
                                        'DOM': pi[dom_pools].sum(axis=1),
                                        'Ecosystem': pi[biomass_pools + dom_pools].sum(axis=1)})

    fi = cbm_output.classifiers.to_pandas().merge(cbm_output.flux.to_pandas(), 
                                                  left_on=["identifier", "timestep"], 
                                                  right_on=["identifier", "timestep"])

    annual_all_emission = pd.DataFrame({'Year': fi['timestep'],
                                        'All_Emissions': fi[total_emission].sum(axis=1)})

    # Calculating gross growth fluxes
    annual_gross_growth = pd.DataFrame({'Year': fi['timestep'],
                                        'Gross_Growth': fi[gross_growth].sum(axis=1)})
    
    # Calculating HWPs carbon
    annual_harvested_carbon = pd.DataFrame({'Year': fi['timestep'],
                                            'Harvested_Carbon': fi[production].sum(axis=1)})


    # Calculating net emissions (All_Emissions - Gross_Growth)
    annual_net_emission = pd.DataFrame({'Year': fi['timestep'],
                                        'Net_Emissions': annual_all_emission['All_Emissions'] - annual_gross_growth['Gross_Growth']})

    df_cs = annual_carbon_stock.groupby('Year').sum()
    df_ae = annual_all_emission.groupby('Year').sum()
    df_gg = annual_gross_growth.groupby('Year').sum()  # Updated to Gross Growth
    df_ne = annual_net_emission.groupby('Year').sum()
    df_hc = annual_harvested_carbon.groupby('Year').sum()

    # Merging all dataframes, now including net emissions and gross growth
    merged_df = pd.merge(pd.merge(pd.merge(pd.merge(df_cs, df_ae, left_index=True, right_index=True, how='outer'),
                                           df_gg, left_index=True, right_index=True, how='outer'),
                                  df_ne, left_index=True, right_index=True, how='outer'),
                         df_hc, left_index=True, right_index=True, how='outer')

    # Calculating stock change
    merged_df['Stock_Change'] = merged_df['Ecosystem'].diff() * (-1)
    merged_df.at[0, 'Stock_Change'] = 0

    # Calculating (Stock Change - Harvested Carbon)
    merged_df['Stock_Change_minus_Harvested_Carbon'] = merged_df['Stock_Change'] - merged_df['Harvested_Carbon']

    # Plotting the graphs
    fig, axs = plt.subplots(8, 1, figsize=(10, 40))

    # Plot 1: Biomass Stock
    axs[0].plot(merged_df.index, merged_df['Biomass'], label='Biomass Stock', color='green')
    axs[0].set_title("Annual Biomass Stock")
    axs[0].set_xlabel("Year")
    axs[0].set_ylabel("Stock (tC)")
    axs[0].set_xlim(0, n_steps)

    # Plot 2: DOM Stock
    axs[1].plot(merged_df.index, merged_df['DOM'], label='DOM Stock', color='brown')
    axs[1].set_title("Annual DOM Stock")
    axs[1].set_xlabel("Year")
    axs[1].set_ylabel("Stock (tC)")
    axs[1].set_xlim(0, n_steps)

    # Plot 3: Ecosystem Stock
    axs[2].plot(merged_df.index, merged_df['Ecosystem'], label='Ecosystem Stock', color='blue')
    axs[2].set_title("Annual Ecosystem Stock")
    axs[2].set_xlabel("Year")
    axs[2].set_ylabel("Stock (tC)")
    axs[2].set_xlim(0, n_steps)

    # Plot 4: All Emissions
    axs[3].plot(merged_df.index, merged_df['All_Emissions'], label='All Emissions', color='red')
    axs[3].set_title("Annual Total Ecosystem Carbon Emissions")
    axs[3].set_xlabel("Year")
    axs[3].set_ylabel("Emissions (tC)")
    axs[3].set_xlim(0, n_steps)

    # Plot 5: Gross Growth
    axs[4].plot(merged_df.index, merged_df['Gross_Growth'], label='Gross Growth', color='purple')
    axs[4].set_title("Annual Gross Growth")
    axs[4].set_xlabel("Year")
    axs[4].set_ylabel("Growth (tC)")
    axs[4].set_xlim(0, n_steps)

    # Plot 6: Net Emissions
    axs[5].plot(merged_df.index, merged_df['Net_Emissions'], label='Net Emissions', color='orange')
    axs[5].set_title("Annual Net Ecosystem Carbon Emissions")
    axs[5].set_xlabel("Year")
    axs[5].set_ylabel("Emissions (tC)")
    axs[5].set_xlim(0, n_steps)

    # Plot 7: Stock Change
    axs[6].plot(merged_df.index, merged_df['Stock_Change'], label='Stock Change', color='cyan')
    axs[6].set_title("Annual Ecosystem Carbon Stock Change")
    axs[6].set_xlabel("Year")
    axs[6].set_ylabel("Stock change (tC)")
    axs[6].set_xlim(0, n_steps)

    # Plot 8: Harvested Carbon
    axs[7].plot(merged_df.index, merged_df['Harvested_Carbon'], label='Harvested Carbon', color='magenta')
    axs[7].set_title("Annual Harvested Carbon")
    axs[7].set_xlabel("Year")
    axs[7].set_ylabel("Harvested stock (tC)")
    axs[7].set_xlim(0, n_steps)

    # Rest of the code remains the same

def compare_ws3_libcbm(
    fm, cbm_output, disturbance_type_mapping, biomass_pools, dom_pools,
    total_emission, gross_growth, plots, filename=None
):
    eco_pools = biomass_pools + dom_pools
    pi = cbm_output.classifiers.to_pandas().merge(
        cbm_output.pools.to_pandas(), 
        left_on=["identifier", "timestep"], 
        right_on=["identifier", "timestep"]
    )
    fi = cbm_output.classifiers.to_pandas().merge(
        cbm_output.flux.to_pandas(),
        left_on=["identifier", "timestep"], 
        right_on=["identifier", "timestep"]
    )

    # Summarize CBM data
    df_cbm = pd.DataFrame({
        'period': pi["timestep"] * 0.1, 
        'biomass_stock': pi[biomass_pools].sum(axis=1),
        'dom_stock': pi[dom_pools].sum(axis=1),
        'eco_stock': pi[eco_pools].sum(axis=1),
        'total_emission': fi[total_emission].sum(axis=1),
        'gross_growth': fi[gross_growth].sum(axis=1)
    }).groupby('period').sum().iloc[1::10, :].reset_index()

    df_cbm['period'] = (df_cbm['period'] + 0.9).astype(int)
    df_cbm['net_emission'] = df_cbm['total_emission'] - df_cbm['gross_growth']

    # Summarize WS3 data
    df_ws3 = pd.DataFrame({
        'period': fm.periods,
        'biomass_stock': [fm.inventory(period, 'biomass') for period in fm.periods],
        'dom_stock': [fm.inventory(period, 'DOM') for period in fm.periods],
        'eco_stock': [fm.inventory(period, 'ecosystem') for period in fm.periods],
        'net_emission': [fm.inventory(period, 'net_emission') for period in fm.periods]
    })

    # Plot
    if plots == "whole":
        # Create a single figure
        fig = plt.figure(figsize=(10, 6))

        # Ecosystem stock
        plt.plot(df_cbm['period'], df_cbm['eco_stock'], label='libcbm Ecosystem Stock')
        plt.plot(df_ws3['period'], df_ws3['eco_stock'], label='ws3 Ecosystem Stock')

        # Biomass stock
        plt.plot(df_cbm['period'], df_cbm['biomass_stock'], label='libcbm Biomass Stock')
        plt.plot(df_ws3['period'], df_ws3['biomass_stock'], label='ws3 Biomass Stock')

        # DOM stock
        plt.plot(df_cbm['period'], df_cbm['dom_stock'], label='libcbm DOM Stock')
        plt.plot(df_ws3['period'], df_ws3['dom_stock'], label='ws3 DOM Stock')

        plt.xlabel('Period')
        plt.ylabel('Stock (tC)')
        plt.ylim(0, None)
        ticks = np.arange(df_cbm['period'].min() - 1, df_cbm['period'].max() + 1, 2)
        plt.xticks(ticks)
        plt.legend()
        plt.tight_layout()

    elif plots == "individual":
        # Create a figure with subplots
        fig, axs = plt.subplots(4, 1, figsize=(10, 16))
        ticks = np.arange(df_cbm['period'].min() - 1, df_cbm['period'].max() + 1, 2)

        # 1) Ecosystem stock
        axs[0].plot(df_cbm['period'], df_cbm['eco_stock'], label='libcbm Ecosystem Stock')
        axs[0].plot(df_ws3['period'], df_ws3['eco_stock'], label='ws3 Ecosystem Stock')
        axs[0].set_xlabel('Period')
        axs[0].set_ylabel('Stock (tC)')
        axs[0].set_xticks(ticks)
        axs[0].legend()

        # 2) Biomass stock
        axs[1].plot(df_cbm['period'], df_cbm['biomass_stock'], label='libcbm Biomass Stock')
        axs[1].plot(df_ws3['period'], df_ws3['biomass_stock'], label='ws3 Biomass Stock')
        axs[1].set_xlabel('Period')
        axs[1].set_ylabel('Stock (tC)')
        axs[1].set_xticks(ticks)
        axs[1].legend()

        # 3) DOM stock
        axs[2].plot(df_cbm['period'], df_cbm['dom_stock'], label='libcbm DOM Stock')
        axs[2].plot(df_ws3['period'], df_ws3['dom_stock'], label='ws3 DOM Stock')
        axs[2].set_xlabel('Period')
        axs[2].set_ylabel('Stock (tC)')
        axs[2].set_xticks(ticks)
        axs[2].legend()


        # 5) Net ecosystem carbon emission
        axs[3].plot(df_cbm['period'], df_cbm['net_emission'], label='libcbm Net Ecosystem Emission')
        axs[3].plot(df_ws3['period'], df_ws3['net_emission'], label='ws3 Net Ecosystem Emission')
        axs[3].set_xlabel('Period')
        axs[3].set_ylabel('Emission (tC)')
        axs[3].set_xticks(ticks)
        axs[3].legend(loc='lower right')

        plt.tight_layout()

    # If a filename is provided, save the figure to PDF
    if filename is not None:
        plt.savefig(filename, format='pdf', bbox_inches='tight')

    # Finally, show the figure
    plt.show()

    return df_cbm, df_ws3   


##############################################################
# Implement simple functions to generate, plug-in, and fix carbon yield curves into ws3 models
# #############################################################

def generate_c_curves(fm, disturbance_type_mapping, pools, fluxes):
    for dtype_key in fm.dtypes:
        fm.dt(dtype_key).last_pass_disturbance = 'fire' if dtype_key[2] == dtype_key[4] else 'harvest'
    sit_config, sit_tables = fm.to_cbm_sit(softwood_volume_yname='swdvol', 
                                           hardwood_volume_yname='hwdvol', 
                                           admin_boundary='British Columbia', 
                                           eco_boundary='Montane Cordillera',
                                           disturbance_type_mapping=disturbance_type_mapping)
    
    df = sit_tables['sit_inventory']
    df = df.iloc[0:0]
    data = []
    for dtype_key in fm.dtypes:
        dt = fm.dt(dtype_key)
        values = list(dtype_key) 
        values += [dt.leading_species, 'FALSE', 0, 1., 0, 0, 'fire', 'fire' if dtype_key[2] == dtype_key[4] else 'harvest']
        data.append(dict(zip(df.columns, values)))
    sit_tables['sit_inventory'] = pd.DataFrame(data)
    
    n_steps = fm.horizon * fm.period_length
    cbm_output = run_cbm(sit_config, sit_tables, n_steps)
    
    pi = cbm_output.classifiers.to_pandas().merge(cbm_output.pools.to_pandas(), 
                                                  left_on=["identifier", "timestep"], 
                                                  right_on=["identifier", "timestep"])
    fi = cbm_output.classifiers.to_pandas().merge(cbm_output.flux.to_pandas(), 
                                                  left_on=["identifier", "timestep"], 
                                                  right_on=["identifier", "timestep"])

    pi['dtype_key'] = pi.apply(lambda r: '%s %s %s %s %s' % (r['theme0'], r['theme1'], r['theme2'], r['theme3'], r['theme4']), axis=1)
    fi['dtype_key'] = fi.apply(lambda r: '%s %s %s %s %s' % (r['theme0'], r['theme1'], r['theme2'], r['theme3'], r['theme4']), axis=1)

    c_curves_p = pi.groupby(['dtype_key', 'timestep'], as_index=True)[pools].sum()
    c_curves_f = fi.groupby(['dtype_key', 'timestep'], as_index=True)[fluxes].sum()
    
    #c_curves_p[pools] = c_curves_p[pools].apply(lambda x: x*(1+pool_corr))
    #c_curves_f[fluxes] = c_curves_f[fluxes].apply(lambda x: x*flux_corr)
    
    return c_curves_p, c_curves_f

def plugin_c_curves(fm, c_curves_p, c_curves_f, pools, fluxes):
    # Dictionary to track registered curves for each dtype_key
    registered_curves = {}

    for dtype_key in fm.dtypes:
        dt = fm.dt(dtype_key)
        mask = ('?', '?', dtype_key[2], '?', dtype_key[4])
        
        for _mask, ytype, curves in fm.yields:
            if _mask != mask: 
                continue  # Only proceed if the mask matches

            print('found match for mask', mask)

            # Initialize the tracking of registered curves for the dtype_key if not already done
            if dtype_key not in registered_curves:
                registered_curves[dtype_key] = set()

            # Register pool curves
            pool_data = c_curves_p.loc[' '.join(dtype_key)]
            for yname in pools:
                if yname not in registered_curves[dtype_key]:  # Check if curve is already registered
                    points = list(zip(pool_data.index.values, pool_data[yname]))
                    curve = fm.register_curve(ws3.core.Curve(yname, 
                                                             points=points, 
                                                             type='a', 
                                                             is_volume=False,
                                                             xmax=fm.max_age,
                                                             period_length=fm.period_length))
                    curves.append((yname, curve))
                    dt.add_ycomp('a', yname, curve)

                    # Mark the curve as registered
                    registered_curves[dtype_key].add(yname)

            # Register flux curves
            flux_data = c_curves_f.loc[' '.join(dtype_key)]
            for yname in fluxes:
                if yname not in registered_curves[dtype_key]:  # Check if curve is already registered
                    points = list(zip(flux_data.index.values, flux_data[yname]))
                    curve = fm.register_curve(ws3.core.Curve(yname, 
                                                             points=points, 
                                                             type='a', 
                                                             is_volume=False,
                                                             xmax=fm.max_age,
                                                             period_length=fm.period_length))
                    curves.append((yname, curve))
                    dt.add_ycomp('a', yname, curve)

                    # Mark the curve as registered
                    registered_curves[dtype_key].add(yname)

               
def compile_events(self, softwood_volume_yname, hardwood_volume_yname, n_yield_vals):
    
    def leading_species(dt):
        """
        Determine if softwood or hardwood leading species by comparing softwood and hardwood
        volume at peak MAI age.
        """
        svol_curve, hvol_curve = dt.ycomp(softwood_volume_yname), dt.ycomp(hardwood_volume_yname)
        tvol_curve = svol_curve + hvol_curve
        x_cmai = tvol_curve.mai().ytp().lookup(0)
        return 'softwood' if svol_curve[x_cmai] > hvol_curve[x_cmai] else 'hardwood'

    for dtype_key in self.dtypes:
        dt = self.dt(dtype_key)
        dt.leading_species = leading_species(dt)
    
    theme_cols = [theme['__name__'] for theme in self._themes]
    columns = theme_cols.copy()
    columns += ['species',
                'using_age_class',
                'min_softwood_age',
                'max_softwood_age',
                'min_hardwood_age',
                'max_hardwood_age',
                'MinYearsSinceDist',
                'MaxYearsSinceDist',
                'LastDistTypeID',
                'MinTotBiomassC',
                'MaxTotBiomassC',
                'MinSWMerchBiomassC',
                'MaxSWMerchBiomassC',
                'MinHWMerchBiomassC',
                'MaxHWMerchBiomassC',
                'MinTotalStemSnagC',
                'MaxTotalStemSnagC',	
                'MinSWStemSnagC',
                'MaxSWStemSnagC',
                'MinHWStemSnagC',
                'MaxHWStemSnagC',
                'MinTotalStemSnagMerchC',
                'MaxTotalStemSnagMerchC',
                'MinSWMerchStemSnagC',
                'MaxSWMerchStemSnagC',
                'MinHWMerchStemSnagC',
                'MaxHWMerchStemSnagC',
                'efficiency',
                'sort_type',
                'target_type',
                'target',
                'disturbance_type',
                'disturbance_year']
    data = {c:[] for c in columns}
    for dtype_key, age, area, acode, period, _ in self.compile_schedule():
        #set_trace()
        for i, c in enumerate(theme_cols): data[c].append(dtype_key[i])
        data['species'].append(self.dt(dtype_key).leading_species)
        data['using_age_class'].append('FALSE')
        #############################################################################
        # might need to be more flexible with age range, to avoid OBO bugs and such?
        data['min_softwood_age'].append(-1)
        data['max_softwood_age'].append(-1)
        data['min_hardwood_age'].append(-1)
        data['max_hardwood_age'].append(-1)
        # data['min_softwood_age'].append(age)
        # data['max_softwood_age'].append(age)
        # data['min_hardwood_age'].append(age)
        # data['max_hardwood_age'].append(age)
        #############################################################################
        for c in columns[11:-6]: data[c].append(-1)
        data['efficiency'].append(1)
        data['sort_type'].append(3) # oldest first (see Table 3-3 in the CBM-CFS3 user guide)
        data['target_type'].append('A') # area target
        data['target'].append(area)
        data['disturbance_type'].append(acode)
        data['disturbance_year'].append((period-1)*self.period_length+1)
    sit_events = pd.DataFrame(data)         
    return sit_events

def cmp_c_cs(fm, path, expr, yname, half_life_solid_wood, half_life_paper, proportion_solid_wood, mask=None):
    """
    Compile objective function coefficient for total system carbon stock indicators (given ForestModel instance, 
    leaf-to-root-node path, and expression to evaluate).
    """
    
    result = 0.
    
    # Calculate decay rates based on half-lives
    k_solid_wood = math.log(2) / half_life_solid_wood  # Decay rate for solid wood products (30-year half-life)
    k_paper = math.log(2) / half_life_paper  # Decay rate for paper (2-year half-life)
    
    # Define the allocation distribution
    proportion_paper = 1 - proportion_solid_wood
    
    # wood density (Kennedy, 1965)
    wood_density = 460

    # carbon content
    carbon_content = 0.5
    
    product_stock_dict = {}  # Dictionary to track product stock for each node across iterations
    
    for t, n in enumerate(path, start=1):

        d = n.data()
        node_id = id(n)  # or another unique identifier specific to your application
        
        # Track the ecosystem carbon stock
        if mask and not fm.match_mask(mask, d['_dtk']): continue
        result = fm.inventory(t, yname, age=d['_age'], dtype_keys=[d['_dtk']])
        
        # Retrieve the last tuple of stocks from the dictionary
        last_stocks = next(reversed(product_stock_dict.values()), (0, 0))
        old_product_stock_solid_wood, old_product_stock_paper = last_stocks
        
        if fm.is_harvest(d['acode']):
            # Calculate new product stock
            new_product_carbon = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False) * wood_density * carbon_content / 1000 # Convert kg to ton
            new_product_stock_solid_wood = new_product_carbon * proportion_solid_wood
            new_product_stock_paper = new_product_carbon * proportion_paper 

            # Apply decay to old stocks and add new stocks
            # Apply decay to all stocks within the same period they're produced
            sum_product_stock_solid_wood = old_product_stock_solid_wood * (1 - k_solid_wood)**10 + new_product_stock_solid_wood
            sum_product_stock_paper = (old_product_stock_paper) * (1 - k_paper)**10 + new_product_stock_paper
        
        else:
            # If not harvesting, simply apply decay to the old product stocks
            sum_product_stock_solid_wood = old_product_stock_solid_wood * (1 - k_solid_wood)
            sum_product_stock_paper = old_product_stock_paper * (1 - k_paper)
            
        # Update product_stock_dict with the new sum product stocks for this node
        product_stock_dict[node_id] = (sum_product_stock_solid_wood, sum_product_stock_paper)

        ecosystem_stock = fm.inventory(t, yname, age=d['_age'], dtype_keys=[d['_dtk']])
        result = ecosystem_stock + sum_product_stock_solid_wood + sum_product_stock_paper
        
    return result

def cmp_c_ce(fm, path, expr, half_life_solid_wood, half_life_paper, proportion_solid_wood, util_rate, displacement_factor_solid, credit, mask=None):
        
        """
        Compile objective function coefficient for net system carbon emission indicators 
        (given ForestModel instance, leaf-to-root-node path, and expression to evaluate).
        """
        
        hwps_solid_pool = 0.
        hwps_paper_pool = 0.
        hwps_residue_pool = 0.
        hwps_solid_emission = 0.
        hwps_paper_emission = 0.
        hwps_residue_emission = 0.
        result = 0.
        
        # Calculate decay rates based on half-lives
        k_solid_wood = math.log(2) / half_life_solid_wood  # Decay rate for solid wood products (30-year half-life)
        k_paper = math.log(2) / half_life_paper  # Decay rate for paper (2-year half-life)
        
        # Define the allocation distribution
        proportion_paper = 1 - proportion_solid_wood
        
        # wood density (Kennedy, 1965)
        wood_density = 460
    
        # carbon content
        carbon_content = 0.5
        
        for t, n in enumerate(path, start=1):
            
            d = n.data()
    
            eco_emission = fm.inventory(t, 'net_emission', age=d['_age'], dtype_keys=[d['_dtk']])
            
            # Track the new product stock
            if fm.is_harvest(d['acode']):             
                # Calculate new product stock
                new_product_stock = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False) * wood_density * carbon_content / 1000  # Convert kg to ton
            else:
                new_product_stock = 0.
            
            # Calculate product stock
            hwps_solid_pool = hwps_solid_pool * (1 - k_solid_wood)**10 + new_product_stock * util_rate * proportion_solid_wood
            hwps_paper_pool = hwps_paper_pool * (1 - k_paper)**10 + new_product_stock * util_rate * proportion_paper
            hwps_residue_pool = new_product_stock * (1.0 - util_rate)
    
            # Calculate prodcut emission
            hwps_solid_emission = hwps_solid_pool * (1- (1 - k_solid_wood)**10)
            hwps_paper_emission = hwps_paper_pool * (1- (1 - k_paper)**10)
            hwps_residue_emission = hwps_residue_pool
            hwps_sum_emission = hwps_solid_emission + hwps_paper_emission + hwps_residue_emission
           
            # Calculate Substitution Effect
            substitution_effect_solid = new_product_stock * util_rate * proportion_solid_wood * displacement_factor_solid * credit  # Emissions avoided by using HWPs
            substitution_effect_paper = 0.
            substitution_effect_residue = 0.
            substitution_effect_sum = substitution_effect_solid + substitution_effect_paper + substitution_effect_residue
            
            # Accumlate the total system carbon stock in each timestep
            result += (eco_emission + hwps_sum_emission - substitution_effect_sum) * -1 * 44 / 12
            # result += (eco_emission + hwps_sum_emission - substitution_effect_sum)
        
        return result

def track_system_stock(fm, half_life_solid_wood, half_life_paper, proportion_solid_wood, util_rate):
    
    hwps_solid_pool = 0.
    hwps_paper_pool = 0.
    hwps_residue_pool = 0.

    sys_pool_list = []
    eco_pool_list = []
    hwps_solid_pool_list = []
    hwps_paper_pool_list = []
    hwps_residue_pool_list = []
    hwps_sum_pool_list = []
    
    # Calculate decay rates based on half-lives
    k_solid_wood = math.log(2) / half_life_solid_wood
    k_paper = math.log(2) / half_life_paper

    # Define the allocation distribution
    proportion_paper = 1 - proportion_solid_wood

    # Constants
    wood_density = 460 #(Kennedy, 1965)
    carbon_content = 0.5

    for period in fm.periods:

        # Track Ecosystem Carbon Stock
        eco_pool = fm.inventory(period, 'ecosystem')
        
        # Calculate new product stocks
        new_product_stock = fm.compile_product(period, 'totvol', acode='harvest') * wood_density * carbon_content / 1000 # Convert kg to ton
        # Calculate product stock
        hwps_solid_pool = hwps_solid_pool * (1 - k_solid_wood)**10 + new_product_stock * util_rate * proportion_solid_wood
        hwps_paper_pool = hwps_paper_pool * (1 - k_paper)**10 + new_product_stock * util_rate * proportion_paper
        hwps_residue_pool = new_product_stock * (1.0 - util_rate)
        hwps_sum_pool = hwps_solid_pool + hwps_paper_pool

        # Calculate total system carbon stock
        sys_pool = eco_pool + hwps_sum_pool

        # Update stock lists for this period
        sys_pool_list.append(sys_pool)
        eco_pool_list.append(eco_pool)
        hwps_solid_pool_list.append(hwps_solid_pool)
        hwps_paper_pool_list.append(hwps_paper_pool)
        hwps_residue_pool_list.append(hwps_residue_pool)
        hwps_sum_pool_list.append(hwps_sum_pool)

    # Prepare data for plotting
    data = {
        'period': fm.periods,
        'solid_wood': hwps_solid_pool_list,
        'paper': hwps_paper_pool_list,
        'residue': hwps_residue_pool_list,
        'sum_product': hwps_sum_pool_list,
        'ecosystem': eco_pool_list,
        'system': sys_pool_list
    }

    df = pd.DataFrame(data)

    return df

def track_system_emission(fm, half_life_solid_wood, half_life_paper, proportion_solid_wood, util_rate, displacement_factor_solid, credit):
    
    hwps_solid_pool = 0.
    hwps_paper_pool = 0.
    hwps_residue_pool = 0.
    hwps_solid_emission = 0.
    hwps_paper_emission = 0.
    hwps_residue_emission = 0.

    sys_emission_list = []
    eco_emission_list = []
    hwps_solid_emission_list = []
    hwps_paper_emission_list = []
    hwps_residue_emission_list = []
    hwps_sum_emission_list = []
    
    # Calculate decay rates based on half-lives
    k_solid_wood = math.log(2) / half_life_solid_wood
    k_paper = math.log(2) / half_life_paper

    # Define the allocation distribution
    proportion_paper = 1-proportion_solid_wood

    # Constants
    wood_density = 460 #(Kennedy, 1965)
    carbon_content = 0.5

    for period in fm.periods:

        # Calculate ecosytem emission
        eco_emission = fm.inventory(period, 'net_emission') * 44 / 12

        # Calculate new product stocks
        new_product_stock = fm.compile_product(period, 'totvol', acode='harvest') * wood_density * carbon_content / 1000

        # Calculate total product stock
        hwps_solid_pool = hwps_solid_pool * (1 - k_solid_wood)**10 + new_product_stock * util_rate * proportion_solid_wood
        hwps_paper_pool = hwps_paper_pool * (1 - k_paper)**10 + new_product_stock * util_rate * proportion_paper
        hwps_residue_pool =  new_product_stock * (1- util_rate)

        # Calculate product emission
        hwps_solid_emission = hwps_solid_pool*(1- (1 - k_solid_wood)**10) * 44 / 12 #Convert from C to CO2
        hwps_paper_emission = hwps_paper_pool*(1- (1 - k_paper)**10) * 44 / 12 #Convert from C to CO2
        hwps_residue_emission = hwps_residue_pool * 44 / 12 #Convert from C to CO2
        
        hwps_sum_emission = hwps_solid_emission + hwps_paper_emission + hwps_residue_emission

        # Calculate Substitution Effect
        substitution_effect_solid = new_product_stock * util_rate * proportion_solid_wood * displacement_factor_solid * credit # Emissions avoided by using HWPs
        substitution_effect_paper = 0.
        substitution_effect_residue = 0.
        substitution_effect_sum = (substitution_effect_solid + substitution_effect_paper + substitution_effect_residue) * 44 / 12 #Convert from C to CO2
        
        # Calculate net system carbon emission
        sys_emission = (eco_emission + hwps_sum_emission - substitution_effect_sum)

        # Update stock lists for this period
        sys_emission_list.append(sys_emission)
        eco_emission_list.append(eco_emission)
        hwps_solid_emission_list.append(hwps_solid_emission)
        hwps_paper_emission_list.append(hwps_paper_emission)
        hwps_residue_emission_list.append(hwps_residue_emission)
        hwps_sum_emission_list.append(hwps_sum_emission)

    # Prepare data for plotting
    data = {
        'period': fm.periods,
        'solid_wood': hwps_solid_emission_list,
        'paper': hwps_paper_emission_list,
        'residue': hwps_residue_emission_list,
        'sum_product': hwps_sum_emission_list,
        'ecosystem': eco_emission_list,
        'system': sys_emission_list
    }

    df = pd.DataFrame(data)
    
    # Plotting
    fig, ax = plt.subplots(1, 6, figsize=(20, 4))  # Adjusted for 5 subplots
    ax[0].bar(df.period, df.solid_wood)
    ax[0].set_title('Solid Wood Product CO2 Emission')
    ax[1].bar(df.period, df.paper)
    ax[1].set_title('Paper Product CO2 Emission')
    ax[2].bar(df.period, df.residue)
    ax[2].set_title('Residue CO2 Emission')
    ax[3].bar(df.period, df.sum_product)
    ax[3].set_title('Total Product CO2 Emission')
    ax[4].bar(df.period, df.ecosystem)
    ax[4].set_title('Net Ecosystem CO2 Emission')
    ax[5].bar(df.period, df.system)
    ax[5].set_title('Net System CO2 Emission')

    for a in ax:
        a.set_ylim(None, None)
        a.set_xlabel('Period')
        a.set_ylabel('CO2 (tons)')

    plt.tight_layout()

    return df

def gen_scenario(fm, name='base', util_rate=0.85, harvest_acode='harvest',
                 cflw_ha={}, cflw_hv={}, 
                 cgen_ha={}, cgen_hv={}, 
                 cgen_gs={}, tvy_name='totvol', cp_name='ecosystem', cf_name='total_emissions', obj_mode='max', mask=None):
    
    from functools import partial
    import numpy as np
    coeff_funcs = {}
    cflw_e = {}
    cgen_data = {}
    acodes = ['null', harvest_acode] # define list of action codes
    vexpr = '%s * %0.2f' % (tvy_name, util_rate) # define volume expression

    #Define constants from product carbon estimation

    if obj_mode == 'max': # maximize harvest volume
        sense = ws3.opt.SENSE_MAXIMIZE 
        zexpr = vexpr
    elif obj_mode == 'min': # maximize harvest volume
        sense = ws3.opt.SENSE_MINIMIZE 
        zexpr = vexpr
    else:
        raise ValueError('Invalid obj_mode: %s' % obj_mode)
    
    # coeff_funcs['z'] = partial(cmp_c_i, yname=cf_name) # define objective function coefficient function for inventory data
    # coeff_funcs['z'] = partial(cmp_c_id, yname=cf_name) # define objective function coefficient function for inventory change data
    coeff_funcs['z'] = partial(cmp_c_z, expr=vexpr) # define objective function coefficient function for harvest volume
    # coeff_funcs['z'] = partial(cmp_c_cs, expr=vexpr) # define objective function coefficient function for total system carbon stock
    # coeff_funcs['z'] = partial(cmp_c_ce, expr=vexpr) # define objective function coefficient function for net system carbon emission
    
    T = fm.periods
    if cflw_ha: # define even flow constraint (on harvest area)
        cname = 'cflw_ha'
        coeff_funcs[cname] = partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=None) 
        cflw_e[cname] = cflw_ha
    if cflw_hv: # define even flow constraint (on harvest volume)
        cname = 'cflw_hv'
        coeff_funcs[cname] = partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=None) 
        cflw_e[cname] = cflw_hv         
    if cgen_ha: # define general constraint (harvest area)
        cname = 'cgen_ha'
        coeff_funcs[cname] = partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=None) 
        cgen_data[cname] = cgen_ha
    if cgen_hv: # define general constraint (harvest volume)
        cname = 'cgen_hv'
        coeff_funcs[cname] = partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=None) 
        cgen_data[cname] = cgen_hv
    if cgen_gs: # define general constraint (growing stock)
        cname = 'cgen_gs'
        coeff_funcs[cname] = partial(cmp_c_ci, yname=tvy_name, mask=None)
        cgen_data[cname] = cgen_gs
    # if cgen_cp: # define general constraint (carbon pools)
    #     cname = 'cgen_cp'
    #     coeff_funcs[cname] = partial(cmp_c_ci, yname=cp_name, mask=None)
    #     cgen_data[cname] = cgen_cp
    # if cgen_cf: # define general constraint (carbon fluxes)
    #     cname = 'cgen_cf'
    #     coeff_funcs[cname] = partial(cmp_c_ci, yname=cf_name, mask=None)
    #     cgen_data[cname] = cgen_cf
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)

def run_scenario(fm, scenario_name='base'):
    #Import Module
    import gurobipy as grb
    
    cflw_ha = {}
    cflw_hv = {}
    cgen_ha = {}
    cgen_hv = {}
    cgen_gs = {}
    # cgen_cp = {}
    # cgen_cf = {}
    
    # define harvest area and harvest volume even-flow constraints
    cflw_ha = ({p:0.05 for p in fm.periods}, 1)
    cflw_hv = ({p:0.05 for p in fm.periods}, 1)

    if scenario_name == 'base': 
        # Base scenario
        print('running base scenario')
    else:
        assert False # bad scenario name
      
    
    p = gen_scenario(fm=fm, 
                     name=scenario_name, 
                     cflw_ha=cflw_ha, 
                     cflw_hv=cflw_hv,
                     cgen_ha=cgen_ha,
                     cgen_hv=cgen_hv,
                     cgen_gs=cgen_gs,)

    p.solver('gurobi')

    fm.reset()
    p.solve()

    if p.status() != ws3.opt.STATUS_OPTIMAL:
        print('Model not optimal.')
        df = None   
    else:
        sch = fm.compile_schedule(p)
        fm.apply_schedule(sch, 
                        force_integral_area=False, 
                        override_operability=False,
                        fuzzy_age=False,
                        recourse_enabled=False,
                        verbose=False,
                        compile_c_ycomps=True)
        
        df = compile_scenario(fm)
        fig, ax = plot_scenario(df)
    
    return df, fig, ax