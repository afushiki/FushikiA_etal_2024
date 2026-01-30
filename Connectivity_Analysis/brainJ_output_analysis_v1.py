import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.ticker as mtick
import seaborn as sns
import re

def rgb_to_hex(red, green, blue):
    return '#%02x%02x%02x' % (red, green, blue)

def create_rabies_df_for_figure_reorder_cells(data_dir, dir_type, atlas_type, filename, cell_counts_RL, ABA_fiber_tract, inj_site, col_df, cell_threshold):
    df1 = []; df2 = []; df2_temp = []; df3 = []
    for i in range(len(filename)):
        df1.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/data' + '/' + filename[i])))
        df1[i] = df1[i].drop(df1[i][df1[i]['Acronym'] == inj_site].index)
    for i in range(len(filename)):
        df1[i] = df1[i].sort_values(by=['Graph_Order'])
        df2_temp.append(df1[i][[cell_counts_RL, 'ID']])
        df2.append(df2_temp[i].rename(columns={cell_counts_RL: 'Cells_' + filename[i].split('_')[0:3][2]}))
    df3 = pd.concat(df2, axis=1)
    filter_col = [col for col in df3 if col.startswith('Cells_')]
    df3 = pd.concat([df3[filter_col], df3.iloc[:, 1]], axis=1)
    df4 = pd.merge(df3, ABA_fiber_tract, left_on='ID', right_on='output_id', how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
    df4 = df4.dropna(axis=1, how='all')
    df4['mean'] = df4[filter_col].mean(numeric_only=True, axis=1)
    df4['median'] = df4[filter_col].median(numeric_only=True, axis=1)
    df4 = df4.merge(col_df, left_on='ID', right_on='output_id')
    df4['Pass_Status'] = np.logical_and.reduce(df4[filter_col] >= cell_threshold, axis=1)
    df5 = df4[df4['Pass_Status'] == True]
    df5 = df5.sort_values(by=['graph_order_changed'])
    df5['acronym+ID'] = df5['acronym'] + ' (ID:' + df5['ID'].astype(str) + ')'
    return df5

def align_brain_region_cells(df_ipsi, df_contra, col_df, atlas_type):
    df_align = pd.merge(df_ipsi, df_contra, how='outer', on=['acronym', 'acronym'])
    filter_col = [col for col in df_ipsi if col.startswith('Cells_')]
    filter_col_x = [sub + '_x' for sub in filter_col]
    filter_col_y = [sub + '_y' for sub in filter_col]
    df_align_x = df_align[filter_col_x + ['ID_x', 'mean_x', 'median_x', 'hex_color_x', 'name_x', 'acronym']]
    df_align_y = df_align[filter_col_y + ['ID_y', 'mean_y', 'median_y', 'hex_color_y', 'name_y', 'acronym']]
    df_align_x.columns = df_align_x.columns.str.rstrip('_x')
    df_align_y.columns = df_align_y.columns.str.rstrip('_y')
    df_both = [df_align_x, df_align_y]
    df_both_lf = []
    for i in range(len(df_both)):
        df_both[i].loc[:, ['hex_color']] = df_both[i].loc[:, ['hex_color']].fillna('#' + '000000')
        df_both_lf.append(pd.melt(df_both[i][filter_col + ['ID', 'acronym']], id_vars='acronym', value_vars=filter_col))
        col_df_2 = col_df.drop(['output_id', 'hex_color', 'name'], axis=1)
        df_both[i] = df_both[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both[i] = df_both[i].sort_values(by=['graph_order_changed'])
        df_both_lf[i] = df_both_lf[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both_lf[i] = df_both_lf[i].sort_values(by=['graph_order_changed'])
    df_both[1]['ID'] = df_both[0]['ID']
    return (df_both, df_both_lf)

def create_proj_df_for_figure_reorder_projections(data_dir, dir_type, atlas_type, filename, density, ABA_fiber_tract, inj_site, col_df, proj_threshold):
    df1 = []; df2 = []; df2_temp = []; df3 = []
    for i in range(len(filename)):
        df1.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/data' + '/' + filename[i])))
        df1[i] = df1[i].drop(df1[i][df1[i]['Acronym'] == inj_site].index)
    for i in range(len(filename)):
        df1[i] = df1[i].sort_values(by=['Graph_Order'])
        df2_temp.append(df1[i][[density, 'ID']])
        df2.append(df2_temp[i].rename(columns={density: 'Projections_' + filename[i].split('_')[0:3][2]}))
    df3 = pd.concat(df2, axis=1)
    filter_col = [col for col in df3 if col.startswith('Projections_')]
    df3 = pd.concat([df3[filter_col], df3.iloc[:, 1]], axis=1)
    df4 = pd.merge(df3, ABA_fiber_tract, left_on='ID', right_on='output_id', how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
    df4 = df4.dropna(axis=1, how='all')
    df4['mean'] = df4[filter_col].mean(numeric_only=True, axis=1)
    df4['median'] = df4[filter_col].median(numeric_only=True, axis=1)
    df4 = df4.merge(col_df, left_on='ID', right_on='output_id')
    df4['Pass_Status'] = np.logical_and.reduce(df4[filter_col] >= proj_threshold, axis=1)
    df5 = df4[df4['Pass_Status'] == True]
    df5 = df5.sort_values(by=['graph_order_changed'])
    return df5

def align_brain_region_projections(df_ipsi, df_contra, col_df, atlas_type):
    df_align = pd.merge(df_ipsi, df_contra, how='outer', on=['acronym', 'acronym'])
    filter_col = [col for col in df_ipsi if col.startswith('Projections_')]
    filter_col_x = [sub + '_x' for sub in filter_col]
    filter_col_y = [sub + '_y' for sub in filter_col]
    df_align_x = df_align[filter_col_x + ['ID_x', 'mean_x', 'median_x', 'hex_color_x', 'name_x', 'acronym']]
    df_align_y = df_align[filter_col_y + ['ID_y', 'mean_y', 'median_y', 'hex_color_y', 'name_y', 'acronym']]
    df_align_x.columns = df_align_x.columns.str.rstrip('_x')
    df_align_y.columns = df_align_y.columns.str.rstrip('_y')
    df_both = [df_align_x, df_align_y]
    df_both_lf = []
    for i in range(len(df_both)):
        df_both[i].loc[:, ['hex_color']] = df_both[i].loc[:, ['hex_color']].fillna('#' + '000000')
        df_both_lf.append(pd.melt(df_both[i][filter_col + ['ID', 'acronym']], id_vars='acronym', value_vars=filter_col))
        col_df_2 = col_df.drop(['output_id', 'hex_color', 'name'], axis=1)
        df_both[i] = df_both[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both[i] = df_both[i].sort_values(by=['graph_order_changed'])
        df_both_lf[i] = df_both_lf[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both_lf[i] = df_both_lf[i].sort_values(by=['graph_order_changed'])
    df_both[1]['ID'] = df_both[0]['ID']
    return (df_both, df_both_lf)

def create_final_bar_graph(graphtype, save_dir, df_both, animal, projection, atlas_type, channel, inj_site, average, cell_or_density, threshold, ylim, figsize):
    if graphtype == 'Cells':
        col_name = 'Cells_'
        var_name = 'cells'
    if graphtype == 'Cells_prop':
        col_name = 'Cells_'
        var_name = 'cells'
    elif graphtype == 'Projections':
        col_name = 'Projections_'
        var_name = 'projections'
    x_axis_names = 'acronym'
    
    title = []; df_long = []
    for i in range(len(df_both)):
        df_long.append(pd.melt(df_both[i], id_vars=['ID', 'hex_color', 'name', x_axis_names, 'graph_order', 'graph_order_changed'], value_vars=[col for col in df_both[i].columns if col_name in col], var_name=var_name, value_name='value'))
        fig = plt.figure(figsize=figsize, dpi=300)
        sns.barplot(data=df_long[i], x=df_long[i][x_axis_names], y=df_long[i]['value'], hue=x_axis_names, palette=list(df_long[i]['hex_color'])[:len(df_both[0].index)], capsize=0.1, err_kws={'color': 'grey', 'linewidth': 0.5})
        sns.stripplot(data=df_long[i], x=df_long[i][x_axis_names], y=df_long[i]['value'], size=8, alpha=0.5, dodge=True, jitter=True, color='black', legend=False)
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_ylim(0, ylim)
        ax.set_xticks(np.arange(0, len(df_both[0].index), 1))
        ax.set_xticklabels(df_both[i][x_axis_names], size=20, rotation=90, fontstyle='normal')
        ax.yaxis.set_major_formatter(ticker.EngFormatter())
        sns.despine(left=True)
        colors = df_both[i]['hex_color']
        for xtick, color in zip(ax.get_xticklabels(), colors):
            xtick.set_color(color)
        title.append(animal + 'dopamine neuron subtype:' + '_' + projection[i] + '_' + col_name)
        plt.title(title[i].replace('_', ' '), fontsize=25)
        plt.xlabel('Brain region', fontsize=25)
        if graphtype == 'Cells':
            plt.ylabel(cell_or_density.replace('_', ' '), fontsize=25)
            plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_reorder.pdf'), bbox_inches='tight', dpi=300)
            df_both[i].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_reorder.csv'), index=False)
            plt.show()
        elif graphtype == 'Cells_prop':
            plt.ylabel('Proportion of total' + '(%)', fontsize=25)
            plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_Proportion' + '_reorder.pdf'), bbox_inches='tight', dpi=300)
            df_both[i].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_Proportion' + '_reorder.csv'), index=False)
            plt.show()
        elif graphtype == 'Projections':
            plt.ylabel(cell_or_density.replace('_', ' ') + '(%)', fontsize=25)
            plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_reorder.pdf'), bbox_inches='tight', dpi=300)
            df_both[i].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animal + channel + average + '_' + str(threshold) + '_' + projection[i] + '_' + cell_or_density + '_reorder.csv'), index=False)
            plt.show()

def create_rabies_df_sum_reorder_cells(data_dir, dir_type, atlas_type, filename, cell_counts_RL, ABA_fiber_tract, inj_site, col_df):
    df1 = []; df2 = []; df2_temp = []; df3 = []
    for i in range(len(filename)):
        df1.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/data' + '/' + filename[i])))
        df1[i] = df1[i].drop(df1[i][df1[i]['Acronym'] == inj_site].index)
    for i in range(len(filename)):
        df1[i] = df1[i].sort_values(by=['Graph_Order'])
        df2_temp.append(df1[i][[cell_counts_RL, 'ID']])
        df2.append(df2_temp[i].rename(columns={cell_counts_RL: 'Cells_' + filename[i].split('_')[0:3][2]}))
    df3 = pd.concat(df2, axis=1)
    filter_col = [col for col in df3 if col.startswith('Cells_')]
    df3 = pd.concat([df3[filter_col], df3.iloc[:, 1]], axis=1)
    df4 = pd.merge(df3, ABA_fiber_tract, left_on='ID', right_on='output_id', how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
    df4 = df4.dropna(axis=1, how='all')
    df4.loc['Sum_one_side'] = df4[filter_col].sum()
    return df4

def create_rabies_df_for_figure_prop_cells(df, cell_counts_RL, col_df, edge_threshold):
    filter_col = [col for col in df if col.startswith('Cells_')]
    df_prop = df.loc[:, filter_col].div(df.loc['Total'], axis=1) * 100
    df_prop['ID'] = df['ID']
    df_prop = df_prop.drop(['Sum_one_side', 'Total'], axis=0)
    df_prop['mean'] = df_prop[filter_col].mean(numeric_only=True, axis=1)
    df_prop['median'] = df_prop[filter_col].median(numeric_only=True, axis=1)
    df_prop = df_prop.merge(col_df, left_on='ID', right_on='output_id')
    df_prop['Side'] = cell_counts_RL.split('_')[2]
    df_prop['Pass_Status'] = np.logical_and.reduce(df_prop[filter_col] >= edge_threshold, axis=1)
    df_prop2 = df_prop[df_prop['Pass_Status'] == True]
    df_prop2 = df_prop2.sort_values(by=['graph_order_changed'])
    return (filter_col, df_prop2)

def create_figures_MP(df_both, df_both_lf, average, ylim, save_dir, animal, projection, atlas_type, density):
    title = []
    fig = plt.figure(figsize=(20, 5), dpi=300)
    sns.barplot(data=df_both[0], x=df_both[0]['acronym'], y=df_both[0][average], palette=df_both[0]['hex_color'], capsize=0.1)
    sns.swarmplot(data=df_both_lf[0], x=df_both_lf[0]['acronym'], y=df_both_lf[0]['value'], color='0', alpha=0.7)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_ylim(0, ylim)
    ax.set_xticklabels(df_both[0]['acronym'], size=25, rotation=90, fontstyle='normal')
    ax.yaxis.set_major_formatter(ticker.EngFormatter())
    sns.despine(left=True)
    colors = df_both[0]['hex_color']
    for xtick, color in zip(ax.get_xticklabels(), colors):
        xtick.set_color(color)
    title.append(animal + '_dopamine neuron:' + '_' + projection[0] + '_' + 'Projections')
    plt.title(title[0].replace('_', ' '), fontsize=25)
    plt.xlabel('Brain region', fontsize=25)
    plt.ylabel(density.replace('_', ' ') + '(%)', fontsize=25)
    plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + projection[0] + '_' + density + '_MP_reorder.pdf'), bbox_inches='tight', dpi=300)
    df_both[0].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + projection[0] + '_' + density + '_MP_reorder.csv'), index=False)
    plt.show()

def create_figures_MP_cells_both(df_both, df_both_lf, average, ylim, save_dir, animal, cells, atlas_type):
    title = []
    fig = plt.figure(figsize=(20, 5), dpi=300)
    sns.barplot(data=df_both[0], x=df_both[0]['acronym'], y=df_both[0][average], palette=df_both[0]['hex_color'], capsize=0.1)
    sns.swarmplot(data=df_both_lf[0], x=df_both_lf[0]['acronym'], y=df_both_lf[0]['value'], color='0', alpha=0.7)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_ylim(0, ylim)
    ax.set_xticklabels(df_both[0]['acronym'], size=25, rotation=90, fontstyle='normal')
    ax.yaxis.set_major_formatter(ticker.EngFormatter())
    sns.despine(left=True)
    colors = df_both[0]['hex_color']
    for xtick, color in zip(ax.get_xticklabels(), colors):
        xtick.set_color(color)
    title.append(animal + '_dopamine neuron:' + '_' + 'Combined' + '_' + 'Cells')
    plt.title(title[0].replace('_', ' '), fontsize=25)
    plt.xlabel('Brain region', fontsize=25)
    plt.ylabel('Total Cells', fontsize=25)
    plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + 'Combined_Total_Cells' + '_MP_reorder_both_ave.pdf'), bbox_inches='tight', dpi=300)
    df_both[0].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + 'Combined_Total_Cells' + '_MP_reorder_both_ave.csv'), index=False)
    plt.show()

def create_figures_MP_projections_both(df_both, df_both_lf, average, ylim, save_dir, animal, projection, atlas_type):
    title = []
    fig = plt.figure(figsize=(20, 5), dpi=300)
    sns.barplot(data=df_both[0], x=df_both[0]['acronym'], y=df_both[0][average], palette=list(df_both[0]['hex_color']), hue=df_both[0]['acronym'], capsize=0.1)
    sns.swarmplot(data=df_both_lf[0], x=df_both_lf[0]['acronym'], y=df_both_lf[0]['value'], color='0', alpha=0.7, size=10)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_ylim(0, ylim)
    ax.set_xticklabels(df_both[0]['acronym'], size=35, rotation=90, fontstyle='normal')
    ax.yaxis.set_major_formatter(ticker.EngFormatter())
    sns.despine(left=True)
    colors = df_both[0]['hex_color']
    for xtick, color in zip(ax.get_xticklabels(), colors):
        xtick.set_color(color)
    title.append(animal + '_dopamine neuron:' + '_' + 'Combined' + '_' + 'Projections')
    plt.title(title[0].replace('_', ' '), fontsize=25)
    plt.xlabel('Striatal region', fontsize=30)
    plt.ylabel('Projection Density(%)', fontsize=30)
    plt.savefig(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + 'Combined_Projection_Density' + '_MP_reorder_both_ave.pdf'), bbox_inches='tight', dpi=300)
    df_both[0].to_csv(os.path.join(save_dir + '/BrainJ_' + atlas_type + '_' + animal + '_' + 'Combined_Projection_Density' + '_MP_reorder_both_ave.csv'), index=False)
    plt.show()

def align_brain_region_projections_both(df_both):
    filter_col = [col for col in df_both if col.startswith('Projection_')]
    df_both = [df_both]
    df_both_lf = []
    for i in range(len(df_both)):
        df_both[i].loc[:, ['hex_color']] = df_both[i].loc[:, ['hex_color']].fillna('#' + '000000')
        df_both_lf.append(pd.melt(df_both[i][filter_col + ['ID', 'acronym']], id_vars='acronym', value_vars=filter_col))
    return (df_both, df_both_lf)

def align_brain_region_cells_both(df_both, col_df):
    filter_col = [col for col in df_both if col.startswith('Total_Cells_')]
    df_both = [df_both]
    df_both_lf = []
    for i in range(len(df_both)):
        df_both[i].loc[:, ['hex_color']] = df_both[i].loc[:, ['hex_color']].fillna('#' + '000000')
        df_both_lf.append(pd.melt(df_both[i][filter_col + ['ID', 'acronym']], id_vars='acronym', value_vars=filter_col))
        col_df_2 = col_df.drop(['output_id', 'hex_color', 'name'], axis=1)
        df_both[i] = df_both[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both[i] = df_both[i].sort_values(by=['graph_order_changed'])
        df_both_lf[i] = df_both_lf[i].merge(col_df_2, left_on='acronym', right_on='acronym')
        df_both_lf[i] = df_both_lf[i].sort_values(by=['graph_order_changed'])
    return (df_both, df_both_lf)

def create_proj_df_for_figure_reorder_projections_both(data_dir, dir_type, atlas_type, filename, ABA_fiber_tract, inj_site, col_df, proj_threshold):
    df1 = []; df2 = []; df2_temp = []; df3 = []; df2_right = []; df2_left = []; df2_right_temp = []; df2_left_temp = []; df2_both = []
    for i in range(len(filename)):
        df1.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/data' + '/' + filename[i])))
        df1[i] = df1[i].drop(df1[i][df1[i]['Acronym'] == inj_site].index)
    for i in range(len(filename)):
        df1[i] = df1[i].sort_values(by=['Graph_Order'])
        df2_right.append(df1[i][['Projection_Density_Right', 'ID']])
        df2_left.append(df1[i][['Projection_Density_Left', 'ID']])
        df2_both.append(pd.merge(df2_right[i], df2_left[i], on='ID'))
        df2_both[i]['Projection_Both_ave'] = (df2_both[i]['Projection_Density_Right'] + df2_both[i]['Projection_Density_Left']) / 2
        df2_both[i] = df2_both[i].drop(['Projection_Density_Right', 'Projection_Density_Left'], axis=1)
    df3 = pd.concat(df2_both, axis=1)
    filter_col = [col for col in df3 if col.startswith('Projection_')]
    df3 = pd.concat([df3.iloc[:, 0], df3[filter_col]], axis=1)
    df4 = pd.merge(df3, ABA_fiber_tract, left_on='ID', right_on='output_id', how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
    df4 = df4.dropna(axis=1, how='all')
    df4['mean'] = df4[filter_col].mean(numeric_only=True, axis=1)
    df4['median'] = df4[filter_col].median(numeric_only=True, axis=1)
    df4 = df4.merge(col_df, left_on='ID', right_on='output_id')
    df4['Pass_Status'] = np.logical_and.reduce(df4[filter_col] >= proj_threshold, axis=1)
    df5 = df4[df4['Pass_Status'] == True]
    df5 = df5.sort_values(by=['graph_order_changed'])
    return df5

def create_proj_df_for_figure_cells_both(data_dir, dir_type, atlas_type, filename, ABA_fiber_tract, inj_site, col_df, cell_threshold):
    df1 = []; df2 = []; df2_temp = []; df3 = []; df2_right = []; df2_left = []; df2_right_temp = []; df2_left_temp = []; df2_both = []
    for i in range(len(filename)):
        df1.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/data' + '/' + filename[i])))
    for i in range(len(filename)):
        df1[i] = df1[i].sort_values(by=['Graph_Order'])
        df2_right.append(df1[i][['Total_Cells_Right', 'ID']])
        df2_left.append(df1[i][['Total_Cells_Left', 'ID']])
        df2_both.append(pd.merge(df2_right[i], df2_left[i], on='ID'))
        df2_both[i]['Total_Cells_Both_ave'] = (df2_both[i]['Total_Cells_Right'] + df2_both[i]['Total_Cells_Left']) / 2
        df2_both[i] = df2_both[i].drop(['Total_Cells_Right', 'Total_Cells_Left'], axis=1)
    df3 = pd.concat(df2_both, axis=1)
    filter_col = [col for col in df3 if col.startswith('Total_Cells_')]
    df3 = pd.concat([df3.iloc[:, 0], df3[filter_col]], axis=1)
    df4 = pd.merge(df3, ABA_fiber_tract, left_on='ID', right_on='output_id', how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
    df4 = df4.dropna(axis=1, how='all')
    df4['mean'] = df4[filter_col].mean(numeric_only=True, axis=1)
    df4['median'] = df4[filter_col].median(numeric_only=True, axis=1)
    df4 = df4.merge(col_df, left_on='ID', right_on='output_id')
    df4['Pass_Status'] = np.logical_and.reduce(df4[filter_col] >= cell_threshold, axis=1)
    df5 = df4[df4['Pass_Status'] == True]
    df5 = df5.sort_values(by=['graph_order'])
    return df5

def create_data_fig_ratio(data_dir, save_dir, dir_type, atlas_type, average, animals, inj_site, cell_threshold, projection):
    filename = [fname for fname in os.listdir(os.path.join(data_dir + dir_type + atlas_type + '/output')) if fname.endswith(os.path.join('_' + str(cell_threshold) + '_' + projection + '_Total_Cells_reorder' + '.csv'))]
    filename = [fname for fname in filename if inj_site in fname]
    filename.sort()
    df = []
    for i in range(len(animals)):
        df.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/output' + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animals[i] + '_C2_' + average + '_' + str(cell_threshold) + '_' + projection + '_Total_Cells_reorder' + '.csv')))
    data = pd.DataFrame({animal: df[i]['median'].to_list() for i, animal in enumerate(animals)} | {'color': df[0]['hex_color'].to_list()} | {'graph_order': df[0]['graph_order'].to_list()} | {'graph_order_changed': df[0]['graph_order_changed'].to_list()}, index=df[0]['acronym'].to_list()).T
    data_filtered = data.loc[:, ~data.loc['DAT_Cre'].isna()]
    numbers_to_exclude = ['23', '4']
    pattern = '|'.join(numbers_to_exclude)
    data_filtered = data_filtered.loc[:, ~data_filtered.columns.str.contains(pattern)]
    data_filtered = data_filtered.loc[:, ~data_filtered.columns.str.contains('SNc|VTA')]
    for i, animal in enumerate(animals):
        n_cells = df[i].columns.str.contains('Cells_').sum()
        new_name = f'{animal} (n={n_cells})' if i != 0 else f'DAT (n={n_cells})'
        data_filtered.rename(index={animal: new_name}, inplace=True)
    data_fig = data_filtered[0:len(animals)]
    data_fig = data_fig.apply(pd.to_numeric, errors='coerce')
    ratio = data_fig.div(data_fig.iloc[0])
    ratio = ratio.iloc[1:]
    color_df = data_filtered.loc['color']
    graph_order = data_filtered.loc['graph_order']
    graph_order_changed = data_filtered.loc['graph_order_changed']
    return (data_fig, ratio, color_df, data_filtered)

def create_data_fig_ratio_proj(data_dir, save_dir, dir_type, atlas_type, average, animals, inj_site, proj_threshold, projection):
    filename = [fname for fname in os.listdir(os.path.join(data_dir + dir_type + atlas_type + '/output')) if fname.endswith(os.path.join('_' + str(proj_threshold) + '_' + projection + '_Projection_Density_reorder' + '.csv'))]
    filename = [fname for fname in filename if inj_site in fname]
    filename.sort()
    df = []
    for i in range(len(animals)):
        df.append(pd.read_csv(os.path.join(data_dir + dir_type + atlas_type + '/output' + '/BrainJ_' + atlas_type + '_' + inj_site + '_' + animals[i] + '_C3_' + average + '_' + str(proj_threshold) + '_' + projection + '_Projection_Density_reorder' + '.csv')))
    data = pd.DataFrame({animal: df[i]['median'].to_list() for i, animal in enumerate(animals)} | {'color': df[0]['hex_color'].to_list()} | {'graph_order': df[0]['graph_order'].to_list()} | {'graph_order_changed': df[0]['graph_order_changed'].to_list()}, index=df[0]['acronym'].to_list()).T
    data_filtered = data.loc[:, ~data.loc['DAT_Cre'].isna()]
    data_filtered = data_filtered.loc[:, ~data_filtered.columns.str.contains('SNc|VTA')]
    for i, animal in enumerate(animals):
        n_projections = df[i].columns.str.contains('Projections_').sum()
        new_name = f'{animal} (n={n_projections})' if i != 0 else f'DAT (n={n_projections})'
        data_filtered.rename(index={animal: new_name}, inplace=True)
    data_fig = data_filtered[0:len(animals)]
    data_fig = data_fig.apply(pd.to_numeric, errors='coerce')
    ratio = data_fig.div(data_fig.iloc[0])
    ratio = ratio.iloc[1:]
    color_df = data_filtered.loc['color']
    graph_order = data_filtered.loc['graph_order']
    graph_order_changed = data_filtered.loc['graph_order_changed']
    return (data_fig, ratio, color_df, data_filtered)