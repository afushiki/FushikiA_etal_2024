"""
DLC cylinder analysis

This module contains analysis helpers for DeepLabCut (DLC) output CSV files.
"""

import csv
import os
from functools import reduce
from itertools import compress, groupby

import numpy as np
import pandas as pd


def DLC_cylinder_individual(datapath, experiment, animal, analysis_results, genotype, age, frame_start, frame_end, 
                            pcutoff, rearing_th, rearing_def_frame, spatial_resolution):
     
    """
    Create per-animal analysis outputs from a DLC export CSV for cylinder arena.
    
    Reads DLC coordinates, filters by likelihood, detects rearing events using a y-threshold, 
    and writes analysis CSV outputs.
    
    Notes:
    - 'spatial_resolution' is in cm/pixel.
    - 'rearing_def_frame' is the minimum consecutive-frame duration for a rearing event.
    """

    fps = 30                    
    DLCscorer = os.path.join(experiment + '/' + animal)
    dataname  = datapath + DLCscorer + '.csv'
    Dataframe = pd.read_csv(os.path.join(dataname), skiprows=0, header=[0,1,2], index_col=0)
    scorer    = Dataframe.columns.get_level_values(0)[0]
    bodyparts = set(Dataframe.columns.get_level_values(1))
    bodyparts = list(bodyparts) 
    Dataframe = Dataframe.iloc[frame_start:frame_end]

    # ------------------------------------------------------------------
    # Filter frames by DLC likelihood ('pcutoff')
    # ------------------------------------------------------------------
    
    # Dataframe selection beyond pcutoff (likelihood)
    Idx1 = Dataframe[scorer]['nose']['likelihood'].values > pcutoff
    Idx2 = Dataframe[scorer]['sensor']['likelihood'].values > pcutoff

    # Intersection of index of likelihood in all bodyparts
    Idxs = [Idx1, Idx2]
    Dataframe = Dataframe[np.logical_and.reduce(Idxs, axis=0)]

    if genotype == 'MitoPark':
        df0 = Dataframe[scorer]['nose'][Dataframe[scorer]['nose']['y'] < rearing_th] 
        df0['Frames'] = df0.index

        # Make bool list: rearing duration greater than or equal to 15 frames (0.5 second)
        bool_list = []
        rearing_list = []       
        for k, g in groupby(enumerate(df0.Frames), lambda x : x[0] - x[1]): #df0=nose
            if len(dict(g).values()) >= rearing_def_frame:
                bool_list.append(True)
            else:
                bool_list.append(False)
        # Make rearing list
        for k, g in groupby(enumerate(df0.Frames), lambda x : x[0] - x[1]):  #df0=nose
            rearing_list.append(list(dict(g).values()))
    else:
        df1 = Dataframe[scorer]['sensor'][Dataframe[scorer]['sensor']['y'] < rearing_th]
        df1['Frames'] = df1.index

        # Make bool list: rearing duration greater than or equal to 15 frames (0.5 second)
        bool_list = []
        rearing_list = []
        for k, g in groupby(enumerate(df1.Frames), lambda x : x[0] - x[1]): #df1=sensor
            if len(dict(g).values()) >= rearing_def_frame:
                bool_list.append(True)
            else:
                bool_list.append(False)
        # Make rearing list
        for k, g in groupby(enumerate(df1.Frames), lambda x : x[0] - x[1]):  #df1=sensor
            rearing_list.append(list(dict(g).values()))

    rearing_frames = list(compress(rearing_list, bool_list))
   
    # For csv outputs
    # ------------------------------------------------------------------
    # Write analysis metadata + results
    # ------------------------------------------------------------------
    list1 = ['DLCscorer', 'frame_start', 'frame_end', 'fps', 'likelihood', 'rearing_th(pixel)','rearing_def_frame', 'spatial_resolution(cm/pixel)', 'number_of_rearing']
    list2 = [DLCscorer, frame_start, frame_end, fps, pcutoff, rearing_th, rearing_def_frame, spatial_resolution, len(rearing_frames)]

    csv_dir = datapath + analysis_results + '_' + str(frame_start) + '_' + str(frame_end) + '_' + str(pcutoff) + '_' + str(rearing_th) + '_' + str(rearing_def_frame) + '/'
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)
    
    if genotype == 'MitoPark':
        filename = csv_dir + 'Analysis_' + experiment + '_' + age + '_' + animal +'.csv'
    else: 
        filename = csv_dir + 'Analysis_' + experiment + '_' + animal +'.csv'

    with open(filename, 'w', newline='') as csvfile:
        writer=csv.writer(csvfile, delimiter=',')
        writer.writerow(list1)
        writer.writerow(list2)        

def DLC_cylinder_summary(datapath, rig, analysis_results, DAY, genotype, manipulation, age, 
                         frame_start, frame_end, pcutoff, rearing_th, rearing_def_frame, spatial_resolution, Control_ID, Test_ID, midbrain):
    
    """
    Create a per-day summary across animals from the individual cylinder arena outputs.
    
    Aggregates per-animal analysis into Control/Test summaries with descriptive statistics,
    and writes summary CSVs.
    """

    # Rename the index name (shorter) at dataframe
    name1_af = '' # no name
    name1_bf = '_' + DAY + '_resnet50_cylinder_modelJul23shuffle1_1030000_filtered'
    name2_af = 'Cy'

    if genotype == 'MitoPark':
        experiment = rig + '_' + genotype + manipulation  # for only MitoPark
        name2_bf = rig + '_' + genotype + '/' # for only MitoPark
    else:
        experiment = rig + '_' + genotype + '_' + manipulation # for other animals
        name2_bf = rig + '_' + genotype + '_' + manipulation + '/' # for other animals
        name3_bf = rig + '_' + genotype + '_' + 'Casp3/' # it's important for control
 
    csv_dir = datapath + analysis_results + '_' + str(frame_start) + '_' + str(frame_end) + '_' + str(pcutoff) + '_' + str(rearing_th) + '_' + str(rearing_def_frame) + '/'
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)
        
    animal = [animal for animal in os.listdir(csv_dir) if animal.startswith('Analysis_Cylinder')&animal.endswith('_filtered.csv')]
    animal.sort()

    df=[]
    for i in range(len(animal)):
        df.append(pd.read_csv(os.path.join(csv_dir + animal[i]), header=[0], index_col=0)) 

    df_init = pd.concat(df, axis=0)
    df_DAY = df_init.filter(regex = DAY, axis=0)

    # rename the index name (shorter)
    df_DAY.index = df_DAY.index.str.replace(name1_bf, name1_af)
    df_DAY.index = df_DAY.index.str.replace(name2_bf, name2_af)

    Control = df_DAY.filter(regex = Control_ID, axis=0)
    Test = df_DAY.filter(regex = Test_ID, axis=0)
    
    filter_ctrl = Control.axes[0].tolist()
    filter_test = Test.axes[0].tolist()

    Control_summary = Control.copy()
    Control_summary.loc['Mean'] = Control_summary.loc[filter_ctrl].mean(numeric_only=True, axis=0)
    Control_summary.loc['Median'] = Control_summary.loc[filter_ctrl].median(numeric_only=True, axis=0)
    Control_summary.loc['Variance'] = Control_summary.loc[filter_ctrl].var(numeric_only=True, axis=0)
    Control_summary.loc['Standard_deviation'] = Control_summary.loc[filter_ctrl].std(numeric_only=True, axis=0)
    Control_summary.loc['Standard_error'] = Control_summary.loc[filter_ctrl].sem(numeric_only=True, axis=0)
    Test_summary = Test.copy()
    Test_summary.loc['Mean'] = Test_summary.loc[filter_test].mean(numeric_only=True, axis=0)
    Test_summary.loc['Median'] = Test_summary.loc[filter_test].median(numeric_only=True, axis=0)
    Test_summary.loc['Variance'] = Test_summary.loc[filter_test].var(numeric_only=True, axis=0)
    Test_summary.loc['Standard_deviation'] = Test_summary.loc[filter_test].std(numeric_only=True, axis=0)
    Test_summary.loc['Standard_error'] = Test_summary.loc[filter_test].sem(numeric_only=True, axis=0)

    # csv outputs
    if genotype == 'DAT':
        Control_summary.to_csv((csv_dir + analysis_results + '_' + midbrain + '_' + DAY + '_Control_summary.csv'))
        Test_summary.to_csv((csv_dir + analysis_results + '_' + midbrain + '_' + DAY + '_Test_summary.csv'))
    else: 
        Control_summary.to_csv((csv_dir + analysis_results + '_' + DAY + '_Control_summary.csv'))
        Test_summary.to_csv((csv_dir + analysis_results + '_' + DAY + '_Test_summary.csv'))
