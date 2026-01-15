"""
DLC square open-field analysis

This module contains analysis helpers for DeepLabCut (DLC) output CSV files.
"""

import csv
import os
from functools import reduce
from itertools import compress, groupby

import numpy as np
import pandas as pd


def DLC_square_individual(datapath, experiment, animal, analysis_results, genotype, age, frame_start, frame_end, 
                          stop_th, stop_frame, move_th, move_frame, pcutoff, speed_outlier, spatial_resolution):
    
    """
    Create per-animal analysis outputs from a DLC export CSV for square open-field.
    
    Reads DLC coordinates, filters by likelihood, computes speed, detects stop/move
    bouts based on thresholds/duration windows, and writes analysis CSV outputs.
    
    Notes:
    - 'spatial_resolution' is in cm/pixel.
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
    Idx1 = Dataframe[scorer]['nose']['likelihood'].values >= pcutoff
    Idx2 = Dataframe[scorer]['rightear']['likelihood'].values >= pcutoff
    Idx3 = Dataframe[scorer]['leftear']['likelihood'].values >= pcutoff
    Idx4 = Dataframe[scorer]['spine1']['likelihood'].values >= pcutoff
    Idx5 = Dataframe[scorer]['spine2']['likelihood'].values >= pcutoff
    Idx6 = Dataframe[scorer]['spine3']['likelihood'].values >= pcutoff
    Idx7 = Dataframe[scorer]['tail']['likelihood'].values >= pcutoff

    # Intersection of index of likelihood in all bodyparts
    Idxs = [Idx1, Idx2, Idx3, Idx4, Idx5, Idx6, Idx7]
    Dataframe = Dataframe[np.logical_and.reduce(Idxs, axis=0)]

    df_list = []
    dist_list = []
    speed_list = dict()
    for bp in bodyparts:
        df = Dataframe[scorer][bp]
        dist = df.diff().fillna(0.)
        dist['Dist_pixel'] = np.sqrt(dist.x**2 + dist.y**2) # pixel
        dist['Dist_cm'] = np.sqrt(dist.x**2 + dist.y**2)*(spatial_resolution)
        dist['Frames'] = dist.index
        dist['Time'] = dist.index/fps # second
        dist['Speed'] = dist.Dist_cm/(1/fps) # cm/s
        df_list.append(df)
        dist_list.append(dist)
        for j in  dist_list: 
            speed_list[bp]=pd.DataFrame(j)

    # Speedlist
    nose_speed     = speed_list['nose']['Speed']
    rightear_speed = speed_list['rightear']['Speed']
    leftear_speed  = speed_list['leftear']['Speed']
    spine1_speed   = speed_list['spine1']['Speed']
    spine2_speed   = speed_list['spine2']['Speed']
    spine3_speed   = speed_list['spine3']['Speed']
    tail_speed     = speed_list['tail']['Speed']

    speed_list['nose']     = speed_list['nose'].drop(speed_list['nose'][speed_list['nose']['Speed'] >= speed_outlier].index)
    speed_list['rightear'] = speed_list['rightear'].drop(speed_list['rightear'][speed_list['rightear']['Speed'] >= speed_outlier].index)
    speed_list['leftear']  = speed_list['leftear'].drop(speed_list['leftear'][speed_list['leftear']['Speed'] >= speed_outlier].index)
    speed_list['spine1']   = speed_list['spine1'].drop(speed_list['spine1'][speed_list['spine1']['Speed'] >= speed_outlier].index)
    speed_list['spine2']   = speed_list['spine2'].drop(speed_list['spine2'][speed_list['spine2']['Speed'] >= speed_outlier].index)
    speed_list['spine3']   = speed_list['spine3'].drop(speed_list['spine3'][speed_list['spine3']['Speed'] >= speed_outlier].index)
    speed_list['tail']     = speed_list['tail'].drop(speed_list['tail'][speed_list['tail']['Speed'] >= speed_outlier].index)

    #---For stop---#
    # Only picked up the frames < stop_threshold(cm/s) in all bodyparts to find stop frames
    # Individual stop dataframe of nose ~ tail
    df0_stop = speed_list['nose'].loc[nose_speed < stop_th]
    df1_stop = speed_list['rightear'].loc[rightear_speed < stop_th]
    df2_stop = speed_list['leftear'].loc[leftear_speed < stop_th]
    df3_stop = speed_list['spine1'].loc[spine1_speed < stop_th]
    df4_stop = speed_list['spine2'].loc[spine2_speed < stop_th]
    df5_stop = speed_list['spine3'].loc[spine3_speed < stop_th]
    df6_stop = speed_list['tail'].loc[tail_speed < stop_th]

    # Add unique column names except Time for merge dataframes below
    df0_stop = df0_stop.rename(columns={c: 'nose_'+c for c in df0_stop.columns if c not in ['Time']})
    df1_stop = df1_stop.rename(columns={c: 'rightear_'+c for c in df1_stop.columns if c not in ['Time']})
    df2_stop = df2_stop.rename(columns={c: 'leftear_'+c for c in df2_stop.columns if c not in ['Time']})
    df3_stop = df3_stop.rename(columns={c: 'spine1_'+c for c in df3_stop.columns if c not in ['Time']})
    df4_stop = df4_stop.rename(columns={c: 'spine2_'+c for c in df4_stop.columns if c not in ['Time']})
    df5_stop = df5_stop.rename(columns={c: 'spine3_'+c for c in df5_stop.columns if c not in ['Time']})
    df6_stop = df6_stop.rename(columns={c: 'tail_'+c for c in df6_stop.columns if c not in ['Time']})

    dfs_stop = [df0_stop, df1_stop, df2_stop, df3_stop, df4_stop, df5_stop, df6_stop]
    df_final_stop = reduce(lambda left,right: pd.merge(left,right,on='Time'), dfs_stop) 

    # Make bool list: stop duration greater than or equal to 15 frames (0.5 second) for all body parts
    bool_stop_list = []
    stop_list = []
    for k, g in groupby(enumerate(df_final_stop.nose_Frames), lambda x : x[0] - x[1]): 
        if len(dict(g).values()) >= stop_frame:
            bool_stop_list.append(True)
        else:
            bool_stop_list.append(False)     
    # Make stop list: picked up stopped frames.
    for k, g in groupby(enumerate(df_final_stop.nose_Frames), lambda x : x[0] - x[1]):
        stop_list.append(list(dict(g).values())) # pick up all consecutive frames 
    stop_frames = list(compress(stop_list, bool_stop_list))
    stop_frames_flatten = [i for j in stop_frames for i in j]

    stop_frames_len_list = []
    for x in stop_frames:
        stop_frames_len_list.append(len(x))
    stop_frames_len_list[:] = [x / fps for x in stop_frames_len_list] # frame to second

    #---For move---#
    # Only picked up the frames > move_threshold(cm/s) in all bodyparts to find move frames
    # Individual move dataframe of nose ~ tail
    df0_move = speed_list['nose'].loc[nose_speed >= move_th]
    df1_move = speed_list['rightear'].loc[rightear_speed >= move_th]
    df2_move = speed_list['leftear'].loc[leftear_speed >= move_th]
    df3_move = speed_list['spine1'].loc[spine1_speed >= move_th]
    df4_move = speed_list['spine2'].loc[spine2_speed >= move_th]
    df5_move = speed_list['spine3'].loc[spine3_speed >= move_th]
    df6_move = speed_list['tail'].loc[tail_speed >= move_th]

    # Add unique column names except Time for merge dataframes below
    df0_move = df0_move.rename(columns={c: 'nose_'+c for c in df0_move.columns if c not in ['Time']})
    df1_move = df1_move.rename(columns={c: 'rightear_'+c for c in df1_move.columns if c not in ['Time']})
    df2_move = df2_move.rename(columns={c: 'leftear_'+c for c in df2_move.columns if c not in ['Time']})
    df3_move = df3_move.rename(columns={c: 'spine1_'+c for c in df3_move.columns if c not in ['Time']})
    df4_move = df4_move.rename(columns={c: 'spine2_'+c for c in df4_move.columns if c not in ['Time']})
    df5_move = df5_move.rename(columns={c: 'spine3_'+c for c in df5_move.columns if c not in ['Time']})
    df6_move = df6_move.rename(columns={c: 'tail_'+c for c in df6_move.columns if c not in ['Time']})

    dfs_move = [df0_move, df1_move, df2_move, df3_move, df4_move, df5_move, df6_move]
    df_final_move = reduce(lambda left,right: pd.merge(left,right,on='Time'), dfs_move) 

    # Make bool list: move duration greater than or equal to 15 frames (0.5 second) for all body parts
    bool_move_list = []
    move_list = []
    for k, g in groupby(enumerate(df_final_move.nose_Frames), lambda x : x[0] - x[1]): 
        if len(dict(g).values()) >= move_frame:
            bool_move_list.append(True)
        else:
            bool_move_list.append(False)     
    # Make move list: picked up moveped frames.
    for k, g in groupby(enumerate(df_final_move.nose_Frames), lambda x : x[0] - x[1]):
        move_list.append(list(dict(g).values())) # pick up all consecutive frames 
    move_frames = list(compress(move_list, bool_move_list))
    move_frames_flatten = [i for j in move_frames for i in j]

    move_frames_len_list = []
    for x in move_frames:
        move_frames_len_list.append(len(x))
    move_frames_len_list[:] = [x / fps for x in move_frames_len_list] # frame to second

    total_immobile_frames = sum([len(i) for i in stop_frames])
    total_mobile_frames   = sum([len(i) for i in move_frames])
    immobile_times = total_immobile_frames/fps
    mobile_times   = total_mobile_frames/fps
    stop_bout_length = 0 if len(stop_frames_len_list) == 0 else (sum(stop_frames_len_list)/len(stop_frames_len_list)) 
    move_bout_length = 0 if len(move_frames_len_list) == 0 else (sum(move_frames_len_list)/len(move_frames_len_list)) 
    immobile_fraction = immobile_times/(immobile_times + mobile_times)*100 
    mobile_fraction = mobile_times/(immobile_times + mobile_times)*100 
    
    # spine2 is a centroid of the animal
    immobile_df_spine2 = speed_list['spine2'].loc[stop_frames_flatten] 
    mobile_df_spine2   = speed_list['spine2'].loc[move_frames_flatten]
    distance_traveled_spine2 = mobile_df_spine2['Dist_cm'].sum()

    if mobile_times !=0:
        average_speed_spine2 = distance_traveled_spine2 / mobile_times
    else: 
        average_speed_spine2 = 0

    no_of_stops = len(stop_frames)
    stop_events = [len(i) for i in stop_frames]
    no_of_stops_1s = len([i for i in stop_events if i >= 30])
    no_of_stops_2s = len([i for i in stop_events if i >= 60]) 
    no_of_stops_3s = len([i for i in stop_events if i >= 90]) 
    
    immobile_speed_list = immobile_df_spine2['Speed'].tolist()
    mobile_speed_list = mobile_df_spine2['Speed'].tolist()
    
    # For csv outputs
    # ------------------------------------------------------------------
    # Write analysis metadata + results
    # ------------------------------------------------------------------
    list1 = ['DLCscorer', 'frame_start', 'frame_end', 'fps','stop_th(cm/s)','stop_frame', 'move_th(cm/s)','move_frame',
             'likelihood(pcutoff)', 'speed_outlier(cm/s)',
             'immobile_times(s)','mobile_times(s)', 'immobile_fraction', 'mobile_fraction',
             'distance_traveled_spine2(cm)', 'average_speed_spine2(cm/s)', 
             'no_of_stops', 'no_of_stops_1s','no_of_stops_2s', 'no_of_stops_3s', 'stop_bout_length(s)','move_bout_length(s)',
             'stop_frames_len_list', 'move_frames_len_list', 
             'mobile_speed_list', 'immobile_speed_list']
    list2 = [DLCscorer, frame_start, frame_end, fps, stop_th, stop_frame, move_th, move_frame, 
             pcutoff, speed_outlier,
             immobile_times, mobile_times, immobile_fraction, mobile_fraction,  
             distance_traveled_spine2, average_speed_spine2, 
             no_of_stops, no_of_stops_1s, no_of_stops_2s, no_of_stops_3s, stop_bout_length, move_bout_length,
             stop_frames_len_list, move_frames_len_list, 
             mobile_speed_list, immobile_speed_list]

    csv_dir = datapath + analysis_results + '_' + str(frame_start) + '_' + str(frame_end) + '_' + str(stop_th) + '_' + str(stop_frame) + '_' + str(move_th) + '_' + str(move_frame) + '_' + str(pcutoff) + '/'
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)
    
    if genotype == 'MitoPark':
        filename = csv_dir + 'Analysis_' + experiment + '_' + age + '_' + animal + '.csv'
    else: 
        filename = csv_dir + 'Analysis_' + experiment + '_' + animal + '.csv'
    
    import sys
    csv.field_size_limit(sys.maxsize)
    with open(filename, 'w', newline='') as csvfile:
        writer=csv.writer(csvfile, delimiter=',')
        writer.writerow(list1)
        writer.writerow(list2)

def DLC_square_summary(datapath, rig, analysis_results, DAY, genotype, manipulation, age, frame_start, frame_end, 
                       stop_th, stop_frame, move_th, move_frame, pcutoff, speed_outlier, spatial_resolution, Control_ID, Test_ID, midbrain):
    
    """
    Create a per-day summary across animals from the individual square open-field outputs.
    
    Aggregates per-animal analysis into Control/Test summaries with descriptive statistics,
    and writes summary CSVs.
    """

    # Rename the index name (shorter) at dataframe
    name1_af = '' # no name
    name1_bf = '_' + DAY + '_resnet50_squareOF_modelApr11shuffle1_1030000_filtered'
    name2_af = 'Sq'

    if genotype == 'MitoPark':
        experiment = rig + '_' + genotype + manipulation  # for only MitoPark
        name2_bf = rig + '_' + genotype + '/' 
    else:
        experiment = rig + '_' + genotype + '_' + manipulation # for other animals
        name2_bf = rig + '_' + genotype + '_' + manipulation + '/' 
        name3_bf = rig + '_' + genotype + '_' + 'Casp3/' # it's important for control
    
    csv_dir = datapath + analysis_results + '_' + str(frame_start) + '_' + str(frame_end) + '_' + str(stop_th) + '_' + str(stop_frame) + '_' + str(move_th) + '_' + str(move_frame) + '_' + str(pcutoff) + '/' 
    
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)
    
    animal = [animal for animal in os.listdir(csv_dir) if animal.startswith('Analysis_SquareOF')&animal.endswith('_filtered.csv')]
    animal.sort()

    df=[]
    for i in range(len(animal)):
        df.append(pd.read_csv(os.path.join(csv_dir + animal[i]), header=[0], index_col=0)) 

    df_init = pd.concat(df, axis=0)
    df_DAY = df_init.filter(regex = DAY, axis=0)
    
    if genotype == 'MitoPark':
    # rename the index name (shorter)
        df_DAY.index = df_DAY.index.str.replace(name1_bf, name1_af)
        df_DAY.index = df_DAY.index.str.replace(name2_bf, name2_af)
    else: 
        # rename the index name (shorter)
        df_DAY.index = df_DAY.index.str.replace(name1_bf, name1_af)
        df_DAY.index = df_DAY.index.str.replace(name2_bf, name2_af)
        df_DAY.index = df_DAY.index.str.replace(name3_bf, name2_af)

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
