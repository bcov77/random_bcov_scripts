#!/usr/bin/env python



import os
import sys
import json
import itertools
import re
import random
import numpy as np
import pandas as pd
import math
import glob
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("ssm_predictor_file", type=str, help="Score file from rosetta. Ends with _000X")
parser.add_argument("af2_pae_file", type=str, help="Score file from af2. Ends with _af2pred")
parser.add_argument("seq_file", type=str, help="file with 2 columns. Protein seq, protein name.")
parser.add_argument("--mpnn_npz_files", type=str, nargs="*", help="npz files from Aditya's thing")
parser.add_argument("--output_file", type=str, default="vssm_scores.sc", help="output name")
parser.add_argument("--partial_ssm", action="store_true")

parser.add_argument("--combined_pae_matrix_npz", type=str, default="", help="Tag -> pae_matrix")
parser.add_argument("--cb_dist_npzs", type=str, nargs="*", help="Only needed if using pae matrix")
parser.add_argument("--position_class_files", type=str, nargs="*", help="Only needed if using pae matrix")

args = parser.parse_args(sys.argv[1:])


# read in scores
df1 = pd.read_csv(args.ssm_predictor_file, sep='\s+')
# sort such that best designs come first
df1 = df1.sort_values('total_score')
# remove nstruct _000X. Replicates will have duplicate descriptions now
df1['description'] = df1['description'].str.slice(None, -5)
# drop duplicates but keep the first one. We're taking the best of each nstruct here
df1 = df1.drop_duplicates('description', keep="first")

# Xinru, I mean, just comment out the next two lines right?

if ( args.af2_pae_file != "None" ):
    df2 = pd.read_csv(args.af2_pae_file, sep="\s+")
    df2['description'] = df2['description'].str.slice(None, -len("_af2pred"))

    full_df = df1.merge(df2, 'inner', 'description')
else:
    full_df = df1.copy()
    full_df['pae_interaction'] = 0
    df2 = []

# At this point, you descriptions need to end with __X__X. If you have __X__X_0001, remove the _0001
full_df['ssm_parent'] = full_df['description'].str.replace("_+([0-9]+)_+[A-Z]$", "")
full_df['ssm_letter'] = full_df['description'].str.extract("_+([A-Z])$")[0]
full_df['ssm_seqpos'] = full_df['description'].str.extract("_+([0-9]+)_+[A-Z]$")[0].astype(float)



print(len(df1), len(df2), len(full_df))


all_seqs = {}
with open(args.seq_file) as f:
    for line in f:
        line = line.strip()
        if ( len(line) == 0 ):
            continue
        sp = line.split()
        assert(len(sp) == 2)
        all_seqs[sp[1]] = sp[0]

if ( not args.mpnn_npz_files is None ):
    npz_records = []
    for npz_file in args.mpnn_npz_files:
        asdf = np.load(npz_file)
        data = asdf['odds'][0]
        ssm_parent = os.path.basename(npz_file).replace(".npz", "")
        alpha = str(asdf['alphabet'])

        for seqpos0 in range(len(data)):
            for il, letter in enumerate(list(alpha)):
                if ( letter == "X" ):
                    continue
                d = {}
                d['aditya_score'] = data[seqpos0,il]
                d['description'] = ssm_parent + "__" + str(seqpos0+1) + "__" + letter

                npz_records.append(d)
    npz_df = pd.DataFrame(npz_records)

    full_df = full_df.merge(npz_df, 'inner', 'description')


if ( args.combined_pae_matrix_npz != "" ):

    pae_matrix_dict = np.load(args.combined_pae_matrix_npz)

    cb_dist_dict = {}
    for file in args.cb_dist_npzs:
        base = os.path.basename(file).replace(".npz", "")
        dat = np.load(file)['cb_dist']
        cb_dist_dict[base] = dat

    class_dict = {}
    for file in args.position_class_files:
        base = os.path.basename(file).replace(".dat", "")
        dat = pd.read_csv(file, sep="\s+")
        class_dict[base] = dat


final_dfs = []


for ssm_parent in all_seqs:
    native_seq = all_seqs[ssm_parent]

    df = full_df[full_df['ssm_parent'] == ssm_parent].copy()

    full_size = len(native_seq)*20
    print("%i / %i entries (%i%% coverage)-- %s"%(len(df), full_size, len(df)/full_size*100, ssm_parent))
    if ( len(df) == 0 ):
        continue


    if ( args.combined_pae_matrix_npz != "" ):
        all_tags = list(pae_matrix_dict)
        these_records = []
        my_re = re.compile(ssm_parent + "__([0-9]+)__([A-Z])")

        class_df = class_dict[ssm_parent]

        for tag in all_tags:
            if ( my_re.match(tag) ):
                groups = my_re.match(tag).groups()
                ssm_seqpos = int(groups[0])
                ssm_letter = groups[1]

                class_row = class_df[class_df['seqpos'] == ssm_seqpos].iloc[0]
                is_interface = class_row['is_interface_core'] or class_row['is_interface_boundary']
                if ( not is_interface ):
                    assert(class_row['is_monomer_core'] or class_row['is_monomer_boundary'] or class_row['is_monomer_surface'] )

                pae_slice = pae_matrix_dict[tag][ssm_seqpos-1]
                distance_from_us = cb_dist_dict[ssm_parent][ssm_seqpos-1]

                monomer_size = len(native_seq)
                monomer_distance_from_us = distance_from_us[:monomer_size]
                target_distance_from_us = distance_from_us[monomer_size:]
                argsort_monomer = np.argsort(monomer_distance_from_us)
                argsort_target = np.argsort(target_distance_from_us) + monomer_size

                if ( is_interface ):
                    nearest5_target_pae_from_mutation = np.mean(pae_slice[argsort_target[:6]])
                    special_region_pae = nearest5_target_pae_from_mutation * 0.125
                else:
                    binder_pae_from_mutation = np.mean(pae_slice[argsort_monomer[:6]])
                    special_region_pae = binder_pae_from_mutation + 0.25


                d = {}
                d['description'] = tag 
                d['special_region_pae'] = special_region_pae

                these_records.append(d)

        pae_matrix_df = pd.DataFrame(these_records)

        old_size = len(df)
        df = df.merge(pae_matrix_df, 'inner', 'description')

        assert(len(df) == old_size)



    if ( len( args.position_class_files) > 0 ):
        class_df = class_dict[ssm_parent].copy()
        class_df['is_interface'] = class_df['is_interface_core'] | class_df['is_interface_boundary']
        class_df['ssm_seqpos'] = class_df['seqpos']
        class_df = class_df[['is_interface', 'ssm_seqpos']]

        old_size = len(df)
        df = df.merge(class_df, 'inner', 'ssm_seqpos')
        assert(len(df) == old_size)

        df.loc[df['is_interface'], 'total_score_monomer'] = 0



    # Assign the native letters. If the assertion throws, it's because your numbering is messed up
    df['native'] = False

    all_natives_present = True
        
    for i in range(len(native_seq)):
        seqpos = i+1
        letter = native_seq[i]
        mask = (df['ssm_letter'] == letter) & (df['ssm_seqpos'] == seqpos)

        if ( args.partial_ssm ):
            if ( mask.sum() == 0 ):
                if ( (df['ssm_seqpos'] == seqpos).sum() != 0):
                    # import IPython
                    # IPython.embed()
                    all_natives_present = False
                    print("Error! Native not present for position %i but mutations are"%(seqpos))
                    df = df[df['ssm_seqpos'] != seqpos].copy()
                # assert((df['ssm_seqpos'] == seqpos).sum() == 0)
                continue


        assert(mask.sum() == 1)
        df.loc[mask, 'native'] = True

    # if ( not all_natives_present ):
    #     continue

    native_df = df[df['native']].copy()
    native_df['native_ddg_no_repack'] = native_df['ddg_no_repack']
    native_df['native_total_score_monomer'] = native_df['total_score_monomer']
    native_df['native_pae_interaction'] = native_df['pae_interaction']
    extra = []
    if ( "aditya_score" in list(df) ):
        native_df['native_aditya_score'] = native_df['aditya_score']
        extra.append("native_aditya_score")

    if ( "native_special_region_pae" in list(df) ):
        extra.append("native_special_region_pae")
        native_df['native_special_region_pae'] = native_df['special_region_pae']


    native_df = native_df[extra + ['native_ddg_no_repack', 'native_total_score_monomer', 
                                                    'native_pae_interaction', 'ssm_seqpos', 'ssm_parent']]

    old_size = len(df)
    df = df.merge(native_df, 'inner', ['ssm_seqpos', 'ssm_parent'])
    assert(len(df) == old_size)


    df['delta_ddg_no_repack'] = df['ddg_no_repack'] - df['native_ddg_no_repack']
    df['delta_total_score_monomer'] = df['total_score_monomer'] - df['native_total_score_monomer']
    df['delta_pae_interaction'] = df['pae_interaction'] - df['native_pae_interaction']

    if ( "aditya_score" in list(df) ):
        df['delta_aditya_score'] = df['aditya_score'] - df['native_aditya_score']

    if ( "native_special_region_pae" in list(df) ):
        df['delta_native_special_region_pae'] = df['native_special_region_pae'] - df['native_native_special_region_pae']




    pae_scaler = 1 #10 / df['native_pae_interaction']

    # comment out pae_interaction
    df['delta_E'] = df['delta_pae_interaction'] * 2.5 * pae_scaler + df['delta_ddg_no_repack'] + df['delta_total_score_monomer'].clip(0, None)

    if ( "aditya_score" in list(df) ):
        df['delta_E'] += - 0.5 * df['delta_aditya_score']

    if ( "native_special_region_pae" in list(df) ):
        df['delta_E'] += 2.5 * pae_scaler * df['native_native_special_region_pae']



    # Take the median of the native values so they're all the same
    # is_native = df['native']
    # avg_E = df[is_native]['my_E'].median()
    # df.loc[is_native, 'my_E'] = avg_E

    # Store the values
    df['kd_ub'] = np.exp( (df['delta_E']) / 0.6)
    df['kd_lb'] = df['kd_ub']
    df['low_conf'] = False
    df['lowest_conc'] = df['kd_lb'].min() / 10
    df['highest_conc'] = df['kd_ub'].max() * 10



    # In[6]:


    # The next script is expecting the native to be called __native
    one_row = df[df['native']].iloc[:1].copy()
    one_row['description'] = re.sub("__.*", "__native", one_row['description'].iloc[0])

    to_dump = pd.concat((one_row, df))
    final_dfs.append(to_dump)


to_dump = pd.concat(final_dfs)    

final_extra = [x.replace("native", "delta") for x in extra]


to_dump = to_dump[['delta_E', 'delta_pae_interaction', 'delta_ddg_no_repack', 'delta_total_score_monomer'] + final_extra + ['native',
                                        'kd_lb', 'kd_ub', 'low_conf', 'lowest_conc', 'highest_conc',
                                        'ssm_parent', 'ssm_seqpos', 'ssm_letter', 'description']]

to_dump.to_csv(args.output_file, index=None, sep=" ")


# In[ ]:





# In[ ]:





# In[ ]:




