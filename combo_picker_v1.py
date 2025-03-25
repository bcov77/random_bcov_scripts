#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import argparse
import warnings
import re
from collections import defaultdict
import itertools

sys.path.append("/home/bcov/sc/random/npose")
from importlib import reload
import voxel_array
reload(voxel_array)
import npose_util
reload(npose_util)
import npose_util as nu
import npose_util_pyrosetta
reload(npose_util_pyrosetta)
import npose_util_pyrosetta as nup

# from numba import njit
# import numba





parser = argparse.ArgumentParser(description="")
parser.add_argument("parent_pdb", type=str, help="The pdb file you wish to perform this analysis on.")
parser.add_argument("--input_affinity_file", type=str, default="", help="The output from estimate_affinity_from_ngs.py. Your SSM mutations must be labeled"
                                                            +" such that they are <parent_name>__<seqpos>__<letter>. And the parent design should"
                                                            +" be labeled <parent_name>__native. See sorting_ngs_data/*_ssm/pooled_counts.list"
                                                            +" for examples.")
parser.add_argument("--input_energies_file", type=str, default="", help="The output from estimate_affinity_from_ngs.py. Your SSM mutations must be labeled"
                                                            +" such that they are <parent_name>__<seqpos>__<letter>. And the parent design should"
                                                            +" be labeled <parent_name>__native. See sorting_ngs_data/*_ssm/pooled_counts.list"
                                                            +" for examples.")

parser.add_argument("--low_conf_kcal_above_native", type=float, default=0.99, help="If we don't know the kd, how much worse do we say it is?")

parser.add_argument("--dont_cache_sap", action="store_true", help="Don't stash the sap values")

parser.add_argument("--kcal_per_better_sap", type=float, default="0.25", help="If sap gets better, how many kcal per sap")
parser.add_argument("--kcal_per_worse_sap", type=float, default="1.5", help="If sap gets worse, how many kcal per sap")
parser.add_argument("--dsap_poison_threshold", type=float, default="2", help="If sap increases by more than this and the adjusted energy is"
                                                                             " > 0. Consider this to be a poison mutation.")
parser.add_argument("--poison_even_if_better", action="store_true", help="Even if the energy is < 0, big dsap residues are still poison")
parser.add_argument("--dna_sequence", type=str, default="", help="DNA sequence or we'll just pick the best codon at each position")

parser.add_argument("--interface_core_weight", type=float, default=1, help="Weight factor for seqpos")
parser.add_argument("--interface_boundary_weight", type=float, default=0.8, help="Weight factor for seqpos")
parser.add_argument("--monomer_core_weight", type=float, default=0.8, help="Weight factor for seqpos")
parser.add_argument("--monomer_boundary_weight", type=float, default=0.5, help="Weight factor for seqpos")
parser.add_argument("--monomer_surface_weight", type=float, default=0.3, help="Weight factor for seqpos")
parser.add_argument("--loop_weight", type=float, default=1.0, help="Weight factor for seqpos loop. Overrides other weight factors")
parser.add_argument("--allow_interface_by_distance", action="store_true", help="Weight factor for seqpos loop. Overrides other weight factors")


parser.add_argument("--min_distance_to_target", type=float, default=0, help="Distance where weight factor is 1")
parser.add_argument("--max_distance_to_target", type=float, default=np.nan, help="Distance where weight factor goes to 0")

parser.add_argument("--within_category_decay", type=float, default=0.25, help="Per pos score: Within each category, first aa gets full score,"
                                                                              " second aa gets score*within_category_decay"
                                                                              " third aa gets score*within_category_decay**2")

parser.add_argument("--between_category_decay", type=float, default=0.5, help="Per pos score: First category gets full score,"
                                                                              " second category gets score*between_category_decay"
                                                                              " third category gets score*between_category_decay**2")

parser.add_argument("--extraneous_aa_threshold", type=float, default=1, help="Per pos score: For each aa with score above this included,"
                                                                             " give extraneous_aa_penalty to codon")

parser.add_argument("--extraneous_aa_penalty", type=float, default=1, help="Per pos score: Score added for each aa above extraneous_aa_threshold")

parser.add_argument("--cys_is_not_extraneous", action="store_true", help="By default, non-native CYS gets the --extraneous_aa_penalty.")

parser.add_argument("--max_diversity", type=float, default=1e7, help="Max diversity that will be output.")

parser.add_argument("--focus_positions", type=str, default="", help="These comma-separated positions get their scores multiplied by --focus_multiplier")
parser.add_argument("--focus_multiplier", type=float, default=5, help="Score multiplier for focus positions")

parser.add_argument("--internal_inv_score_resl", type=float, default=100, help="Score resolution for calculations. Probably don't change this."
                                                                            " Only change if you're doing something weird and script is slow. Try 10")

parser.add_argument("--user_override", type=str, default="", help="Comma separated list of positions and mutations that get assigned --user_override_score."
                                                                " Example. 4DEQ,7NEV,15LR ")
parser.add_argument("--user_override_score", type=float, default=-10, help="Score assigned for user_override mutations")
parser.add_argument("--user_override_is_bonus", action="store_true", help="Add override score as bonus instead of rewrite")
parser.add_argument("--user_override_is_clip_bonus", action="store_true", help="Add score like bonus unless >0 in which case make it bonus")

parser.add_argument("--sketch_kd_threshold", type=float, default=1, help="Mutations with sketch_kd worse than this will have their energies set to 0")


parser.add_argument("--global_shift", type=float, default=0, help="Shift all scores by this")

args = parser.parse_args(sys.argv[1:])

print("Loading pyrosetta")

from pyrosetta import *
from pyrosetta.rosetta import *

# setting the dunbrack prob to 0.6 is severe. We only need to build non-clashing rotamers though for sap
init("-mute all -beta_nov16 -dunbrack_prob_buried 0.6 -dunbrack_prob_nonburied 0.6 -dunbrack_prob_buried_semi 0.6 -dunbrack_prob_nonburied_semi 0.6")


def fix_scorefxn(sfxn):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    sfxn.set_energy_method_options(opts)
scorefxn = get_fa_scorefxn()
fix_scorefxn(scorefxn)

focus_positions = []
if ( args.focus_positions != "" ):
    focus_positions = [int(x) for x in args.focus_positions.split(",")]


sfxn_fast = core.scoring.ScoreFunctionFactory.create_score_function("none")
sfxn_fast.set_weight(core.scoring.fa_rep, 0.55)
sfxn_fast.set_weight(core.scoring.fa_atr, 1.00)
sfxn_fast.set_weight(core.scoring.fa_sol, 1.00)
sfxn_fast.set_weight(core.scoring.fa_dun_rot, 0.76)
sfxn_fast.set_weight(core.scoring.fa_dun_dev, 0.76)
sfxn_fast.set_weight(core.scoring.fa_dun_semi, 0.76)

df = None

if ( args.input_affinity_file != "" ):
      if ( args.input_energies_file != "" ):
            print("You only need one of: --input_affinity_file --input_energies_file")
            assert(False)

      kd_df = pd.read_csv(args.input_affinity_file, sep="\s+")

      assert("kd_lb" in list(kd_df))
      assert("kd_ub" in list(kd_df))
      assert("lowest_conc" in list(kd_df))
      assert("highest_conc" in list(kd_df))
      assert("low_conf" in list(kd_df))
      assert("description" in list(kd_df))

      kd_df['description'] = kd_df['description'].str.replace("([^_])_([0-9]+)_([A-Z])$", r"\1__\2__\3")

      kd_df['kd_center'] = np.sqrt( kd_df['kd_lb'].clip(kd_df['lowest_conc']/10, kd_df['highest_conc']*1e8)
                                  * kd_df['kd_ub'].clip(kd_df['lowest_conc']/10, kd_df['highest_conc']*1e8)
                                  )

      kd_df['energy'] = 0.6*np.log( kd_df['kd_center'] ) 

      kd_df.loc[kd_df['low_conf'], 'energy'] = np.nan




      if args.sketch_kd_threshold < 1:
        energy_df = kd_df[['energy', 'sketch_kd_level', 'description']].copy()
      else:
        energy_df = kd_df[['energy', 'description']].copy()

if ( args.input_energies_file != "" ):
      energy_df = pd.read_csv(args.input_energies_file, sep="\s+")
      if ( "energy" not in list(energy_df) or "description" not in list(energy_df)):
            print("--input_energies_file should have 2 columns. 'energy' and 'description'")
            assert(False)
      energy_df = energy_df[['energy', 'description']].copy()

energy_df = energy_df.drop_duplicates('description')


abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))


def get_consensus(letters):
    counts = defaultdict(lambda : 0, {})
    for letter in letters:
        counts[letter] += 1

    maxx_letter = 0
    maxx = 0
    for key in counts:
        if ( counts[key] > maxx ):
            maxx = counts[key]
            maxx_letter = key
    return maxx_letter


# this is 1 indexed with the start and end with loops converted to nearby dssp
# and HHHHHH turns identified
def better_dssp2(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()
    the_dssp = list(the_dssp)

    n_consensus = get_consensus(the_dssp[3:6])

    the_dssp[1] = n_consensus
    the_dssp[2] = n_consensus
    the_dssp[3] = n_consensus
    the_dssp[4] = n_consensus
    the_dssp[5] = n_consensus

    c_consensus = get_consensus(the_dssp[-5:-2])

    the_dssp[-1] = c_consensus
    the_dssp[-2] = c_consensus
    the_dssp[-3] = c_consensus
    the_dssp[-4] = c_consensus
    the_dssp[-5] = c_consensus

    the_dssp = "".join(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" or abego == "E" and seqpos > 5 and seqpos < len(the_dssp)-5 ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp

def load_pose_data(pdbs):

    chainA = core.select.residue_selector.ChainSelector("A")
    chainB = core.select.residue_selector.ChainSelector("B")
    interface_on_A = core.select.residue_selector.NeighborhoodResidueSelector(chainB, 10.0, False)
    interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 10.0, False)
    interface_by_distance = core.select.residue_selector.OrResidueSelector(interface_on_A, interface_on_B)
    interface_by_vector = core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_A, interface_on_B)
    interface_by_vector.cb_dist_cut(11)
    interface_by_vector.cb_dist_cut(5.5)
    interface_by_vector.vector_angle_cut(75)
    interface_by_vector.vector_dist_cut(9)

    pose_dats = []
    sequences = {}
    poses = []

    # We're loading the poses to figure out their sequence as well as to determine
    # The core/boundary/surface and interface/monomer
    for pdb in pdbs:
        print("    " + pdb)

        name = os.path.basename(pdb)
        if ( name.endswith(".gz") ):
            name = name[:-len(".gz")]
        if ( name.endswith(".pdb") ):
            name = name[:-len(".pdb")]

        pose = pose_from_file(pdb)
        monomer = pose.split_by_chain()[1]
        sequence = monomer.sequence()
        dssp = "x" + better_dssp2(monomer)

        sc_neighbors = core.select.util.SelectResiduesByLayer()
        sc_neighbors.use_sidechain_neighbors( True )
        sc_neighbors.compute(pose, "", True)

        atomic_depth = core.scoring.atomic_depth.AtomicDepth( pose, 2.3, False, 0.5 )
        atomic_depth_monomer = core.scoring.atomic_depth.AtomicDepth( monomer, 2.3, False, 0.5 )
        type_set = pose.residue(1).type().atom_type_set()

        probe_size = 2.8
        per_atom_sasa = core.id.AtomID_Map_double_t()
        rsd_sasa = utility.vector1_double()
        core.scoring.calc_per_atom_sasa(pose, per_atom_sasa, rsd_sasa, 2.8, False)

        scorefxn(pose)
        scorefxn(monomer)

        interface_subset = interface_by_vector.apply(pose)
        interface_by_dist_subset = interface_by_distance.apply(pose)

        for seqpos in range(1, monomer.size()+1):

            data = {"ssm_parent":name, "ssm_seqpos":seqpos}

            data['sc_neighbors'] = sc_neighbors.rsd_sasa(seqpos)
            data['is_loop'] = dssp[seqpos] == "L"
            data['by_vector'] = interface_subset[seqpos]
            data['by_dist'] = interface_by_dist_subset

            res = pose.residue(seqpos)
            monomer_res = monomer.residue(seqpos)
            data['depth'] = atomic_depth.calcdepth(res.atom(res.nbr_atom()), type_set)
            data['depth_monomer'] = atomic_depth_monomer.calcdepth(monomer_res.atom(monomer_res.nbr_atom()), type_set)

            data['ddg'] = 2*(pose.energies().residue_total_energy(seqpos) - monomer.energies().residue_total_energy(seqpos))

            sc_sasa = 0
            for i in range(res.first_sidechain_atom(), res.nheavyatoms()+1):
                sc_sasa += per_atom_sasa(seqpos, i)
                for j in range(res.attached_H_begin(i), res.attached_H_end(i)+1):
                    sc_sasa += per_atom_sasa(seqpos, j)

            data['sc_sasa'] = sc_sasa

            pose_dats.append(data)

        sequences[name] = sequence

        poses.append(pose)

    pose_df = pd.DataFrame(pose_dats)

    # These should be option flags. But this calculation is so complicated that unless you're looking at the code, you're not going
    #  to specify this stuff correctly.
    # So either change the hardcoded stuff, or make the option flags yourself. (Or even change the calculation here)
    has_ddg = np.abs(pose_df['ddg']) > 1
    if args.allow_interface_by_distance:
        has_ddg |= pose_df['by_dist']
    is_core = ((pose_df['sc_neighbors'] > 5.2) | (pose_df['depth'] - pose_df['depth_monomer'] > 1 )) & (pose_df['sc_sasa'] < 5)

    pose_df['is_interface_core'] = has_ddg & is_core
    pose_df['is_interface_boundary'] = ( has_ddg | pose_df['by_vector'] ) & ~pose_df['is_interface_core']
    pose_df['is_monomer_core'] = ~has_ddg & (pose_df['sc_neighbors'] >= 5.2) & ~pose_df['is_interface_boundary']
    pose_df['is_monomer_boundary'] = ~has_ddg & (pose_df['sc_neighbors'] < 5.2) & (pose_df['sc_neighbors'] >= 2.0) & ~pose_df['by_vector']
    pose_df['is_monomer_surface'] = ~has_ddg & (pose_df['sc_neighbors'] < 2.0) & ~pose_df['by_vector']

    # Make sure all positions are in exactly 1 category
    assert( np.all( pose_df[['is_interface_core', 'is_interface_boundary', 'is_monomer_core', 'is_monomer_boundary',
                            'is_monomer_surface']].sum(axis=1) == 1) )

    pose_df['ssm_seqpos'] = pose_df['ssm_seqpos'].astype(str)

    return sequences, pose_df, poses


def break_disulfide(pose, seqpos):
    if ( pose.residue(seqpos).name1() != "C" ): return
    try:
        other = core.conformation.get_disulf_partner(pose.conformation(), seqpos)
    except:
        return

    try:
        core.conformation.break_disulfide(pose.conformation(), seqpos, other)
    except:
        pass


def get_sap_scores(pose, name):
    fname = name + ".sap"

    og_pose = pose.split_by_chain()[1]
    sequence = og_pose.sequence()

    true_sel = core.select.residue_selector.TrueResidueSelector()


    if ( not args.dont_cache_sap and os.path.exists(fname) ):
          sap_df = pd.read_csv(fname, sep="\s+")
          return sap_df

    records = []

    for pos in range(1, og_pose.conformation().chain_end(1)+1):
        pose = og_pose.clone()
        if ( pose.residue(pos).name1() == "C" ):
            break_disulfide( pose, pos )
        for name1 in "ACDEFGHIKLMNPQRSTVWY":
  
                
            protocols.toolbox.pose_manipulation.repack_this_residue(pos, pose, sfxn_fast, False, name1)
  
            sap = core.pack.guidance_scoreterms.sap.calculate_sap(pose, true_sel, true_sel, true_sel)
  
            description = name + "__" + str(pos) + "__" + name1

            records.append({'sap_score':sap, 'description':description})

    sap_df = pd.DataFrame(records)
    if ( not args.dont_cache_sap ):
        sap_df.to_csv(fname, sep=" ", index=None)

    return sap_df


name = os.path.basename(args.parent_pdb).replace(".gz", "").replace(".pdb", "")

print("Loading pose")
sequences, pose_df, poses = load_pose_data([args.parent_pdb])
sequence = sequences[name]
pose = poses[0]

sys.stdout.write("Calculating saps...")
sys.stdout.flush()
sap_df = get_sap_scores(pose, name)
sys.stdout.write("Done\n")
sys.stdout.flush()



monomer_size = len(sequence)
distance_based_factors = np.ones(monomer_size, np.float)

if ( not np.isnan(args.max_distance_to_target) ):
    npose = nup.npose_from_pose(pose)
    cas = nu.extract_atoms(npose, [nu.CA])[:,:3]

    monomer_cas = cas[:monomer_size]
    target_cas = cas[monomer_size:]
    closest_target = np.linalg.norm( monomer_cas[:,None] - target_cas[None,:], axis=-1 ).min(axis=-1)
    assert(len(closest_target) == monomer_size)

    span = args.max_distance_to_target - args.min_distance_to_target

    distance_based_factors = ( 1 - (closest_target - args.min_distance_to_target)/span).clip(0, 1)

dist_records = []
for seqpos0, factor in enumerate(distance_based_factors):
    dist_records.append({"ssm_seqpos":str(seqpos0+1), "dist_factor":factor})
dist_df = pd.DataFrame(dist_records)


print("Finalizing scores")

# cases:
# 1. __native exists and nothing else             -- fill __native to rest
# 2. __native exists and some/all others exist    -- fill __native to missing, mark rest
# 3. __native missing and some/all others exist   -- fill first to missing, mark rest

# steps:
# A. Identify fill row
# B. Filling missing natives
# C. Mark natives


def nat_name(sequence, name, position):
    letter = sequence[position-1]
    return name + "__" + str(position) + "__" + letter

native_row = None

energy_df['native'] = False

maybe_native_rows = energy_df[energy_df['description'].str.contains("^" + name + "_+native$")]
if ( len(maybe_native_rows) > 0 ):
    if ( len(maybe_native_rows) > 1 ):
        print("Error: You have multiple __native rows in your energies file?")
        assert(False)
    native_row = maybe_native_rows.iloc[0].copy()
else:
    for seqpos in range(1, len(sequence)+1):
        maybe = energy_df[energy_df['description'] == nat_name(sequence, name, seqpos)]
        if ( len(maybe) == 0 ):
            continue
        native_row = maybe.iloc[0]
        break
    if ( native_row is None ):
        print("Error: You don't have any mutations corresponding to native?")
        assert(False)

assert( not native_row is None )

native_row['native'] = True


new_rows = []
for seqpos in range(1, len(sequence)+1):
    nat_desc = nat_name(sequence, name, seqpos)
    maybe = energy_df[energy_df['description'] == nat_desc]
    if ( len(maybe) == 0 ):
        new_row = native_row.copy()
        new_row['description'] = nat_desc
        new_rows.append(new_row)
    else:
        energy_df.loc[maybe.index[0], 'native'] = True

extra_energy_df = pd.DataFrame(new_rows)

energy_df = pd.concat((energy_df, extra_energy_df)).reset_index(drop=True)
if args.sketch_kd_threshold < 1:
    energy_df.loc[energy_df['sketch_kd_level'] < args.sketch_kd_threshold, 'energy'] = np.nan
energy_df['energy'] = energy_df['energy'].fillna(native_row['energy'] + args.low_conf_kcal_above_native)


df = energy_df.merge(sap_df, 'right', 'description')

parts = df['description'].str.extract(r"(.*[^_])_+([0-9]+)_+([A-Z])$")
df['ssm_parent'] = parts[0]
df['ssm_seqpos'] = parts[1].astype(str)
df['ssm_letter'] = parts[2]

if df['energy'].isnull().any():

    missing_letters = df[df['energy'].isnull()]['ssm_letter'].unique()
    if "".join(list(missing_letters)) != 'C':
        print("Missing rows in dataframe other than C")
        assert(False)

    df['energy'] = df['energy'].fillna(10000)
    df['native'] = df['native'].fillna(False)


assert(len(df) == len(sap_df))


df = df.merge(pose_df, 'inner', ['ssm_seqpos', 'ssm_parent'])

assert(len(df) == len(sap_df))



parent_df = df[df['native']][['energy', 'sap_score', 'ssm_seqpos', 'ssm_parent']].copy()
parent_df['parent_energy'] = parent_df['energy']
parent_df['parent_sap_score'] = parent_df['sap_score']
parent_df = parent_df[['parent_energy', 'parent_sap_score', 'ssm_seqpos', 'ssm_parent']]


df = df.merge(parent_df, 'inner', ['ssm_seqpos', 'ssm_parent'])

assert(len(df) == len(sap_df))

df['delta_energy'] = df['energy'] - df['parent_energy']
df['delta_sap'] = df['sap_score'] - df['parent_sap_score']


sap_up = df['delta_sap'] > 0

df.loc[sap_up, 'delta_energy'] += df[sap_up]['delta_sap'] * args.kcal_per_worse_sap
df.loc[~sap_up, 'delta_energy'] += df[~sap_up]['delta_sap'] * args.kcal_per_better_sap


df['position_weight_factor'] = 1
df.loc[df['is_interface_core'], 'position_weight_factor'] = args.interface_core_weight
df.loc[df['is_interface_boundary'], 'position_weight_factor'] = args.interface_core_weight
df.loc[df['is_monomer_core'], 'position_weight_factor'] = args.monomer_core_weight
df.loc[df['is_monomer_boundary'], 'position_weight_factor'] = args.monomer_boundary_weight
df.loc[df['is_monomer_surface'], 'position_weight_factor'] = args.monomer_surface_weight
df.loc[df['is_loop'], 'position_weight_factor'] = args.loop_weight                          # keep this last because it's not mutually exclusive

old_size = len(df)
df = df.merge(dist_df, 'inner', 'ssm_seqpos')
assert(len(df) == old_size)
df['position_weight_factor'] *= df['dist_factor']
df['position_weight_factor'] = df['position_weight_factor'].clip(1e-10, None)


if ( args.user_override != "" ):
    for group in args.user_override.split(","):
        groups = re.match("([0-9]+)([A-Z]+)", group)
        if ( groups is None ):
            print("Bad user_override specification", group)
            assert(False)
        seqpos = int(groups.group(1))
        letters = groups.group(2)
        for letter in letters:
            mask = (df['ssm_seqpos'].astype(int) == seqpos) & (df['ssm_letter'] == letter)
            
            if ( args.user_override_is_clip_bonus ):
                df.loc[mask, 'delta_energy'] = df.loc[mask, 'delta_energy'].clip(None, 0) + args.user_override_score / df[mask]['position_weight_factor']
            elif ( args.user_override_is_bonus ):
                df.loc[mask, 'delta_energy'] += args.user_override_score / df[mask]['position_weight_factor']
            else:
                # we're cancelling out the region-specific weight here
                df.loc[mask, 'delta_energy'] = args.user_override_score / df[mask]['position_weight_factor']


for focus_position in focus_positions:
    df.loc[df['ssm_seqpos'].astype(int) == focus_position, 'position_weight_factor'] *= args.focus_multiplier

df['delta_energy'] += args.global_shift


if ( args.poison_even_if_better ):
    df['is_poison'] = (df['delta_sap'] > args.dsap_poison_threshold)
else:
    df['is_poison'] = (df['delta_sap'] > args.dsap_poison_threshold) & (df['delta_energy'] > 0)

df.loc[df['is_poison'], 'delta_energy'] = np.inf

dump_df = df[['delta_energy', 'native', 'description']].copy()
dump_df['kd_ub'] = np.exp(dump_df['delta_energy']/0.6)
dump_df['kd_lb'] = np.exp(dump_df['delta_energy']/0.6)
dump_df['low_conf'] = False
dump_df['avid_doesnt_agree'] = False
dump_df['lowest_conc'] = dump_df['kd_ub'].min()
dump_df['highest_conc'] = 1000000
dump_df.loc[np.isinf(dump_df['delta_energy']), 'low_conf'] = True
dump_native = dump_df[dump_df['native']].iloc[:1]
dump_native['description'] = name + "__native"
dump_df = pd.concat((dump_df, dump_native))
dump_df.to_csv(name + "_heatmap.sc", sep=" ", index=None, na_rep='Nan')





codon_to_aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

aa_to_codons = defaultdict(list)
for codon in codon_to_aa:
    aa = codon_to_aa[codon]
    aa_to_codons[aa].append(codon)

codon_frequency_yeast = {
"TTT": 0.59,"ATT": 0.46,"TCT": 0.26,"ACT": 0.35,
"TTC": 0.41,"ATC": 0.26,"TCC": 0.16,"ACC": 0.22,
"TTA": 0.28,"ATA": 0.27,"TCA": 0.21,"ACA": 0.3,
"TTG": 0.29,"ATG": 1,"TCG": 0.1,"ACG": 0.13,
"TAT": 0.56,"AAT": 0.59,"TGT": 0.63,"AGT": 0.16,
"TAC": 0.44,"AAC": 0.41,"TGC": 0.37,"AGC": 0.11,
"TAA": 0.48,"AAA": 0.58,"TGA": 0.29,"AGA": 0.48,
"TAG": 0.24,"AAG": 0.42,"TGG": 1,"AGG": 0.21,
"CTT": 0.13,"GTT": 0.39,"CCT": 0.31,"GCT": 0.38,
"CTC": 0.06,"GTC": 0.21,"CCC": 0.15,"GCC": 0.22,
"CTA": 0.14,"GTA": 0.21,"CCA": 0.41,"GCA": 0.29,
"CTG": 0.11,"GTG": 0.19,"CCG": 0.12,"GCG": 0.11,
"CAT": 0.64,"GAT": 0.65,"CGT": 0.15,"GGT": 0.47,
"CAC": 0.36,"GAC": 0.35,"CGC": 0.06,"GGC": 0.19,
"CAA": 0.69,"GAA": 0.71,"CGA": 0.07,"GGA": 0.22,
"CAG": 0.31,"GAG": 0.29,"CGG": 0.04,"GGG": 0.12,
}
aa_best_to_worst_codon = {}
for aa in aa_to_codons:
    codons = aa_to_codons[aa]
    freqs = [codon_frequency_yeast[x] for x in codons]
    arg_sorted = np.argsort(-np.array(freqs))

    best_to_worst = [codons[arg_sorted[i]] for i in range(len(arg_sorted))]
    aa_best_to_worst_codon[aa] = best_to_worst


# Just pick the best codon at each position, but don't repeat 2 in a row if you can help it
def sketchy_codon_optimizer(sequence):
    dna_sequence = ""
    last_codon = ""
    for letter in sequence:
        codons = aa_best_to_worst_codon[letter]
        if ( len(codons) == 1 ):
            codon = codons[0]
        else:
            codon = codons[0]
            if ( codon == last_codon ):
                codon = codons[1]
        last_codon = codon
        dna_sequence += codon
    return dna_sequence




categories = {
    'polar':'HQNST',
    'aromatic':'FWY',
    'packing':'ILMV',
    'folding':'AGP',
    'negative':'DE',
    'positive':'RK'
}

all_codons = {
'A':['A'],
'T':['T'],
'G':['G'],
'C':['C'],
'R':['A','G'],
'Y':['C','T'],
'M':['A','C'],
'K':['G','T'],
'S':['G','C'],
'W':['A','T'],
'H':['A','C','T'],
'B':['G','C','T'],
'V':['A','C','G'],
'D':['A','G','T'],
'N':['A','C','G','T']
}

dna_sequence = args.dna_sequence
if ( dna_sequence == "" ):
    print("==== No DNA sequence provided, using sketchy_codon_optimizer() ====")
    dna_sequence = sketchy_codon_optimizer(sequence)


# print("Preparing dynamic programming")


# score per position:
#  sum of max( delta_energy ) from each category
# +inf if poison allowed


aa_order = 'CPGAVIMLFYWSTNQDERKH*'
aa_to_offset = {}
for iletter, letter in enumerate(aa_order):
    aa_to_offset[letter] = iletter


category_masks = []
for cat, letters in categories.items():
    mask = np.zeros(len(aa_order), bool)
    for letter in letters:
        mask[aa_to_offset[letter]] = True
    category_masks.append(mask)

# figure out which aas a given degenerate codon is going to produce
degen_codon_to_mask_and_diveristy = {}
for codon1 in all_codons:
    codons1 = all_codons[codon1]
    for codon2 in all_codons:
        codons2 = all_codons[codon2]
        for codon3 in all_codons:
            codons3 = all_codons[codon3]

            diversity = len(codons1)*len(codons2)*len(codons3)
            sum_mask = np.zeros(len(aa_order), int)
            for c1, c2, c3 in itertools.product(codons1, codons2, codons3):
                letter = codon_to_aa[c1+c2+c3]
                sum_mask[aa_to_offset[letter]] += 1

            mask = sum_mask > 0
            assert(sum_mask.sum() == diversity)

            degen_codon_to_mask_and_diveristy[codon1 + codon2 + codon3] = (mask, sum_mask, diversity)

all_degen_order = list(degen_codon_to_mask_and_diveristy)
base_degen = len(all_degen_order)

score_table = np.zeros((len(sequence), len(aa_order)))
score_table[:] = np.nan

aa_offsets = df['ssm_letter'].map(aa_to_offset).values
seqpos_offsets = df['ssm_seqpos'].astype(int).values - 1
denergies = df['delta_energy'].values * df['position_weight_factor'].values

# Put the aa energies into the score table
score_table[tuple((seqpos_offsets, aa_offsets ))] = denergies

# Mark all stop codons as extraneous
score_table[:, aa_to_offset['*']] = args.extraneous_aa_threshold

# Mark all non-native CYS as extraneous
if ( not args.cys_is_not_extraneous ):
    for seqpos0 in range(len(sequence)):
        if ( sequence[seqpos0] != "C" ):
            score_table[seqpos0, aa_to_offset['C']] = args.extraneous_aa_threshold

assert( not np.any( np.isnan(score_table)))


total_combinations = 1

# Get all valid degenerate codon choices at each position
per_position_choices = []
for seqpos0 in range(len(sequence)):

    # loop through literally every single denerate codon at this position
    this_position_choices = []
    for idegen, degen_codon in enumerate(all_degen_order):

        # allowed aa mask, and total codons added
        mask, sum_mask, diversity = degen_codon_to_mask_and_diveristy[degen_codon]

        native_aa = sequence[seqpos0]
        native_aa_offset = aa_to_offset[native_aa]

        # native aa must be present
        if ( not mask[native_aa_offset] ):
            continue

        # get aa scores for this seqpos
        our_scores = score_table[seqpos0].copy()

        # figure out extaneous before we set stuff to 0
        is_extraneous = our_scores >= args.extraneous_aa_threshold
        is_extraneous[native_aa_offset] = False

        # set all non-present aas to 0 and mark as not extraneous
        our_scores[~mask] = 0
        is_extraneous[~mask] = False

        # skip this codon if a poison is present
        if ( np.isinf(our_scores).any() ):
            continue

        # find score for each category
        cat_scores = []
        for cat_mask in category_masks:

            # all scores for this category
            sorted_cat_scores = sorted(our_scores[cat_mask])

            # apply --within_category_decay as we sum this category
            this_score = 0
            mult = 1
            for ind_score in sorted_cat_scores:
                if ( ind_score < 0 ):
                    this_score += ind_score * mult
                mult *= args.within_category_decay

            # only add to category scores if below 0
            if ( this_score < 0 ):
                cat_scores.append(this_score)

        # sum categories using --between_category_decay
        cat_scores = sorted(cat_scores)
        score = 0
        mult = 1
        for cat_score in cat_scores:
            score += cat_score*mult
            mult *= args.between_category_decay

        # for every extraneous aa this will produce (including synonymous sequences) apply the penalty
        score += sum_mask[is_extraneous].sum() * args.extraneous_aa_penalty

        # add this codon to the list of choices at this positon. Convert score to int
        if ( score < 0 ):
            this_position_choices.append((int(score*args.internal_inv_score_resl), diversity, idegen))

    per_position_choices.append(np.array(this_position_choices))
    total_combinations *= max(len(this_position_choices), 1)


# re-order the choices so that we only have the best choice at each diversity

all_diversities = []
all_scores = []
all_idegens = []

max_possible_score = 0

for seqpos0 in range(len(sequence)):

    these_choices = per_position_choices[seqpos0]

    diversities = []
    scores = []
    idegens = []

    diversities.append(1)
    scores.append(0)
    idegens.append(base_degen)

    best_score = 0
    if ( len(these_choices) > 0 ):
        ascending_diversity = np.argsort(these_choices[:,1])

        last_diversity = 1
        for i in ascending_diversity:
            score, diversity, idegen = these_choices[i]

            if ( score >= best_score ):
                continue

            if ( diversity == diversities[-1] ):
                diversities[-1] = diversity
                scores[-1] = score
                idegens[-1] = idegen
            else:
                diversities.append(diversity)
                scores.append(score)
                idegens.append(idegen)

            best_score = score

    diversities = np.array(diversities, int)
    scores = np.array(scores, int)
    idegens = np.array(idegens, int)

    # Expliticly changing the scores to positive here
    scores *= -1

    all_diversities.append(diversities)
    all_scores.append(scores)
    all_idegens.append(idegens)

    max_possible_score += best_score

max_possible_score = abs(max_possible_score)


print("Debug: max_possible_score", max_possible_score)

dp_diversity = np.zeros((len(sequence), max_possible_score+1), np.float64)
dp_idegen = np.zeros((len(sequence), max_possible_score+1), int)
dp_prev_idx = np.zeros((len(sequence), max_possible_score+1), int)

# max_value = np. 

dp_diversity[:] = np.inf
dp_idegen[:] -1
# dp_prev_idx[:] = -1


# idk, numba caused this to give different outputs. Don't feel like debugging plus that's super sketch
# print("Dynamic programming")
# @njit(cache=True)
def dynamic_programming( dp_diversity, dp_idegen, dp_prev_idx, all_diversities, all_scores, all_idegens ):

# seqpos loop of dynamic programming
    for seqpos0 in range(len(sequence)):

        diversities = all_diversities[seqpos0]
        scores = all_scores[seqpos0]
        idegens = all_idegens[seqpos0]


        for i_idegen in range(len(idegens)):
            diversity = diversities[i_idegen]
            score = scores[i_idegen]
            idegen = idegens[i_idegen]

            range_ub = max_possible_score+1
            if ( seqpos0 == 0 ):
                range_ub = 1

            for iprev_score in range(range_ub):
                if ( seqpos0 == 0 ):
                    prev_div = 1
                else:
                    prev_div = dp_diversity[seqpos0-1, iprev_score]

                if ( np.isinf( prev_div ) ):
                    continue

                new_div = diversity * prev_div
                new_score = score + iprev_score


                if ( new_div < dp_diversity[seqpos0, new_score] ):
                    dp_diversity[seqpos0, new_score] = new_div
                    dp_idegen[seqpos0, new_score] = idegen
                    dp_prev_idx[seqpos0, new_score] = iprev_score


# stupid numba update
# typed_all_diversities = numba.typed.List()
# [typed_all_diversities.append(x) for x in all_diversities]
# typed_all_scores = numba.typed.List()
# [typed_all_scores.append(x) for x in all_scores]
# typed_all_idegens = numba.typed.List()
# [typed_all_idegens.append(x) for x in all_idegens]

# dynamic_programming( dp_diversity, dp_idegen, dp_prev_idx, typed_all_diversities, typed_all_scores, typed_all_idegens )
dynamic_programming( dp_diversity, dp_idegen, dp_prev_idx, all_diversities, all_scores, all_idegens )


final_diversities = []
final_scores = []
best_diversity = np.inf

for score in reversed(range(max_possible_score+1)):

    this_diversity = dp_diversity[-1, score]

    if ( this_diversity < best_diversity ):
        final_diversities.append(this_diversity)
        final_scores.append(score)
        best_diversity = this_diversity

final_diversities = list(reversed(final_diversities))
final_scores = list(reversed(final_scores))

def traceback(final_score, do_error=False):
    reversed_idegens = []

    cur_score = final_score
    seqpos0 = len(sequence)-1

    while seqpos0 >= 0:
        reversed_idegens.append(dp_idegen[seqpos0, cur_score])

        cur_score = dp_prev_idx[seqpos0, cur_score]
        seqpos0 -= 1

    assert(len(reversed_idegens) == len(sequence))

    final_dna = ""
    for seqpos0, idegen in enumerate(reversed(reversed_idegens)):
        assert(idegen >= 0)
        this_degen = dna_sequence[seqpos0*3:(seqpos0+1)*3]
        if ( idegen < base_degen ):
            this_degen = all_degen_order[idegen]
        else:
            if ( do_error and codon_to_aa[this_degen] != sequence[seqpos0] ):
                print("DNA -- protein mismatch at position %i. AA=%s DNA=%s"%(seqpos0+1, sequence[seqpos0], codon_to_aa[this_degen]))
                assert( False )
        final_dna += this_degen

    return final_dna



def fancy_print_dna(dna):
    assert(len(dna) % 3 == 0)

    short_parts = []

    for seqpos0 in range(len(dna)//3):
        codon = dna[seqpos0*3:(seqpos0+1)*3]

        bases1 = all_codons[codon[0]]
        bases2 = all_codons[codon[1]]
        bases3 = all_codons[codon[2]]

        allowed_aas = set()

        for b1 in bases1:
            for b2 in bases2:
                for b3 in bases3:
                    letter = codon_to_aa[b1+b2+b3]
                    allowed_aas.add(letter)

        to_print = ""
        for letter in aa_order:
            if ( letter in allowed_aas ):
                to_print += letter

        if ( len(to_print) > 1 ):
            short_parts.append("%i"%(seqpos0+1) + to_print)

        print("%5i %s"%(seqpos0+1, to_print))

    return ",".join(short_parts)



our_choice = np.searchsorted(final_diversities, args.max_diversity)
if ( our_choice == len(final_diversities) ):
    our_choice -= 1


lb = max(0, our_choice-10)
ub = min(len(final_scores)-1, our_choice+10)

print("Nearby choices: ")
print(" diversity -- score (max delta kcal/mol)")
for choice in range(lb, ub+1):
    prefix = "   " if choice != our_choice else "***"
    print(prefix + "%6.1e -- %.2f"%(final_diversities[choice], final_scores[choice]/-args.internal_inv_score_resl))


print("Selected: ")

final_dna = traceback(final_scores[our_choice], do_error=True)
short_str = fancy_print_dna(final_dna)




def score_dna( dna, sequence ):
    assert(len(dna) % 3 == 0)

    overall_diversity = 1
    overall_frac_non_extraneous = 1
    score = 0

    short_parts = []

    for seqpos0 in range(len(dna)//3):
        codon = dna[seqpos0*3:(seqpos0+1)*3]

        mask, sum_mask, diversity = degen_codon_to_mask_and_diveristy[codon]

        overall_diversity *= diversity

        our_scores = score_table[seqpos0]
        is_extraneous = our_scores >= args.extraneous_aa_threshold
        is_extraneous[aa_to_offset[sequence[seqpos0]]] = False

        extraneous_sum = sum_mask[is_extraneous].sum()
        good_sum = diversity - extraneous_sum

        overall_frac_non_extraneous *= good_sum / diversity


        wanted_mask = mask & ~is_extraneous
        this_str = ""
        for il, letter in enumerate(aa_order):
            if ( wanted_mask[il] ):
                this_str += letter
        if ( len(this_str) > 1 ):
            short_parts.append("%i"%(seqpos0+1) + this_str)


    return 0, overall_frac_non_extraneous, overall_diversity,",".join(short_parts)



score, frac_non_extraneous, diversity, short_str = score_dna( final_dna, sequence )
print(short_str)

assert(np.isclose(diversity, final_diversities[our_choice]))

print("")
print("DNA Diversity: %6.1e"%diversity)
print("Percent desired aa: %.2f%%"%(frac_non_extraneous*100))
print("")

print("")
print(short_str)
print("")

print("DNA:")
print(final_dna)












