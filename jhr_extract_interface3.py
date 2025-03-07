#!/usr/bin/env python
from __future__ import division

# This program accepts arguments like this:

#./remove_superfluous_trp.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./remove_superfluous_trp.py -in:file:silent my.silent

import os
import sys
import math

import distutils.spawn
import os
import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools

from pyrosetta import *
from pyrosetta.rosetta import *

sys.path.append("/home/bcov/sc/random/npose")
import npose_util_pyrosetta as nup
import npose_util as nu

import numpy as np
from collections import defaultdict
import time
import argparse
import itertools
import subprocess
import time
import re

import sklearn.decomposition
import scipy.spatial.transform

# import pyRMSD.RMSDCalculator

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all"
    " -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc"
    )




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("-patchdock_res", type=str, default="66,68,73,106")
# parser.add_argument("-patchdock_res2", type=str, default="74,88,90,101")
parser.add_argument("-cp_per_res_cut", type=float, default=2)
parser.add_argument("-min_helix_cp", type=float, default=30)
parser.add_argument("-poly_letter", type=str, default="A")
parser.add_argument("-target_poly_letter", type=str, default="")
parser.add_argument("-min_helix_size", type=int, default=6)
parser.add_argument("-min_strand_size", type=int, default=3)
parser.add_argument("-align_target", type=str, default='')
parser.add_argument("-dump_all_helices", action="store_true")
parser.add_argument("-com_of_patch_res", action="store_true")
# parser.add_argument("-allow_strand", action="store_true")
parser.add_argument("-ideal_helix_pdb", type=str, default="/mnt/home/bcov/util/ideal_helix/ideal_helix_40.pdb")
parser.add_argument("-ideal_strand_pdb", type=str, default="/mnt/home/bcov/util/ideal_helix/ideal_strand_20.pdb")
parser.add_argument("-max_clashes_till_stop", type=int, default=4, help="n,ca,c,cb")


args = parser.parse_args(sys.argv[1:])
args.allow_strand = True

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")

align_pose = None
if args.align_target != '':
    align_pose = pose_from_file(args.align_target)

ideal_helix = pose_from_file(args.ideal_helix_pdb)
ideal_strand = pose_from_file(args.ideal_strand_pdb)

scorefxn = get_fa_scorefxn()
scorefxn_fa_atr = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_fa_atr.set_weight(core.scoring.fa_atr, 1)

scorefxn_none = core.scoring.ScoreFunctionFactory.create_score_function("none")


chainA = core.select.residue_selector.ChainSelector("A")
chainB = core.select.residue_selector.ChainSelector("B")
interface_on_A = core.select.residue_selector.NeighborhoodResidueSelector(chainB, 10.0, False)
interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 10.0, False)
big_interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 14.0, False)
interface_by_vector = core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_A, interface_on_B)
interface_by_vector.cb_dist_cut(11)
interface_by_vector.cb_dist_cut(5.5)
interface_by_vector.vector_angle_cut(75)
interface_by_vector.vector_dist_cut(9)



A_or_big_B = core.select.residue_selector.OrResidueSelector(chainA, big_interface_on_B)


xml = f"""

<RESIDUE_SELECTORS>
    <Chain name="chainA" chains="A"/>
    <Chain name="chainB" chains="B"/>
    <Slice name="patchdock_res" indices="{args.patchdock_res}" selector="chainB" />
</RESIDUE_SELECTORS>


<FILTERS>
    <ContactMolecularSurface name="contact_patch_res" apolar_target="true"
        distance_weight="0.5" target_selector="patchdock_res" binder_selector="chainA" confidence="0" />
    <ContactMolecularSurface name="contact_patch_res_quick" apolar_target="true"
        quick="1" distance_weight="0.5" target_selector="patchdock_res" binder_selector="chainA" confidence="0" />

</FILTERS>

"""


objs = protocols.rosetta_scripts.XmlObjects.create_from_string(xml)

patchdock_sel = objs.get_residue_selector("patchdock_res")

contact_patch_filt = objs.get_filter("contact_patch_res")
contact_patch_quick_filt = objs.get_filter("contact_patch_res_quick")

patchdock_res_1 = np.array([int(x) for x in args.patchdock_res.split(',')])

# contact_patch1_filt_sq5 = objs.get_filter("contact_patch_res1_sq5")
# contact_patch2_filt_sq5 = objs.get_filter("contact_patch_res2_sq5")
# contact_patch1_quick_filt_sq5 = objs.get_filter("contact_patch_res1_quick_sq5")
# contact_patch2_quick_filt_sq5 = objs.get_filter("contact_patch_res2_quick_sq5")

scorefxn_bb_hbond = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_bb_hbond.set_weight(core.scoring.hbond_sr_bb, 1)
scorefxn_bb_hbond.set_weight(core.scoring.hbond_lr_bb, 1)

scorefxn_bb_hbond2 = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_bb_hbond2.set_weight(core.scoring.hbond_sr_bb, 1)
scorefxn_bb_hbond2.set_weight(core.scoring.hbond_lr_bb, 1)
scorefxn_bb_hbond2.set_weight(core.scoring.hbond_bb_sc, 1)

    
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

fix_scorefxn(scorefxn_bb_hbond)
fix_scorefxn(scorefxn_bb_hbond2)
# fix_scorefxn(scorefxn)


def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string


def add_to_score_map(og_map, to_add, prefix, suffix=""):
    for name, score in list(to_add.items()):    # this iterator is broken. use list()
        og_map[prefix + name + suffix] = score

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose


def which_chain(pose, resi):
    for i in range(1, pose.num_chains()+1):
        if ( pose.conformation().chain_begin(i) <= resi and pose.conformation().chain_end(i) >= resi):
            return i
    assert(False)

def get_monomer_score(pose, scorefxn):
    pose = pose.split_by_chain()[1]
    return scorefxn(pose)


def get_filter_by_name(filtername):
    the_filter = objs.get_filter(filtername)

    # Get rid of stochastic filter
    if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
        the_filter = the_filter.subfilter()

    return the_filter

def add_filter_to_results(pose, filtername, out_score_map):
    filter = get_filter_by_name(filtername)
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    if (isinstance(filter, protocols.simple_filters.ShapeComplementarityFilter)):
        value = filter.compute(pose)
        out_score_map[filtername] = value.sc
        out_score_map[filtername+"_median_dist"] = value.distance
    else:
        value = filter.report_sm(pose)
        out_score_map[filtername] = value
    print("============End report for " + filtername + "=======================")

def score_with_these_filters(pose, filters, out_score_map):
    for filtername in filters:
        add_filter_to_results(pose, filtername, out_score_map)

def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

def atid(resnum, atno):
    return core.id.AtomID( atno, resnum )

abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))

def dump_region(pose, start, end, name):
    residues = utility.vector1_unsigned_long()
    for i in range(start, end ):
        residues.append(i)

    to_dump = core.pose.Pose()
    core.pose.pdbslice(to_dump, pose, residues)
    pdbinfo = core.pose.PDBInfo( to_dump )
    to_dump.pdb_info(pdbinfo)
    to_dump.dump_pdb(name)


def get_per_atom_sasa(pose):
    atoms = core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
    surf_vol = core.scoring.packing.get_surf_vol( pose, atoms, 2.8)
    # print(surf_vol.tot_surf)
    # print(surf_vol.surf(2, 1))  # this is per atom sasa (residue 2, atom 1)
    return surf_vol

# this is 1 indexed with the start and end as xx
# and HHHHHH turns identified
def better_dssp(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()
    the_dssp = list(the_dssp)
    the_dssp[1] = "x"
    the_dssp[-1] = "x"
    the_dssp[2] = "x"
    the_dssp[-2] = "x"
    the_dssp[3] = "x"
    the_dssp[-3] = "x"
    the_dssp = "".join(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp


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


# this is 1 indexed with the start and end with loops converted to nearby dssp
# and HHHHHH turns identified
def better_dssp3(pose, length=-1, force_consensus=None, consensus_size=6):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()[:length]
    the_dssp = list(the_dssp)

    n_consensus = get_consensus(the_dssp[3:consensus_size+1])
    if ( not force_consensus is None ):
        n_consensus = force_consensus

    for i in range(1, consensus_size+1):
        the_dssp[i] = n_consensus

    c_consensus = get_consensus(the_dssp[-(consensus_size):-2])
    if ( not force_consensus is None ):
        c_consensus = force_consensus

    for i in range(1, consensus_size+1):
        the_dssp[-i] = c_consensus

    the_dssp = "".join(the_dssp)

    # print(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( (abego == "B" or abego == "E") and seqpos > consensus_size and seqpos < len(the_dssp)-consensus_size ):
                this_dssp = "L"

        my_dssp += this_dssp

    # print(my_dssp)

    return my_dssp



# assumes dssp starts with X and removes it
def get_ss_elements2(dssp):
    assert(dssp[0] == "x")
    ss_elements = []

    offset = 0
    ilabel = -1
    for label, group in itertools.groupby(dssp):
        ilabel += 1
        this_len = sum(1 for _ in group)
        next_offset = offset + this_len

        ss_elements.append( (label, offset, next_offset-1))

        offset = next_offset
    return ss_elements[1:]



def superposition_transform_with_weight(input_coords_move, input_coords_to, rotation_matrix, move_com, ref_com, atom_weight):

    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()

    assert(len(input_coords_move) == len(input_coords_to))
    assert(len(input_coords_move) == len(atom_weight))

    for i in range(len(input_coords_move)):
        for j in range(atom_weight[i]):
            coords_move.append(input_coords_move[i+1])
            coords_to.append(input_coords_to[i+1])

    protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )


# align two things using sequence and accepting gaps
def pymol_align( move_pose, to_pose, sel_move=None, sel_to=None, atoms=["N", "CA", "C"], throw_away=0.1, extend_penalty=-0.1,
        to_pose_upweight_mask=None, move_pose_upweight_mask=None, return_all_move_to=False ):

    if ( to_pose_upweight_mask is None ):
        to_pose_upweight_mask = np.zeros(to_pose.size()+1, np.bool)

    if ( move_pose_upweight_mask is None ):
        move_pose_upweight_mask = np.zeros(move_pose.size()+1, np.bool)

    if ( not sel_move is None ):
        move_res = np.array(list(core.select.get_residues_from_subset(sel_move)))
    else:
        move_res = np.array(list(range(1, move_pose.size()+1)))

    if ( not sel_to is None ):
        to_res = np.array(list(core.select.get_residues_from_subset(sel_to)))
    else:
        to_res = np.array(list(range(1, to_pose.size()+1)))

    seq_move = "x" + move_pose.sequence()
    seq_to = "x" + to_pose.sequence()

    seq_move = "".join(np.array(list(seq_move))[move_res])
    seq_to = "".join(np.array(list(seq_to))[to_res])

    from Bio import pairwise2
    alignment = align_move, align_to, idk1, idk2, idk3 = pairwise2.align.globalxs(seq_move,seq_to, -2, extend_penalty, penalize_end_gaps=False)[0]
    # print(align_move, align_to)
    # print(alignment)

    all_move_to = []
    move_to_pairs = []
    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()
    atom_weight = []

    i_move = 0
    i_to = 0
    for i in range(len(align_move)):
        if ( align_move[i] == align_to[i] ):

            seqpos_move = move_res[i_move]
            seqpos_to = to_res[i_to]

            move_to_pairs.append((seqpos_move, seqpos_to))
            all_move_to.append((seqpos_move, seqpos_to))


            weight = 1
            if ( to_pose_upweight_mask[seqpos_to] ):
                weight = 5
            if ( move_pose_upweight_mask[seqpos_move] ):
                weight = 5


            for atom in atoms:
                coords_move.append(move_pose.residue(seqpos_move).xyz(atom))
                coords_to.append(to_pose.residue(seqpos_to).xyz(atom))
                atom_weight.append(weight)
        else:

            if ( align_move[i] != "-" and align_to[i] != "-" ):

                seqpos_move = move_res[i_move]
                seqpos_to = to_res[i_to]

                all_move_to.append((seqpos_move, seqpos_to))



        if ( align_move[i] != "-" ):
            i_move += 1
        if ( align_to[i] != "-" ):
            i_to += 1

    move_pose_copy = move_pose.clone()


    rmsd = 0

    distances = []

    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        superposition_transform_with_weight(coords_move, coords_to, rotation_matrix, move_com, ref_com, atom_weight)
        # protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )

        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance
                distances.append(distance)

        rmsd /= len(move_to_pairs)*len(atoms)
        rmsd = np.sqrt(rmsd)

    move_pose = move_pose_copy

    distances = np.array(distances)

    # print("Initial RMSD: %.3f"%rmsd)

    cutoff = np.percentile(distances, 100 - throw_away * 10)
    # print("Cutoff %.3f"%cutoff)

    mask = distances <= cutoff

    # print(mask.sum(), len(mask))

    coords_move_old = list(coords_move)
    coords_to_old = list(coords_to)
    atom_weight_old = list(atom_weight)
    # move_to_pairs_old = move_to_pairs

    # move_to_pairs = []
    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()
    atom_weight = []

    for i in range(len(coords_move_old)):
        if ( not mask[i] ):
            continue
        coords_move.append(coords_move_old[i])
        coords_to.append(coords_to_old[i])
        atom_weight.append(atom_weight_old[i])
        # move_to_pairs.append(move_to_pairs_old[i])

    # print(len(coords_move), len(coords_move_old))

    rmsd = 0
    imask = -1
    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        superposition_transform_with_weight( coords_move, coords_to, rotation_matrix, move_com, ref_com, atom_weight )

        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                imask += 1
                if ( not mask[imask] ):
                    continue
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance

        rmsd /= imask
        rmsd = np.sqrt(rmsd)

        zero = numeric.xyzVector_double_t(0, 0, 0)
        xform = nup.vector_to_xform( zero - ref_com ) @ nup.matrix_to_xform( rotation_matrix ) \
                    @ nup.vector_to_xform( move_com )


    print("Final RMSD: %.3f over %i atoms"%(rmsd, mask.sum()))


    if ( return_all_move_to ):
        return rmsd, move_to_pairs, move_pose, xform, all_move_to
    else:
        return rmsd, move_to_pairs, move_pose, xform


def pose_from_silent_lines(structure, tag):
    vec = utility.vector1_std_string()
    vec.append(tag)

    stream = std.istringstream(structure)

    sfd = core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd.read_stream(stream, vec, True, "fake")

    pose = core.pose.Pose()
    sfd.get_structure(tag).fill_pose(pose)

    return pose

import pyrosetta
def ros():
    return pyrosetta.rosetta
def pyro():
    return pyrosetta

class RosettaPacker:

    def __init__(self):

        self.chainA = ros().core.select.residue_selector.ChainSelector("A")
        self.chainB = ros().core.select.residue_selector.ChainSelector("B")
        self.interface_on_A = ros().core.select.residue_selector.NeighborhoodResidueSelector(self.chainB, 10.0, False)
        self.interface_on_B = ros().core.select.residue_selector.NeighborhoodResidueSelector(self.chainA, 10.0, False)
        self.AB_interface = ros().core.select.residue_selector.OrResidueSelector( self.interface_on_A, self.interface_on_B )
        self.Not_interface = ros().core.select.residue_selector.NotResidueSelector( self.AB_interface )
        self.chainA_not_interface = ros().core.select.residue_selector.AndResidueSelector( self.Not_interface, self.chainA )
        self.chainB_not_interface = ros().core.select.residue_selector.AndResidueSelector( self.Not_interface, self.chainB )

        self.scorefxn_insta = pyro().get_fa_scorefxn()
        for term in self.scorefxn_insta.get_nonzero_weighted_scoretypes():
            name = ros().core.scoring.name_from_score_type(term)

            if ( "_dun" in name ):
                continue
            if ( "rama" in name ):
                continue
            if ( "p_aa_pp" in name ):
                continue
            if ( "fa_rep" in name ):
                continue
            if ( "fa_atr" in name ):
                continue
            if ( "fa_sol" in name ):
                continue
            if ( "hbond" in name ):
                continue
            if ( "pro_close" in name ):
                continue

            self.scorefxn_insta.set_weight(term, 0)

        self.scorefxn_insta_soft = self.scorefxn_insta.clone()
        self.scorefxn_insta_soft.set_weight(ros().core.scoring.fa_rep, 0.15)


        self.scorefxn_none = ros().core.scoring.ScoreFunctionFactory().create_score_function("none")
        self.scorefxn_atr = ros().core.scoring.ScoreFunctionFactory().create_score_function("none")
        self.scorefxn_atr.set_weight(ros().core.scoring.fa_atr, 1)
        self.scorefxn_beta = pyro().get_fa_scorefxn()


    # pack with only dunbrack, vdw, and hbonds
    #  elec is super slow so we can't have that
    def insta_pack(self, pose):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.RestrictToRepacking() )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_insta )
        packer.task_factory( tf )

        packer.apply( pose )

    def beta_pack(self, pose):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.RestrictToRepacking() )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_beta )
        packer.task_factory( tf )

        packer.apply( pose )


    # use the packer to do this to prevent issues with disulfides and to ensure we get variant types right
    def thread_seq(self, pose, new_seq):

        old_seq = pose.sequence()

        locked_subset = ros().utility.vector1_bool( pose.size() )

        tf = ros().core.pack.task.TaskFactory()
        for seqpos in range(1, pose.size()+1 ):
            old_letter = old_seq[seqpos-1]
            new_letter = new_seq[seqpos-1]

            if ( old_letter == new_letter ):
                locked_subset[seqpos] = True
                continue

            restrict_aa = ros().core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            restrict_aa.aas_to_keep( new_letter )

            subset = ros().utility.vector1_bool( pose.size() )
            subset[seqpos] = True
            tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )

        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          locked_subset ) )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_none )
        packer.task_factory( tf )

        packer.apply( pose )


    def fast_hydrophobic_interface(self, pose):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.RestrictToRepackingRLT(),
                          self.chainB, False ))
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainA_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        restrict_aa = ros().core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        restrict_aa.aas_to_keep( "GAFILMVW" )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( restrict_aa,
                          self.interface_on_A, False ))

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_insta_soft )
        packer.task_factory( tf )

        packer.apply( pose )



rosetta_packer = RosettaPacker()


def array1_to_subset(array1):
    subset = utility.vector1_bool(len(array1)-1)
    for i in range(1, len(subset)+1):
        subset[i] = array1[i]
    return subset

def subset_to_array1(subset):
    size = len(subset)
    try:
        subset[0]
        size -= 1
    except:
        pass
    array1 = np.zeros(size+1, dtype=bool)
    for i in range(1, size+1):
        array1[i] = subset[i]
    return array1


def rigid_body_subset(pose, subset, translation=np.array([0, 0, 0]), rotation=np.identity(3), com=None):

    if ( np.sum(list(subset)) == 0 ):
        return

    ids = utility.vector1_core_id_AtomID()
    new_atoms = utility.vector1_numeric_xyzVector_double_t()

    atoms = []
    for seqpos in range(1, pose.size()+1):
        if ( not subset[seqpos] ):
            continue
        res = pose.residue(seqpos)

        for iatom in range(1, res.natoms()+1):

            idd = core.id.AtomID(iatom, seqpos)
            ids.append(idd)
            atoms.append(nup.from_vector(res.xyz(iatom)))

    atoms = np.array(atoms)
    if ( com is None ):
        com = np.mean( atoms, axis=0 )


    atoms = ( rotation @ ( atoms - com ).T ).T + com + translation

    for i in range(len(atoms)):
        new_atoms.append(nup.to_vector(atoms[i]))

    pose.batch_set_xyz(ids, new_atoms)


def subset_cms(pose, subset, monomer_size, away_vec, quick=False, sq5=False ):
    pose = pose.clone()

    subset = subset_to_array1(subset)
    subset = ~subset
    subset[monomer_size+1:] = False

    rigid_body_subset(pose, subset, translation=away_vec*30)

    if ( sq5 ):
        if ( quick ):
            cms1 = contact_patch_quick_filt_sq5.report_sm(pose)
        else:
            cms1 = contact_patch_filt_sq5.report_sm(pose)
    else:
        if ( quick ):
            cms1 = contact_patch_quick_filt.report_sm(pose)
        else:
            cms1 = contact_patch_filt.report_sm(pose)

    return cms1


# returns regions of repeating elements
# (value, start, stop)
# stop is the last element of the region for slicing using stop+1
def repeating_regions(vector, only_keep=None):
    offset = 0
    regions = []
    for value, group in itertools.groupby(vector):
        this_len = len(list(group))
        next_offset = offset + this_len
        if ( only_keep is None or only_keep == value ):
            regions.append( (value, offset, next_offset-1))
        offset = next_offset

    return regions



def delete_residues_smart(pose, start, end):

    ft = pose.fold_tree()
    if ( start == 1 ):
        ft.reorder(pose.size())
    else:
        ft.reorder(1)
    pose.fold_tree(ft)

    pose.delete_residue_range_slow(start, end)
    pose.conformation().detect_disulfides()


def delete_mask(pose, mask):

    assert(pose.size()+1 == len(mask))
    delete_regions = repeating_regions(mask, only_keep=True)

    for _, start, end in reversed(delete_regions):
        delete_residues_smart(pose, start, end)



def random_rigid(translation_scale, rotation_scale, in_rotation_unit=None):
    random_units = np.random.normal(size=(2, 3))
    random_units /= np.linalg.norm( random_units, axis=-1 )[:,None]

    rotation_unit = random_units[0]
    translation_unit = random_units[1]

    if ( not in_rotation_unit is None ):
        rotation_unit = in_rotation_unit

        rotation_scale *= np.sign( np.random.random() - 0.5 )


    rotation_scale = np.random.random()*rotation_scale
    translation_scale = np.random.random()*translation_scale

    rotation_matrix = scipy.spatial.transform.Rotation.from_rotvec(rotation_unit * np.radians( rotation_scale ) ).as_matrix()


    return translation_unit * translation_scale, rotation_matrix



def contact_patch_per_res(pose, start, end, monomer_size, away_unit):


    size = end-start+1

    # forward direction
    cp_forwards = []
    cp_backwards = []
    for i_trim in range(size):
        mask = np.zeros(pose.size()+1, bool)
        mask[start+i_trim:end+1] = True
        cms1 = subset_cms(pose, mask, monomer_size, away_unit, quick=True)
        cp_forwards.append(cms1 )

        mask = np.zeros(pose.size()+1, bool)
        mask[start:end+1-i_trim] = True
        cms1 = subset_cms(pose, mask, monomer_size, away_unit, quick=True)
        cp_backwards.append(cms1)

    cp_forwards.append(0)
    cp_backwards.append(0)

    cp_forwards = np.array(cp_forwards)
    cp_backwards = np.array(cp_backwards)

    # i - i+1 -- should mostly be positive asside from glitchy behavior
    del_forwards = (cp_forwards[:-1] - cp_forwards[1:]).clip(0, None)
    del_backwards = (cp_backwards[:-1] - cp_backwards[1:]).clip(0, None)

    per_res = (del_forwards + del_backwards[::-1]) / 2

    return per_res


def mean_pca(pcas):

    all_by_dot = np.sum( pcas[:,None] * pcas[None,:] )

    imost_parallel = np.argmax( np.abs(all_by_dot.sum(axis=0) ) )

    most_parallel = pcas[imost_parallel]

    for i in range(len(pcas)):
        pcas[i] = orient_a_vector_to_b( pcas[i], most_parallel )

    assert( np.all( np.sum( most_parallel[None,:] * pcas, axis=-1 ) >= 0 ) )

    mean_pca = np.mean(pcas, axis=0)
    mean_pca /= np.linalg.norm(mean_pca)

    return mean_pca




def orient_a_vector_to_b(a, b):
    
    if ( ( a * b ).sum() < 0 ):
        a = a * -1

    return a


def get_pca_and_orient_forwards(cas):

    pca = sklearn.decomposition.PCA()
    pca.fit(cas)
    axis = pca.components_[0].copy()
    
    axis /= np.linalg.norm(axis)
    
    # now we orient the PCAs so that they all point forwards
    #  this handles sharp kinks better than dotting with the previous pca
    first_to_last = cas[-1] - cas[0]

    axis = orient_a_vector_to_b( axis, first_to_last )

    return axis


def align_two_segs(move_pose, to_pose, move_start, move_end, to_start, to_end):

    coords_move = utility.vector1_numeric_xyzVector_double_t()
    for seqpos in range(move_start, move_end+1):

        res = move_pose.residue(seqpos)
        coords_move.append(res.xyz("N"))
        coords_move.append(res.xyz("CA"))
        coords_move.append(res.xyz("C"))

    coords_to = utility.vector1_numeric_xyzVector_double_t()
    for seqpos in range(to_start, to_end+1):

        res = to_pose.residue(seqpos)
        coords_to.append(res.xyz("N"))
        coords_to.append(res.xyz("CA"))
        coords_to.append(res.xyz("C"))

    rotation_matrix = numeric.xyzMatrix_double_t()
    move_com = numeric.xyzVector_double_t()
    ref_com = numeric.xyzVector_double_t()

    protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )

    protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)


def get_pairing_bounds(pairing):
    # all they had to do was implement begin2() end2()...

    ostream = std.ostringstream()
    pairing.show_internals(ostream)
    pair_string = ostream.str()
    begin1 = None
    end1 = None
    begin2 = None
    end2 = None
    for line in pair_string.split('\n')[:2]:
        if line.startswith("pairing1"):
            sp = line.strip().split()
            begin1 = int(sp[1])
            end1 = int(sp[-1])
            if end1 < begin1:
                begin1,end1 = end1,begin1
        if line.startswith("pairing2"):
            sp = line.strip().split()
            begin2 = int(sp[1])
            end2 = int(sp[-1])
            if end2 < begin2:
                begin2,end2 = end2,begin2
    assert(begin1 is not None)
    assert(begin2 is not None)
    if begin2 < begin1:
        begin1,end1,begin2,end2 = begin2,end2,begin1,end1

    return begin1, end1, begin2, end2

def strand_compare(start1, end1, start2, end2):
    if start1 >= start2 and start1 <= end2:
        return True, min(start1, start2), max(end1, end2)
    if end1 >= start2 and end1 <= end2:
        return True, min(start1, start2), max(end1, end2)
    if start1 <= start2 and end1 >= end2:
        return True, min(start1, start2), max(end1, end2)
    return False, None, None

def insert_strand(all_strands, start1, end1):
    same = False
    for i in range(len(all_strands)):
        start, end = all_strands[i]
        same, new_start, new_end = strand_compare(start, end, start1, end1)
        if same:
            all_strands.pop(i)
            insert_strand(all_strands, new_start, new_end)
            break
    if not same:
        all_strands.append((start1, end1))

def index_strand(all_strands, start1, end1):
    for i in range(len(all_strands)):
        start, end = all_strands[i]
        same, new_start, new_end = strand_compare(start, end, start1, end1)
        if same:
            return i

    assert(False)




def get_part_of_isheet(ss_elems, pose):
    pose = pose.split_by_chain()[1]

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()

    sps = dssp.strand_pairing_set()

    # first enumerate all strands

    all_strands = []
    for ipair in range(1, sps.size()+1):
        pairing = sps.strand_pairing(ipair)
        start1, end1, start2, end2 = get_pairing_bounds(pairing)

        insert_strand(all_strands, start1, end1)
        insert_strand(all_strands, start2, end2)

    # begin_to_istrand = {}
    begins = []
    ends = []
    all_strands2 = []
    for begin, end in sorted(all_strands, key=lambda x: x[0]):
        begins.append(begin)
        ends.append(end)
        all_strands2.append((begin, end))
    begins = np.array(begins)
    ends = np.array(ends)
    all_strands = all_strands2
    # for istrand, begin in enumerate(begins):
    #     begin_to_istrand[begin] = istrand

    # next figure out who's paired with who
    is_paired = np.zeros((len(all_strands), len(all_strands)), bool)

    # import IPython
    # IPython.embed()

    for ipair in range(1, sps.size()+1):
        pairing = sps.strand_pairing(ipair)
        start1, end1, start2, end2 = get_pairing_bounds(pairing)

        istrand1 = index_strand(all_strands, start1, end1)
        istrand2 = index_strand(all_strands, start2, end2)

        is_paired[istrand1,istrand2] = True
        is_paired[istrand2,istrand1] = True

    # define sheets
    import scipy.sparse.csgraph
    n_sheets, istrand_in_isheet = scipy.sparse.csgraph.connected_components(is_paired)

    # make a res mask for which sheet
    res_in_isheet = np.zeros(pose.size()+1)
    res_in_isheet[:] = np.nan
    for istrand in range(len(begins)):
        begin = begins[istrand]
        end = ends[istrand]
        isheet = istrand_in_isheet[istrand]
        res_in_isheet[begin:end+1] = isheet

    isheet_sizes = np.array([(res_in_isheet == isheet).sum() for isheet in range(n_sheets)])

    part_of_isheet = []
    for h, begin, end in ss_elems:
        if h != "E":
            part_of_isheet.append(np.nan)
            continue
        options, counts = np.unique(res_in_isheet[begin:end+1], return_counts=True)
        mask = ~np.isnan(options)
        options = options[mask]
        counts = counts[mask]
        if len(options) > 1:
            options = [options[np.argmax(counts)]]
        # if len(options) == 0:
        #     import IPython
        #     IPython.embed()

        assert(len(options) == 1)
        part_of_isheet.append(options[0])

    # import IPython
    # IPython.embed()

    # print('isheet', part_of_isheet)

    return np.array(part_of_isheet, dtype=float), n_sheets



def get_cross_chain_bb_hbonds(pose, cutoff=-0.25):

    scorefxn_bb_hbond( pose )

    monomer_size = pose.conformation().chain_end(1)

    cross_bb_hbonds = np.zeros(pose.size(), bool)

    gr = pose.energies().energy_graph()
    for i in range(1, monomer_size+1):
        for j in range(monomer_size, pose.size()+1):
            if gr.find_edge(i, j):
                score = gr.find_edge(i, j).dot( scorefxn_bb_hbond.weights() )
                if score < cutoff:
                    cross_bb_hbonds[i-1] = True
                    cross_bb_hbonds[j-1] = True

    return cross_bb_hbonds


def bb_to_target_hbond_vectors(pose, cutoff=-0.1):

    scorefxn_bb_hbond2(pose)

    hbset = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

    monomer_size = pose.conformation().chain_end(1)

    bb_hbond_vectors = [[] for x in range(monomer_size)]

    for hbond in hbset.hbonds():
        acc_res = hbond.acc_res()
        don_res = hbond.don_res()

        acc_is_binder = acc_res <= monomer_size
        don_is_binder = don_res <= monomer_size

        if acc_is_binder == don_is_binder:
            continue

        seqpos = None
        if acc_is_binder:
            seqpos = acc_res
            if not hbond.acc_atm_is_backbone():
                continue
            binder_atom = pose.residue( acc_res ).xyz( hbond.acc_atm() )
            target_atom = pose.residue( don_res ).xyz( hbond.don_hatm() )
        else:
            seqpos = don_res
            if not hbond.don_hatm_is_backbone():
                continue
            binder_atom = pose.residue( don_res ).xyz( hbond.don_hatm() )
            target_atom = pose.residue( acc_res ).xyz( hbond.acc_atm() )

        binder_xyz = nup.from_vector(binder_atom)
        target_xyz = nup.from_vector(target_atom)

        unit = target_xyz - binder_xyz
        unit /= np.linalg.norm(unit)

        bb_hbond_vectors[seqpos-1].append(unit)

    return bb_hbond_vectors



the_locals = None

def worst_possible_asp(pose, name_no_suffix, out_score_map, out_string_map, suffix):

    if align_pose is not None:
        rmsd, move_to_pairs, pose, xform = pymol_align(pose, align_pose, sel_move=chainB.apply(pose), sel_to=chainB.apply(align_pose))

    away_dist = 30

    npose = nup.npose_from_pose(pose)
    cas = nu.extract_atoms(npose, [nu.CA])[:,:3]

    monomer_size = pose.conformation().chain_end(1)
    target_size = pose.size() - monomer_size

    sequence1 = "x" + pose.sequence()

    patchdock_res_target_1 = patchdock_res_1 + monomer_size

    cas_patchdock = cas[patchdock_res_target_1-1]


    binder_npose = npose[:nu.R*monomer_size]
    binder_cas = nu.extract_atoms(binder_npose, [nu.CA])
    target_npose = npose[nu.R*monomer_size:]

    binder_com = np.mean(binder_cas, axis=0)[:3]

    dssp = better_dssp3(pose.split_by_chain()[1])
    ss_elems = get_ss_elements2(dssp)
    early_part_of_isheet, n_sheets = get_part_of_isheet(ss_elems, pose)

    # com of each strand so that the sheet axis can be found later
    sheet_coms = []
    sheet_3d_unit = []
    early_i_am_sheet_start = np.zeros(len(ss_elems), bool)
    early_i_am_sheet_end = np.zeros(len(ss_elems), bool)
    for isheet in range(n_sheets):
        mask = early_part_of_isheet == isheet
        these_coms = []
        for ielem in range(len(ss_elems)):
            if not mask[ielem]:
                continue
            h, start, end = ss_elems[ielem]
            these_coms.append( np.mean(binder_cas[start-1:end,:3], axis=0))
        sheet_coms.append(these_coms)

        # import IPython
        # IPython.embed()

        if len(these_coms) == 0:
            axis = np.array([np.nan, np.nan, np.nan])
            sheet_3d_unit.append(axis)
            continue

        pca = sklearn.decomposition.PCA()
        pca.fit(these_coms)
        axis = pca.components_[0].copy()
        axis /= np.linalg.norm(axis)
        sheet_3d_unit.append(axis)

        com = np.mean(these_coms, axis=0)
        projection = np.sum( (np.array(these_coms) - com) * axis, axis=-1)
        sheet_order = np.argsort(projection)

        start_ielem = np.where(mask)[0][sheet_order[0]]
        end_ielem = np.where(mask)[0][sheet_order[-1]]

        early_i_am_sheet_start[start_ielem] = True
        early_i_am_sheet_end[end_ielem] = True




    close_res = np.array([False] + list(interface_on_A.apply(pose)))


    cross_bb_hbonds = get_cross_chain_bb_hbonds(pose)

    bb_hbond_vectors = bb_to_target_hbond_vectors(pose)


    helices = []
    part_of_isheet = []
    i_am_sheet_start = []
    i_am_sheet_end = []
    for ielem, (h, start, end) in enumerate(ss_elems):
        if h == 'L':
            continue
        if h == 'E' and not args.allow_strand:
            continue
        helices.append([start,end])
        part_of_isheet.append(early_part_of_isheet[ielem])
        i_am_sheet_start.append(early_i_am_sheet_start[ielem])
        i_am_sheet_end.append(early_i_am_sheet_end[ielem])
    helices = np.array(helices)
    part_of_isheet = np.array(part_of_isheet)
    i_am_sheet_start = np.array(i_am_sheet_start)
    i_am_sheet_end = np.array(i_am_sheet_end)

    elem_pcas = []
    for start, end in helices:
        elem_pcas.append(get_pca_and_orient_forwards(cas[start-1:end]))
    elem_pcas = np.array(elem_pcas)

    helix_has_res = []
    for start, end in helices:
        helix_has_res.append( close_res[start:end+1].sum() )

    helix_has_res = np.array(helix_has_res)

    helix_close_enough = helix_has_res >= 2
    if ( args.dump_all_helices ):
        helix_close_enough[:] = True
    helices = helices[helix_close_enough]
    part_of_isheet = part_of_isheet[helix_close_enough]
    elem_pcas = elem_pcas[helix_close_enough]
    i_am_sheet_start = i_am_sheet_start[helix_close_enough]
    i_am_sheet_end = i_am_sheet_end[helix_close_enough]


    extra_mask = np.zeros(pose.size()+1, bool)
    extra_mask[1:monomer_size+1] = True
    extra_mask[monomer_size+1:] = False

    for start, end in helices:
        extra_mask[start:end+1] = False


    binder_com = np.mean(binder_npose, axis=0)[:3]
    target_com = np.mean(target_npose, axis=0)[:3]
    if ( args.com_of_patch_res ):
        target_com = np.mean(cas_patchdock, axis=0)[:3]
        

    away_unit = binder_com - target_com
    away_unit /= np.linalg.norm(away_unit)

    all_removed = np.zeros(pose.size()+1, bool)

    rigid_body_subset(pose, extra_mask, translation=away_unit * away_dist)
    all_removed |= extra_mask


    cms_remove_mask = np.zeros(pose.size()+1, bool)

    new_helices = []
    new_strand_in_isheet = []
    new_i_am_sheet_start = []
    new_i_am_sheet_end = []
    new_pcas = []
    for ihelix, (start, end) in enumerate(helices):

        our_isheet = part_of_isheet[ihelix]
        we_are_strand = not np.isnan(our_isheet)

        cp_per_res = contact_patch_per_res(pose, start, end, monomer_size, away_unit)
        print(start, end)
        print(cp_per_res)
        print(cp_per_res.sum())

        for i in range(len(cp_per_res)):
            # if ( cp_per_res[i] > args.cp_per_res_cut ):
            #     continue
            seqpos = i + start

            if ( args.poly_letter != "" ):
                protocols.toolbox.pose_manipulation.repack_this_residue(seqpos, pose, scorefxn, False, args.poly_letter)



        for i in range(len(cp_per_res)):

            seqpos = i + start

            if cross_bb_hbonds[seqpos-1]:
                cp_per_res[i] += 30


        new_start = None
        new_end = None

        # trim until we hit a residue with enough contact patch
        for i in range(len(cp_per_res)):
            real_i = i
            seqpos = real_i + start

            ok = cp_per_res[real_i] > args.cp_per_res_cut

            if ( ok ):
                new_start = seqpos
                break

            # cms_remove_mask[seqpos] = True

        # trim reverse until we hit a residue with enough contact patch
        for i in range(len(cp_per_res)):
            real_i = len(cp_per_res) - i - 1
            seqpos = real_i + start

            ok = cp_per_res[real_i] > args.cp_per_res_cut

            if ( ok ):
                new_end = seqpos
                break

        # if there's only 1 residue and that residue is R or K, discard this motif
        if ( (cp_per_res > args.cp_per_res_cut).sum() == 1 ):
            who_is = np.where(cp_per_res > args.cp_per_res_cut)[0][0] + start
            if ( sequence1[who_is] in "RK" ):
                new_start = None
                new_end = None

        if ( cp_per_res.sum() < args.min_helix_cp ):
            new_start = None
            new_end = None
            # cms_remove_mask[seqpos] = True
        if ( args.dump_all_helices ):
            new_start = start
            new_end = end

        if ( new_start is None ):
            assert(new_end is None)
        else:

            size = new_end - new_start + 1
            min_size = args.min_strand_size if we_are_strand else args.min_helix_size
            if ( size < min_size ):
                missing = min_size - size
                add_start = int(np.ceil(missing/2))
                add_end = int(np.floor(missing/2))
                new_start -= add_start
                new_end += add_end
                if ( new_start < start ):
                    shift = start - new_start
                    new_start += shift
                    new_end += shift
                if ( new_end > end ):
                    shift = new_end - end
                    new_end -= shift
                    new_start -= shift
                if ( new_start < start ):
                    assert(start-end+1 < args.min_helix_size)
                    new_start = None
                    new_end = None

            if ( new_start is not None ):
                new_helices.append((new_start, new_end))
                new_strand_in_isheet.append(our_isheet)
                new_pcas.append(elem_pcas[ihelix])
                new_i_am_sheet_start.append( i_am_sheet_start[ihelix] )
                new_i_am_sheet_end.append( i_am_sheet_end[ihelix] )

                for seqpos in range(start, end+1):
                    if ( seqpos < new_start or seqpos > new_end ):
                        cms_remove_mask[seqpos] = True

        if ( new_start is None ):
            cms_remove_mask[start:end+1] = True


    new_strand_in_isheet = np.array(new_strand_in_isheet)
    new_i_am_sheet_start = np.array(new_i_am_sheet_start)
    new_i_am_sheet_end = np.array(new_i_am_sheet_end)


    if ( len(new_helices) == 0 ):
        return None

    new_mean_pca = mean_pca(np.array(new_pcas))

    new_coms = []
    for start, end in new_helices:
        com = np.mean(cas[start-1:end], axis=0)
        new_coms.append(com)

    origin = np.mean(new_coms, axis=0)

    motif_com = origin.copy()

    # y-axis is the rejection of away_vec from z_axis

    z_axis = new_mean_pca

    projection = np.sum( z_axis * away_unit ) * z_axis
    rejection = away_unit - projection

    y_axis = rejection / np.linalg.norm(rejection)
    x_axis = np.cross( y_axis, z_axis)

    def to_xy_coords(pt, ray=False):

        if ( not ray ):
            pt = pt - origin

        x = np.sum( pt * x_axis )
        y = np.sum( pt * y_axis )

        return np.array([x, y])

    def to_z_coord(pt, ray=False):

        if ( not ray ):
            pt = pt - origin

        z = np.sum( pt * z_axis )

        return z


    assert(np.linalg.norm(to_xy_coords(z_axis, True)) < 0.01)
    assert(np.isclose(to_xy_coords(x_axis, True)[0], 1))
    assert(np.isclose(to_xy_coords(y_axis, True)[1], 1))


    xy_coms = []
    for com in new_coms:
        xy_coms.append(to_xy_coords(com))

    helix_aligned_with_z = []
    for pca in new_pcas:
        aligned = np.sum( pca * z_axis ) > 0
        helix_aligned_with_z.append(aligned)

    helix_sizes = []
    helix_z_bounds = []
    for start, end in new_helices:
        helix_sizes.append(end-start+1)
        z_bounds = np.array( [to_z_coord(cas[start-1]), to_z_coord(cas[end])] )
        helix_z_bounds.append(z_bounds)

    needed_isheets = np.unique(new_strand_in_isheet)
    needed_isheets = list(needed_isheets[~np.isnan(needed_isheets)])
    final_strand_in_isheet = []
    for old_isheet in new_strand_in_isheet:
        if np.isnan(old_isheet):
            final_strand_in_isheet.append(np.nan)
            continue
        final_strand_in_isheet.append(needed_isheets.index(old_isheet))

    isheet_units = []
    sheet_can_extend = []
    for isheet in needed_isheets:
        isheet = int(np.rint(isheet))
        coms = sheet_coms[isheet]
        twod_coms = [to_xy_coords(com) for com in coms]

        pca = sklearn.decomposition.PCA()
        pca.fit(twod_coms)
        axis = pca.components_[0].copy()
        
        axis /= np.linalg.norm(axis)

        threed_unit = to_xy_coords(sheet_3d_unit[isheet], ray=True)
        if np.sum( axis * threed_unit ) < 0:
            axis *= -1

        isheet_units.append(axis)


        strand_mask = new_strand_in_isheet == isheet
        is_start = new_i_am_sheet_start[strand_mask]
        is_end = new_i_am_sheet_end[strand_mask]

        assert( is_start.sum() <= 1)
        assert( is_end.sum() <= 1)

        extend_start = False
        extend_end = False

        if not is_start.any():
            extend_start = True # the start sheet got deleted
        else:
            who_is_start = np.where(strand_mask)[0][is_start][0]

            strand_start, strand_end = new_helices[who_is_start]
            any_bb_to_target = False
            for seqpos in range(strand_start, strand_end+1):
                if len(bb_hbond_vectors[seqpos-1]) > 0:
                    any_bb_to_target = True

            if not any_bb_to_target:
                extend_start = True

        if not is_end.any():
            extend_end = True # the end sheet got deleted
        else:
            who_is_end = np.where(strand_mask)[0][is_end][0]

            strand_start, strand_end = new_helices[who_is_end]
            any_bb_to_target = False
            for seqpos in range(strand_start, strand_end+1):
                if len(bb_hbond_vectors[seqpos-1]) > 0:
                    any_bb_to_target = True

            if not any_bb_to_target:
                extend_end = True


        # import IPython
        # IPython.embed()

        sheet_can_extend.append((extend_start, extend_end))



    xy_str = ",".join("[%.3f,%.3f]"%(x, y) for x, y in xy_coms)
    aligned_str = ",".join(("True" if aligned else "False") for aligned in helix_aligned_with_z)
    z_bounds = ",".join("[%.3f,%.3f]"%(n, c) for n, c in helix_z_bounds)
    sizes = ",".join("%i"%size for size in helix_sizes)
    strand_isheet_str = ",".join(['np.nan' if np.isnan(x) else str(int(x)) for x in final_strand_in_isheet ])
    isheet_units_str = ",".join("[%.3f,%.3f]"%(x, y) for x, y in isheet_units)
    sheet_can_extend_str = ",".join("[%s,%s]"%(s, e) for s, e in sheet_can_extend)


    # import IPython
    # IPython.embed()


    rigid_body_subset(pose, cms_remove_mask, translation=away_unit * away_dist)
    all_removed |= cms_remove_mask

    removed_com = np.mean(cas[all_removed[1:]], axis=0)[:3]


    delete_mask(pose, all_removed)

    og_contig_parts = []

    og_sizes = []
    prev_end = 0
    final_helices = []
    final_monomer_size = 0
    for start, end in new_helices:
        og_sizes.append(start - prev_end - 1)
        prev_end = end

        assert( not np.any(all_removed[start:end+1] ))
        adjust = all_removed[:start].sum()

        start -= adjust
        end -= adjust

        final_helices.append((start, end))

        final_monomer_size += end-start+1

        og_contig_parts.append("A%i-%i"%(start, end))

    og_sizes.append(monomer_size - prev_end)


    out_string_map['og_contig'] = ','.join(og_contig_parts)


    if ( args.target_poly_letter != "" ):
        for seqpos in range(final_monomer_size+1, pose.size()+1):
            protocols.toolbox.pose_manipulation.repack_this_residue(seqpos, pose, scorefxn, False, args.target_poly_letter)

    assert(pose.size()-final_monomer_size == target_size)


    # target_pose = pose.clone()
    # delete_residues_smart(target_pose, 1, final_monomer_size)
    # target_atoms, target_radii = nup.get_atoms_and_radii(target_pose)
    # target_atoms = np.array(target_atoms)
    # target_radii = np.array(target_radii)
    # heavy_mask = target_radii > 1.3
    # target_atoms = target_atoms[heavy_mask]
    # target_radii = target_radii[heavy_mask]

    # clashgrid = nu.clashgrid_from_points(target_atoms, target_radii, 0.5)


    if ( os.path.exists("test.pdb") ):
        os.remove("test.pdb")

    align_helix_size = 6
    align_strand_size = 3

    can_add_to_helix = []
    for ihelix, (start, end) in enumerate(final_helices):
        print(start, end)

        we_are_strand = not np.isnan(final_strand_in_isheet[ihelix])

        target_pose = pose.clone()
        remove_mask = np.zeros((pose.size()+1), bool)
        remove_mask[start:end+1] = True
        delete_mask(target_pose, remove_mask)

        target_atoms, target_radii = nup.get_atoms_and_radii(target_pose)
        target_atoms = np.array(target_atoms)
        target_radii = np.array(target_radii)
        heavy_mask = target_radii > 1.3
        target_atoms = target_atoms[heavy_mask]
        target_radii = target_radii[heavy_mask]

        clashgrid = nu.clashgrid_from_points(target_atoms, target_radii, 0.5)

        ideal_pose = ideal_strand.clone() if we_are_strand else ideal_helix.clone()
        align_size = align_strand_size if we_are_strand else align_helix_size

        # nterm
        test_ideal_pose = ideal_pose.clone()
        ideal_size = test_ideal_pose.size()

        align_two_segs(test_ideal_pose, pose, ideal_size-align_size+1, ideal_size, start, start + align_size-1)

        # a = target_pose.clone()
        # a.append_pose_by_jump(test_ideal_helix, 1)
        # core.io.pdb.add_to_multimodel_pdb(a, "test.pdb", "asdf")

        total_clashes = 0
        for icheck in range(ideal_size-align_size):
            seqpos = ideal_size - align_size - icheck
            res = test_ideal_pose.residue(seqpos)
            atoms = np.array([
                    nup.from_vector(res.xyz("N")),
                    nup.from_vector(res.xyz("CA")),
                    nup.from_vector(res.xyz("C")),
                    nup.from_vector(res.xyz("CB")),
                    ])
            n_clashes = clashgrid.clash_check(atoms, 100)
            n_clashes += total_clashes
            if ( n_clashes >= args.max_clashes_till_stop ):
                break
        this_can_add = [icheck]


        # cterm
        test_ideal_pose = ideal_pose.clone()

        align_two_segs(test_ideal_pose, pose, 1, align_size, end-align_size+1, end)

        # a = target_pose.clone()
        # a.append_pose_by_jump(test_ideal_helix, 1)
        # core.io.pdb.add_to_multimodel_pdb(a, "test.pdb", "asdf")

        total_clashes = 0
        for icheck in range(ideal_size-align_size):
            seqpos = align_size + icheck
            res = test_ideal_pose.residue(seqpos)
            atoms = np.array([
                    nup.from_vector(res.xyz("N")),
                    nup.from_vector(res.xyz("CA")),
                    nup.from_vector(res.xyz("C")),
                    nup.from_vector(res.xyz("CB")),
                    ])
            n_clashes = clashgrid.clash_check(atoms, 100)
            n_clashes += total_clashes
            if ( n_clashes >= args.max_clashes_till_stop ):
                break
        this_can_add.append(icheck)

        can_add_to_helix.append(this_can_add)

    can_add = ",".join("[%i,%i]"%(n,c) for (n,c) in can_add_to_helix)

    extra = ""
    if args.allow_strand:
        extra = ';strand_isheet=[%s];isheet_units=[%s];sheet_can_extend=[%s]'%(strand_isheet_str, isheet_units_str, sheet_can_extend_str)

    twod_string = 'xy_str=[%s];aligned_w_z=[%s];motif_sizes=[%s];motif_can_add=[%s]'%(xy_str, aligned_str, sizes, can_add)
    twod_string += extra
    print("Twod_string:", twod_string)



    out_score_map['n_helices'] = np.isnan(final_strand_in_isheet).sum()
    out_score_map['n_strands'] = (~np.isnan(final_strand_in_isheet)).sum()
    out_score_map['n_sheets'] = len(isheet_units)


    pdb_info = core.pose.PDBInfo(pose)

    for i in range(1, final_monomer_size+1):
        seqpos = i
        pdb_info.number(seqpos, i)
        pdb_info.chain(seqpos, "A")

    for i in range(1, target_size+1):
        seqpos = i + final_monomer_size 
        pdb_info.number(seqpos, i)
        pdb_info.chain(seqpos, "B")

    binder_contig_parts = ["A%i-%i"%(start,end) for start,end in final_helices]

    binder_contig = ",".join(binder_contig_parts) + ",0"
    target_contig = "B%i-%i"%(1, target_size)

    contig = binder_contig + " " + target_contig

    print("JHR_CONTIG:", contig)

    pdb_info.add_reslabel(1, "JHR_CONTIG:" + contig.replace(" ", "_"))
    pdb_info.add_reslabel(1, "TWOD_STRING:" + twod_string)

    pdb_info.add_reslabel(1, "BINDER_COM:" + '_'.join(['%.3f'%x for x in binder_com]))
    pdb_info.add_reslabel(1, "MOTIF_COM:" + '_'.join(['%.3f'%x for x in motif_com]))
    pdb_info.add_reslabel(1, "REMOVED_COM:" + '_'.join(['%.3f'%x for x in removed_com]))

    pose.pdb_info(pdb_info)


    og_contig_parts = []
    assert(len(og_sizes)-1 == len(binder_contig_parts))
    for i in range(len(binder_contig_parts) + 1):
        og_size = og_sizes[i]
        if ( og_size > 0 ):
            og_contig_parts.append("%i-%i"%(og_size, og_size))

        if ( i < len(binder_contig_parts) ):
            og_contig_parts.append(binder_contig_parts[i])

    og_contig = ",".join(og_contig_parts) + ",0 " + target_contig 
    print("OG Contig:", og_contig)






    return pose






############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())



# ckpt = 0
# if (os.path.exists("ckpt")):
#     with open("ckpt") as f:
#         try:
#             ckpt = int(f.read().strip())
#         except:
#             pass


num = -1
for ipdb, pdb in enumerate(pdbs):
    t0 = time.time()
    print("Attempting pose: " + pdb)


    # if ( ckpt > ipdb ):
    #     print("Checkpoint: Already finished: " + pdb)
    #     continue
    # with open("ckpt", "w") as f:
    #     f.write(str(ipdb))

    # try:
    for k in [1]:
        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")

        sfd = core.io.raw_data.ScoreFileData("score.sc")

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()


        out_pose = worst_possible_asp(pose, name_no_suffix, score_map, string_map, "")


        core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
        core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)

        sfd.write_pose( pose, score_map, name_no_suffix, string_map)
        if (out_pose != None):

            if ( silent == "" ):
                out_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                silent_name = "out.silent"
                sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(out_pose, name_no_suffix)
                for key in score_map:
                    struct.add_energy(key, score_map[key])
                for key in string_map:
                    struct.add_string_value(key, string_map[key])
                sfd_out.add_structure(struct)
                sfd_out.write_all(silent_name, False)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



# if ( silent != "" ):
#     sfd_out.write_all("out.silent", False)




















