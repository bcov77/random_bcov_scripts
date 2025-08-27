#!/usr/bin/env python

# A script that will estimate the binding affinity from experiments similar to the
# ones performed in Cao et. al 2021

# These values were originally refered to as Yeast-KD values, and were changed to SC50 in the
#  late stages of preparation. KD in scripts may be thought of as SC50.

# Brian Coventry 2021
#   -- Based in concept on the EC50 calculator script from Rocklin 2017, Science

# The KD used here is the KD of a yeast cell in a sorting experiment
#  (i.e. the concentration of target where 50% of that yeast cell will be collected)




# Actual operating principle of this script
#
# - Try to assess P( KD == X | sorting_data )
# - for each sort: P( KD == X | this_sort )
# - Overall distribution is found by multiplying each sort's P() distribution together
# - Find a KD upper-bound and lower-bound from the resulting distribution
#
# -- The process works by looking at the number of cells in the parent pool, the 
#      number of cells in the child pool, and the fraction of cells collected.
#      By combining all of this information, one can derive the % of this design
#      that was collected by the cell sorter. From there, you can use a binomial
#      distribution to ask questions about how likely is this if the KD == X.

# This script can also make point-estimates from single sorts


# Can you use this script?
#
# Here are critical experimental assumptions that are made. If these are followed, you can probably
#  analyze your data with this script. If they aren't, the results may not be valid:
#
# 1. You are performing a sorting analysis where you collect the binding population
# 2. Your gate is set such that nearly none of the cells are collected when your
#         binder target is not present (negative control). See Figure 1
# 3. Your gate is set to collect *everything* that shows signal above the negative control
#         See Figure 1
# 4. You are fully ngs sequencing every pool of interest
# 5. Your selection step is non-competitive (i.e. Strong binders won't prevent weak binders
#         from being collected) (Enrichment from one round to the next is fine, but 
#         sorting at 1pM in a small volume, taking specifically the top 1% of outputs, or running
#         your experiment in a column are no-gos)
# 6. You are recording the % of cells collected
# 7. Each of your cell-lines (designs) is expressing high enough such that 20% would get
#         collected if you just sorted on expression


# Figure 1. How to set your gate after you run your negative control
#             (Here, negative control means there is no fluorescently labeled target)

#    #                                                                                              
#    #                                                    `........-::::::::/ooooo-                 
#    #                                              `NNddddyyyyyyyyyooooooooo/://Ms                 
#    #                                              /M/                         -Mo                 
#    #                                              yM`                         /M-                 
#    #                                              Nd                          /M-                 
#    #                                             -Mo                          /M-                 
#    #                                             sM.    Set your gate like    yN.                 
#    #                                             mm           this            hN                  
#    #                                            .Ms                           hN                  
#  B #                                            +M:                           hN                  
#  i #                                            hN`                           My                  
#  n #                This gap is too big.       `Nh                            Ms                  
#  d #                 In reality, you           :M+                            Ms                  
#  i #                 should collect about      sM.                           .Mo                  
#  n #                 0.0002 of your            md                            /M:                  
#  g #                  negative control       .Mo              ..::/++syhddmdmm/                  
#    #                                  \       +M/-://ooyyhddmddhysoo//:-.                         
#  A #                                   \___   hNyss++/::..`             ``                        
#  x #                                           `      `..-::+++ssydddNNNMM:                       
#  i #               hdddddddddddddddddddddd:   +yhddNMMMMMMMMMMMMMMMMMMMMm-                        
#  s #               -mMMMMMMMMMMMMMMMMMMMN:    `yMMMMMMMMMMMMMMMMMMMMMMMs`                         
#    #                `yMMMMMMMMMMMMMMMMMN-       :NMMMMMMMMMMMMMMMMMMMm/                           
#  ( #                  /MMMMMMMMMMMMMMMm.         `sMMMMMMMMMMMMMMMMMh`                            
#  l #                   .mMMMMMMMMMMMMd`            -dMMMMMMMMMMMMMN+                              
#  o #                     sMMMMMMMMMMh`              `oMMMMMMMMMMMm-                               
#  g #                      :NMMMMMMMy                  .hMMMMMMMMs                                 
#  ) #                       .dMMMMMo                     /NMMMMm-                                  
#    #                         sMMM/                       .yMMy`                                   
#    #                          :y:                          -:                                     
#    #                                                                                              
#    #                                                                                              
#    #                                                                                              
#    #                                                                                              
#    #                                                                                              
#    #                                                                                              
#    #ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd/ 
#                                  Expression Axis (log)

# https://www.text-image.com/convert/pic2ascii.cgi



# Assumptions this script makes:
#
# 1. Each cell-line (design) expresses perfectly such that 100% of cells would be collected
#      if you only sorted by expression
#   -- Why: This is done to avoid trying to fit two different numbers to each design
#   -- Effect: Higher vs lower expressing designs will get about a 10-fold difference in KD
#                favoring the better expressing design (think 75% expression vs 25%)
#
# 2. Binding levels 20% and higher are all equivalent
#   -- Why: We know assumption 1 isn't correct. So we just consider anything over 20% to be saturating
#   -- Effect: KDs are fit more accuratly. We don't get really confused on binders that only express
#                 30% and therefore saturate at 30%.
#
# 3. Doubly-transformed yeast cells are the largest source of noise
#   -- Why: Based on assumption 4, it's hard to believe a source of noise could be larger
#   -- Effect: Probably pretty accurate, see 4.
#
# 4. We can avoid doubly-transformed yeast cell issues by setting a floor on the number
#      of counts we expect to see. Only counts higher than this are significant.
#      Specifically, we set this floor to most-copied-cell-since-transformation.
#      In practice, this number is:  (cells_pool / cells_expression).max() * cell_copies_before_expression_sort
#   -- Why: Doubly-transformed yeast will carry non-binding design DNA through the
#             experiment. If we assume that only 1 double-transformed yeast cell in
#             the transformation is the cause of each of these non-binding DNA sequences
#             that get pulled through the experiment (i.e. the same piece of DNA didn't
#             get lucky twice), then, the most number of cells we would expect this
#             non-binding DNA sequence to have would be equal to the number of copies
#             of that initial doubly-transformed yeast cell. Then, we just take the
#             worst case scenario for each pool and call that the floor
#   -- Effect: Doubly-transformed yeast issues are largely eliminated. However,
#             designs with low transformation efficiency may be considered non-binders
#
# 



# Weird things this script does:
#
# 5. Model the cell sorter as a binomial distribution
#   -- Why: This allows us to ask: P( KD == X | data ). Specifically,
#             we say, if there were 100 cells, and the KD says that 2% should
#             be collected, what are the odds of seeing 5 cells collected?
#
# 6. Linearly interpolate the x parameter of binomial distribution binom(x, p, n)
#   -- Why: The data isn't pefect integers. Sometimes you have 0.1 cells. This obviously
#             isn't great, but it follows the binomial better than estimating it with a
#             normal distribution does.
#
# 7. Calculate confidence intervals in a weird way
#   -- Why: Because 8. Let's pretend we're doing a 95% confidence interval. Instead of taking 
#             the cumsum(P(KD)) and looking at 0.025 and 0.975. This script just looks for the two
#             places where the P(KD) graph crosses 0.025
#
# 8. Boost the P(KD) graph by the best P(KD). (i.e. P(KD) /= P(KD).max()**2 )
#   -- Why: What do you do when P(KD).max() == 1e-20? Typically when this happens, the P(KD)
#             distribution has a very sharp peak, and if you just normalized it such that
#             P(KD).max() == 1, then, you're going to get a very narrow confidence interval.
#             But clearly, the confidence interval should be wide because you had P(KD).max() == 1e-20
#             and the assumptions in this script are falling apart. 
#             The adjustment here normalizes P(KD).max() == 1, and then does that same adjustment again
#             to boost the P(KD).max() == 1e20. Using 7. above to calculate the confidence interval,
#             we will now end up with a very wide interval, which fits with our confidence level. 







# Input data format:
# 
# This script takes 2 primary pieces information plus random parameters:
#
# ./estimate_affinity_from_ngs.py pooled_counts.list sorts.csv
#
# pooled_counts.list -- A file containing your experimental counts (whitespace separated)
#   -- Must contain a column for each pool to be analyzed
#   -- Must contain a "description" column with the name of each design
#
#   -- The absolute sum of each column does not matter. sorts.csv collected_cells takes care of this
#   -- The column titles must match the labels in sorts.csv
#
#   -- Example header: bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8 description
#
# sorts.csv -- A file containing the details of your sorting experiment (comma separated)
#   -- Sort.csv must contain these columns:
#      -- pool_name -- the name of this pool in pooled_counts.list
#      -- parent_pool -- the name of the pool that entered the cell sorter to produce this pool.
#                           If the parent pool doesn't exist in pooled_counts.list, leave this blank
#      -- concentration -- the concentration of the target protein. Historically in nM (but it doesn't matter)
#      -- avidity -- { avid, "" }. If you sorted with avidity, set this to avid. Blank otherwise
#      -- special -- { expression, other, "" }. A keyword to signal this wasn't a normal sort. This script only
#                       recognizes expression and "", anything else causes the column to be ignored. You must
#                       have an expression column however. (The sort where you only sort by expression)
#      -- special_num -- Not used by this script. Just leave it blank
#      -- collection_fraction -- The fraction of cells the cell sorter collected. Specifically, the % of
#                        cells collected by your final gate. (So, ignore anything that had bad SSC or FSC)
#                        On our cell sorter, this is called % parent.
#      -- collected_cells -- How many cells did you collect? It's ok to estimate a bit here.
#      -- non_standard_gate -- Leave this blank. If you did something weird, set any value here and this
#                                row will get rejected.
#      -- notes -- You can type anything you like here.
#
# Here is an example header for sorts.csv:
#   pool_name,parent_pool,concentration,avidity,special,special_num,collection_fraction,collected_cells,non_standard_gate,notes
#
# If you make a sorts.csv file with that as the only line and then open it with Excel. You'll be in a good starting place.
#
#
#
#
# If you want to see examples of properly formatted input files. See the sorting_ngs_data/ section of Cao 2021




# Output parameters:
#
# - kd_lb -- The lowerbound estimate on the confidence interval of the KD value (in the same units as sorts.csv)
# - kd_ub -- The upperbound estimate on the confidence interval of the KD value (in the same units as sorts.csv)
# - low_conf -- The data had too few counts to analyze. We can't say much about this design.
# - avid_doesnt_agree -- If your data has avid rounds, this means that the predicted KD from the avidity sorts
#                          is in total disagreement with the non-avid sorts. Be very wary of designs with 
#                          this flag set. Although some may be real, others may be non-specific garbage or
#                          doubly-transformed yeast cells.
# - avid_lb -- Same as kd_lb but for the avidity sorts (no correction for the concentration effect of avidity is made. That's up to you)
# - avid_ub -- Same as kd_ub but for the avidity sorts (no correction for the concentration effect of avidity is made. That's up to you)
# - 4000_nm_binder -- See --point_estimates. Whether or not this design looks like it's a KD < 3000_nm binder based on a single
#                        sort closest to 1000nm.
# - 400_nm_binder -- See --point_estimates. Whether or not this design looks like it's a KD < 300_nm binder based on a single
#                        sort closest to 100nm.
# - description -- The name of the design





import sys
import os
import argparse
import scipy.stats
import numpy as np
import pandas as pd



parser = argparse.ArgumentParser(description="Open this file with a text editor for very important information!")
parser.add_argument("pooled_counts.list", type=str)
parser.add_argument("sorts.csv", type=str)
parser.add_argument("--point_estimates", type=str, default="4000,400", help="Comma separated KDs to make point-estimates at."
                                                                        +" Will look for binders with 20% collection frequency at a concentration"
                                                                        +" 3x lower than the KD. If no such sort exists, will look for"
                                                                        +" a sort with a nearby concentration and adjust the collection frequency"
                                                                        +" according to the definition of KD (frac_bound = conc / (KD + conc) )")
parser.add_argument("--cell_copies_before_first_sort", type=float, default=1, help="Between transformation and your first sort with target,"
                                                                                +" how many copies of each cell happened? (i.e. this is"
                                                                                +" 2**cell_divisions)")
parser.add_argument("--min_selection_rate", default=0.0001, type=float, help="What sort of background collection rate do you expect from"
                                                                            +" your experiment?")
parser.add_argument("--collection_rate_saturation", default=0.20, type=float, help="Above what collection rate do we consider a design"
                                                                                +" to be saturated? See Assumptions 1. and 2.")
parser.add_argument("--max_parent_cells", default=1000, type=float, help="Max parent_pool cells used for single-design calculation. Counts are simply"
                                                                        +" scaled if they exceed this. We know some of the assumptions are wrong,"
                                                                        +" this prevents the binomial from giving incredibly sharp outputs for"
                                                                        +" designs with lots of counts.")
parser.add_argument("--mask_kd_with_avid_sorts", action="store_true", help="If the avid and non-avid sorts don't agree, by default, the "
                                                                        +" avid_doesnt_agree flag is set. Using this option goes further"
                                                                        +" and prevents kd_lb and kd_ub from extending to stronger binding"
                                                                        +" than what the avidity allowed. You'll get more robust outputs"
                                                                        +" with this at the cost of losing some real binders.")
parser.add_argument("--confidence_interval_width", default=-6, type=float, help="This script will produce a 1-exp(confidence_interval_width) "
                                                                                +"confidence interval. The default of -6 produces a "
                                                                                +"99.7% confidence interval.")
parser.add_argument("--kd_fold_step", default=1.1, type=float, help="How small of bins should we use for KD-estimates. "
                                                                                +"Default steps of 1.1-fold.")
parser.add_argument("--p_dna_is_sequenced", default=0.90, type=float, help="A rather insignificant parameter. This is only used"
                                                                         +" to prevent divide-by-0 for cells with 0 counts.")
parser.add_argument("--target", default="", type=str, help="Add a target column to your output.")

parser.add_argument("--dump_round_information", action="store_true", help="Output counts and fraction-passed-the-gate.")
parser.add_argument("--override_every_round_with_avid", action="store_true", help="Add this if you did every round with avidity")
parser.add_argument("--output_file", default="affinity_estimates.sc", type=str, help="The output file.")


parser.add_argument("--sketch_kd_threshold", default=1/20, type=float, help="sketch_kd_level is the fraction of the enrichment of this design to the"
                                                                            " enrichment of the most enriched binder with the same kd. This parmeter"
                                                                            " sets the threshold to be considered sketch_kd.")


args = parser.parse_args(sys.argv[1:])


point_estimates = []
for est_str in args.point_estimates.split(","):
    if ( len(est_str) == 0 ):
        continue
    try:
        point_estimates.append(float(est_str))
    except:
        sys.exit("Bad value for --point_estimates: \"%s\""%est_str)



def get_sort_name(row):
    special = row['special']
    if ( special == "expression" ):
        return special
    if ( special == "other_target/competitor" ):
        return special
    if ( special == "no_binding" ):
        return special
    if ( special == "naive" ):
        return special
    if ( special == "poly" ):
        return special
    if ( special == "neg_ctrl" ):
        return special
    if ( special == "protease" ):
        return special
    conc = row['concentration']
    target = ""
    if ( special == "other_target"):
        conc = row['special_num']
        target = "_other_target"
    elif ( special != "" ):
        return special

    avid = ""
    if ( row['avidity'] == "avid" ):
        avid = "_w/avid"
    name = "S%i_%inm%s%s"%(row['round'], conc, avid, target)

    
    return name



min_prob = -1e8

# Is linearly interpolating the binomial function the right thing to do?
# "right" is a very subjective term here
def my_binom_llh_cdf(x, n, p):
    
    x = x.clip(0, n)
    
    lb = scipy.stats.binom.logcdf(np.floor(x).astype(int), n=n, p=p).clip(min_prob, None)
    ub = scipy.stats.binom.logcdf(np.ceil(x).astype(int), n=n, p=p).clip(min_prob, None)
    frac = x - np.floor(x)
    
    return frac * ub + (1-frac) * lb


# Same concept as my_binom_llh_cdf except we come from the other side
# my_binom_llh_cdf + my_reverse_binom_llh_cdf ~= 1
def my_reverse_binom_llh_cdf(x, n, p):
    
    x = (n - x).clip(0, n)
    p = 1 - p
    
    lb = scipy.stats.binom.logcdf(np.floor(x).astype(int), n=n, p=p).clip(min_prob, None)
    ub = scipy.stats.binom.logcdf(np.ceil(x).astype(int), n=n, p=p).clip(min_prob, None)
    frac = x - np.floor(x)
    
    return frac * ub + (1-frac) * lb


# Null hypothesis is that x == n * p
# Creating a p-value that this is true by measuring the tails of the distribution
# Basically, we want to multiply the cdf by 2.
# However, you have to get the cdf from the side that's less than 0.5
# (If you try to go from the wrong side, you lose precision on numbers like 0.999999999)
#
# left_clipped and right_clipped are how we deal with data where the experimental p
#  was outside the significant bounds (i.e. it was lower than the min_selection_rate
#   or greater than collection_rate_saturation)
# If your p got clipped, your output llh is 0
def my_logp_equal(x, n, p, left_clipped=None, right_clipped=None):

    forward = my_binom_llh_cdf(x, n, p)
    reverse = my_reverse_binom_llh_cdf(x, n, p)
    
    if ( not right_clipped is None ):
        reverse[right_clipped] = -np.log(2)
        
    if ( not left_clipped is None ):
        forward[left_clipped] = -np.log(2)
    
    llh = np.minimum( forward, reverse ).clip(None, -np.log(2))
        
    llh += np.log(2) # like multiplying a probability by 2
    
    return llh



# For a given sorting round. Determine the highest concentration
# that was sorted at
def get_highest_conc_at_round(sorts, roundd, expression_parent, quiet=False):
    highest_conc = 0
    highest_bc = None
    for bc in sorts[(sorts['round'] == roundd) & (sorts['expression_parent'] == expression_parent)].index:
        if ( sorts.loc[bc]['special'] != "" ):
            continue
        concentration = sorts.loc[bc]['concentration']
        avid = sorts.loc[bc]['avidity'] == "avid"
        if ( avid ):
            concentration *= 100
        if ( concentration > highest_conc ):
            highest_bc = bc
            highest_conc = concentration
            
    if ( not highest_bc is None):
        return highest_bc
    
    
    for bc in sorts[ (sorts['round'] == roundd) & (sorts['expression_parent'] == expression_parent)].index:
        if ( sorts.loc[bc]['special'] != "other_target" ):
            continue
        concentration = sorts.loc[bc]['special_num']
        avid = sorts.loc[bc]['avidity'] == "avid"
        if ( avid ):
            concentration *= 100
        if ( concentration > highest_conc ):
            highest_bc = bc
            highest_conc = concentration
    
    if ( not highest_bc is None):
        if (not quiet):
            print("Warning: Using other_target for highest conc at round", roundd)
        return highest_bc
    
    return highest_bc
    
    

# definition of KD
def kd_to_frac(kd, conc):
    return conc / (kd + conc)


def assign_rounds_and_names(sorts):
    # Assing round by starting from expression and moving forwards
    sorts['expression_parent'] = ""
    sorts.loc[sorts['special'] == 'naive', 'round'] = -1
    expression_mask = sorts['special'] == 'expression'
    sorts.loc[expression_mask, 'round'] = 0
    sorts.loc[expression_mask, 'expression_parent'] = sorts[expression_mask]['pool_name']
    missing = True
    iterr = 0
    while ( missing ):
        missing = False
        iterr +=1
        for idx, row in sorts.iterrows():
            parent = row['parent_pool']
            if ( parent == "" ):
                continue
            if ( np.isnan(sorts.loc[parent]['round']) ):
                missing = True
                continue
            sorts.loc[idx, 'round'] = sorts.loc[parent]['round'] + 1
            if ( row['expression_parent'] == "" ):
                sorts.loc[idx, 'expression_parent'] = sorts.loc[parent]['expression_parent']
            
        if ( iterr > 100 ):
            sys.exit("Could not assign rounds. Parent structure does not make sense. All pools must be derived from expression"
                    +" (With the exception of a naive pool)")


    sorts['name'] = sorts.apply(get_sort_name, axis=1)



# Totally non-scientific function to ensure you formatted sorts.csv correctly
def pretty_print_sorts(sorts):
    def get_hier_size(hier):
        names = 0
        for name, child_hier in hier:
            names += 1 + get_hier_size(child_hier)
        return names
    def row_sorter(row):
        if ( row['special'] == "" ):
            if ( row['avidity'] == 'avid' ):
                return "%020.2f"%(1e8 - row['concentration'])
            else:
                return "%020.2f"%(1e9 + 1e8 - row['concentration'])
        if ( row['special'] == 'other_target'):
            return "%020.2f"%(1e10 + 1e8-row['special_num'])
        return row['special']
    def build_display_hierarchy(sorts, parent):
        children = sorts[sorts['parent_pool'] == parent ]
        if ( len(children) == 0 ):
            return []
        child_hiers = []
        child_hier_sizes = []
        for child in children['pool_name']:
            hier_child  = build_display_hierarchy(sorts, child)
            child_hiers.append([child, hier_child])
            child_hier_sizes.append(get_hier_size(hier_child))
        size_0s = [x for x, y in zip(child_hiers, child_hier_sizes) if y == 0]
        size_0s = sorted(size_0s, key=lambda item: row_sorter(sorts.loc[item[0]]) )
        this_hier = []
        for i in range(len(size_0s)):
            this_hier.append(size_0s[i])
        if ( len(size_0s) < len(child_hiers) ):
            size_other, other_sizes = list(zip(*[(x,y) for x, y in zip(child_hiers, child_hier_sizes) if y != 0]))
            arg_sorted = np.argsort(other_sizes)
            for i in arg_sorted:
                this_hier.append(size_other[i])
        return this_hier
    display_hierarchy = build_display_hierarchy(sorts, "")
    delta_indent = 2
    def display_it(sorts, hier, indent, fmt="%-40s"):
        longest = 0
        for name, inner in hier:
            row = sorts.loc[name]
            to_display = (" "*indent + row['name'])
            frac = "      " if np.isnan(row['collection_fraction']) else "%.4f"%row['collection_fraction']
            num = "        " if np.isnan(row['collected_cells']) else "%8i"%row['collected_cells']
            notes = "" if row['notes'] == "" else " -- " + row['notes']
            print("  %4s"%row['pool_name'] + fmt%to_display + "-- %s -- %s"%(frac, num) + notes )
            display_it(sorts, inner, indent+delta_indent, fmt)
    longest = (sorts['round']*delta_indent + sorts['name'].str.len()).max() + 4
    display_it(sorts, display_hierarchy, 2, fmt="%%-%is"%longest)










############################# MAIN ###########################################





sorts = pd.read_csv(args.__getattribute__("sorts.csv"), sep=",")
counts = pd.read_csv(args.__getattribute__("pooled_counts.list"), sep="\s+")
sorts.index = sorts['pool_name']

# change nan to empty string
sorts['parent_pool'] = sorts['parent_pool'].fillna("")
sorts['special'] = sorts['special'].fillna("")
sorts['avidity'] = sorts['avidity'].fillna("")
sorts['non_standard_gate'] = sorts['non_standard_gate'].fillna("")
sorts.loc[sorts['non_standard_gate'] != "", 'special'] += "_ns_gate"
sorts['notes'] = sorts['notes'].fillna("").astype(str)
sorts.loc[sorts['avidity'] == 'avi', 'avidity'] = 'avid'
sorts['round'] = np.nan
if ( "num_cells" in list(sorts) and "collected_cells" not in list(sorts)):
    sorts['collected_cells'] = sorts['num_cells']

if ( "frac_matching_seqs" in list(sorts) ):
    sorts['frac_matching_seqs'] = sorts['frac_matching_seqs'].fillna(1)
else:
    sorts['frac_matching_seqs'] = 1


if ( (sorts['special'] == 'expression').sum() == 0 ):
    sys.exit("You must have a value in sorts.csv where the special column == \"expression\"")

if ( ~np.all((sorts['avidity'] == 'avid') | (sorts['avidity'] == "") )):
    sys.exit("Valid choices for the avidity column are \"avid\" and \"\"")

if ( sorts['collection_fraction'].max() > 1 ):
    sys.exit("collection_fraction should be in fraction! ( 0.3333 ) not percent!")


# Give each sort a name and assign round 1, 2, 3, etc
assign_rounds_and_names(sorts)

print("")
print("Parsed input data")
# Sanity check your inputs
pretty_print_sorts(sorts)


# _frac is the fraction of each design collected by the cell sorter
# _cells is the number of cells for each design collected by the cell sorter
for bc in sorts.index:
    if ( counts[bc].sum() == 0 ):
        print("Error! There are 0 total counts for %s! Drop this from sorts.csv if you don't want it."%bc)
        sys.exit(1)
    norm_bc = counts[bc] / counts[bc].sum() * sorts.loc[bc, 'frac_matching_seqs']
    parent = sorts.loc[bc]['parent_pool']
    if ( parent != ""):
        counts[bc + "_frac"] = norm_bc / (counts[parent] / counts[parent].sum() * sorts.loc[parent, 'frac_matching_seqs'] ) \
                                    * sorts.loc[bc]['collection_fraction']

    counts[bc + "_cells"] = norm_bc * sorts.loc[bc]['collected_cells']


# Here we are back-calculating how many cells you sorted for R1
#  this ends up getting used for the doubly-transformed yeast calculations
for expression_parent in sorts['expression_parent'].unique():
    if ( expression_parent == "" ):
        continue
    
    r1_bc = get_highest_conc_at_round(sorts, 1, expression_parent)

    if ( r1_bc is None ):
        print("Expression sort %s doesn't have any children? Will throw an error later if this matters."%expression_parent)
        counts[expression_parent + "_cells"] = np.nan
        continue 

    r1_row = sorts.loc[r1_bc]
    
    counts[expression_parent + "_sorted_cells"] = \
                    (counts[expression_parent] / counts[expression_parent].sum() * sorts.loc[expression_parent, 'frac_matching_seqs'] * 
                                        (r1_row['collected_cells'] / r1_row['collection_fraction'])  )


print("")
print("Minimum cutoff for significant collected-cells and counts at sort (doubly transformed plasmid cuts)")
print("                                      (i.e. if counts are less than this, design might not be real)")
print("               cells      counts")

# This is the enrichment for the design with the most counts
#  We do it this way instead of looking for the most enriched design so we don't divide by 0 (or 1)
sorts['max_collection_enrich'] = 0
for bc in sorts.index:
    row = sorts.loc[bc]
    if ( row['round'] < 1 or np.isnan(row['round'])):
        continue
    expression_parent = row['expression_parent']
    counts_has_expression = counts[counts[expression_parent] != 0]
    biggest_i = counts_has_expression[bc].argmax()
    biggest_row = counts_has_expression.iloc[biggest_i]
    
    sorts.loc[bc, 'max_collection_enrich'] = (biggest_row[bc + "_cells"] / biggest_row[expression_parent + '_sorted_cells']
                                      * args.cell_copies_before_first_sort )
    row = sorts.loc[bc]

    print("%8s -- %8.0f -- %8.0f -- %s"%(bc, row['max_collection_enrich'], row['max_collection_enrich']*
            counts[bc].sum()/row['collected_cells'], row['name']))


# Call these the bounds. If you exceed them, you end up at 0 or inf

lowest_concentration = sorts[sorts['special'] == ""]['concentration'].min()
highest_concentration = sorts[sorts['special'] == ""]['concentration'].max()
low_kd = lowest_concentration / 10
high_kd = highest_concentration * 1e8

# The x-axis for KD calculations
kd_x = np.exp(np.arange( np.log(low_kd), np.log(high_kd), np.log(args.kd_fold_step)))


# if there's a 90% chance of getting sequenced, what's the expectation value of number
# of cells if you have 0 counts. This prevents divide-by-0
# With p_dna_is_sequenced == 0.9. min_allowable == 0.1
maybe_the_counts_were = np.arange(0, 1000)
p_the_counts_were = scipy.stats.binom.pmf(0, p=args.p_dna_is_sequenced, n=maybe_the_counts_were)
min_allowable_cells = (p_the_counts_were * maybe_the_counts_were).sum() / p_the_counts_were.sum()
    
avid_kd_llh = None
kd_llh = None

for avid in ['avid', ""]:
    print("")
    print("Calculating log-likelihoods: %s pass"%("avidity" if avid == "avid" else "non-avidity"))
    # llhs = []
    llh = np.zeros((len(counts), len(kd_x)))
    
    for bc in sorts.index:
        row = sorts.loc[bc]
        if ( not row['avidity'] == avid ):
            continue
        if ( row['special'] != "" ):
            continue

        print("   %s -- %s"%(row['pool_name'], row['name']))


        parent_bc = sorts.loc[bc]['parent_pool']


        # get parent and child cells ready
        # parent_cells entered the cell sorter and child_cells came out
        #
        # we need to make sure a few things are true:
        #  1. child_cells.sum() == sorts['collected_cells'] * sorts['frac_matching_seqs']
        #  2. parent_cells.sum() == sorts['collected_cells'] * sorts['frac_matching_seqs'] / sorts['collection_fraction']
        #  3. parent_cells.min() ~= min_allowable_cells
        
        child_cells = counts[bc + "_cells"].copy()
        total_parent_cells = row['collected_cells'] / row['collection_fraction']
        parent_cells = counts[parent_bc] / counts[parent_bc].sum() * sorts.loc[parent_bc, 'frac_matching_seqs'] * total_parent_cells

        # adjust 0s to be non-zero
        parent_cells = parent_cells.clip(min_allowable_cells, None)
        parent_cells = parent_cells / parent_cells.sum() * total_parent_cells * sorts.loc[parent_bc, 'frac_matching_seqs']

        assert(np.isclose(child_cells.sum() / row['frac_matching_seqs'], row['collected_cells']))
        

        # Get the probability of selection (p_sel) based on KD, min_sel_rate, and doubly-transformed floor

        p_sel_kd = kd_to_frac(kd_x, row['concentration'])
        
        # This is to account for doubly-transformed. We basically set the min_sel_rate
        #  high enough that you could get max_collection_enrich cells out as a null solution
        p_sel_doubly = ( row['max_collection_enrich'] / parent_cells.values ).clip(0, 1)
        
        # The p_sel value that is not reflective of binding
        p_sel_not_kd = p_sel_doubly.clip(args.min_selection_rate, 1)

        # Take the maximum of the two rates
        p_sel = np.maximum(p_sel_kd[None,:], p_sel_doubly[:,None])
    
        # If the predicted pass rate is > 0.20, and the actual pass rate is 0.20
        # we call it a perfect success because the max pass rate depends on other things
        right_clipped = p_sel > args.collection_rate_saturation
        
        # If the predicted pass rate is lower than the non_kd floor, we don't care how much lower it is
        # we give it the same value regardless of how much lower it is
        left_clipped = p_sel < p_sel_not_kd[:,None] * 1.05
        
        # Actually do the clipping. If p_sel_not_kd > collection_rate_saturation, the data is total garbage
        p_sel = p_sel.clip(p_sel_not_kd[:,None].clip(0, args.collection_rate_saturation),  args.collection_rate_saturation)
        

    
        # Now we start messing with the numbers to make the calculations more stable
        
        # limit the max number of cells to 1000 so we don't freak out the binomials
        to_max_parent_cells = (args.max_parent_cells / parent_cells).clip(0, 1)
        parent_cells *= to_max_parent_cells
        child_cells *= to_max_parent_cells
    
        # We're setting up the calculation so that it's:
        #  Given that there were n parent cells and the p(selection) was p
        #   what are the odds that we saw x cells (or worse) from this experiment
        #  Think of the classic z-test with null hypothesis x_0 == x
        p = p_sel
        n = np.ceil(parent_cells).values.astype(int)[:,None]
        x = child_cells.values[:,None]
        
    
        this_logp = my_logp_equal(x, n, p, 
                                  left_clipped=left_clipped,
                                  right_clipped=right_clipped
            )
        

        # llhs.append(this_logp)
        
        llh += this_logp



        # store this inside the loop so we don't store an empty to avid if there is none
        if ( avid == "avid"):
            avid_kd_llh = llh
        else:
            kd_llh = llh

if ( kd_llh is None ):
    if ( args.override_every_round_with_avid ):
        kd_llh = avid_kd_llh.copy()
    else:
        sys.exit("You did every sort with avidity? This script can't handle that situation.")


# This is the part of the script that is really bad. Bad from an assumption standpoint

low_conf = np.zeros(len(counts), np.bool)


avid_kd_includes_this = np.ones(kd_llh.shape, np.bool)
if ( not avid_kd_llh is None ):

    # The "best" logp point on each trace
    max_avid = avid_kd_llh.max(axis=-1)

    # assumption 7 and assumption 8
    # We boost the llh graph up and then find the two zero crossing points
    shifted_avid = avid_kd_llh - 2*max_avid[:,None] - args.confidence_interval_width
    avid_above_0 = (shifted_avid > 0)

    # this triggers when you have super low starting counts and the confidence interval includes
    #   all posible kd values
    avid_entirely_above_0 = avid_above_0.all(axis=-1) | ( avid_above_0[:,0] & avid_above_0[:,-1])

    # Find the first and last place where the shifted curve crosses 0
    avid_lb_i = avid_above_0.argmax(axis=-1)
    avid_ub_i = np.cumsum(avid_above_0.astype(int), axis=-1).argmax(axis=-1)

    avid_lb = kd_x[avid_lb_i]
    avid_ub = kd_x[avid_ub_i]

    avid_lb[avid_lb_i == 0] = 0
    avid_ub[avid_ub_i == 0] = 0
    avid_lb[avid_lb_i == len(kd_x)-1] = np.inf
    avid_ub[avid_ub_i == len(kd_x)-1] = np.inf

    avid_lb[avid_entirely_above_0] = 0
    avid_ub[avid_entirely_above_0] = np.inf

    # We're using the avidity to mask out regions of KD that are allowed
    # You aren't allowed to have a KD better than Avid, but you're allowed
    #  to be worse.

    avid_kd_includes_this &= np.cumsum(avid_above_0.astype(int), axis=-1).clip(0, 1).astype(bool)
    avid_kd_includes_this[:,-1] = True


    low_conf[avid_entirely_above_0] = True
else:
    avid_lb = np.zeros(len(counts), np.float)
    avid_ub = np.zeros(len(counts), np.float)
    avid_lb[:] = np.nan
    avid_ub[:] = np.nan


# The best logp on each trace
max_prob = kd_llh.max(axis=-1)


# assumption 7 and assumption 8
# We boost the llh graph up and then find the two zero crossing points
shifted_prob = kd_llh - 2*max_prob[:,None] - args.confidence_interval_width
prob_above_0 = shifted_prob > 0

# Mask of positions where avidity and non-avid sorts agree with each other
avid_and_prob_agree = prob_above_0 & avid_kd_includes_this

if ( args.mask_kd_with_avid_sorts ):
    if ( avid_kd_llh is None ):
        print("Warning! Can't do --mask_kd_with_avid_sorts if there are no avid sorts ")
    else:
        prob_above_0 = avid_and_prob_agree

# this triggers when you have super low starting counts and the confidence interval includes
#   all posible kd values
# Include the final | because after using avid_and_prob agree, you could end up with none of the
# probabilities above 0
prob_entirely_above_0 = ( prob_above_0.all(axis=-1) | ( prob_above_0[:,0] & prob_above_0[:,-1])
                        | ~prob_above_0.any(axis=-1) )


# Output value. True if there are no kd values with avidity and non-avidity agree
avid_doesnt_agree = ~avid_and_prob_agree.any()


# only assign low_conf from kd if there are no avid sorts
if ( avid_kd_llh is None ):
    low_conf[prob_entirely_above_0] = True

kd_lb_i = prob_above_0.argmax(axis=-1)
kd_ub_i = np.cumsum(prob_above_0.astype(int), axis=-1).argmax(axis=-1)

kd_lb = kd_x[kd_lb_i]
kd_ub = kd_x[kd_ub_i]

kd_lb[kd_lb_i == 0] = 0
kd_ub[kd_ub_i == 0] = 0
kd_lb[kd_lb_i == len(kd_x)-1] = np.inf
kd_ub[kd_ub_i == len(kd_x)-1] = np.inf

kd_lb[prob_entirely_above_0] = np.inf
kd_ub[prob_entirely_above_0] = np.inf

# One more check for low_conf. If the CI range is massive, we throw out datapoints
low_conf[ (kd_ub / kd_lb.clip(lowest_concentration, None) > 1000) & ( kd_lb < lowest_concentration ) ] = True


lowest_conc = np.zeros((len(kd_lb)))
highest_conc = np.zeros((len(kd_lb)))
lowest_conc[:] = lowest_concentration
highest_conc[:] = highest_concentration


records = {
    "kd_lb":kd_lb,
    "kd_ub":kd_ub,
    "low_conf":low_conf,
    "avid_doesnt_agree":avid_doesnt_agree,
    "avid_lb":avid_lb,
    "avid_ub":avid_ub,
    "lowest_conc":lowest_conc,
    "highest_conc":highest_conc
}




################# Point estimate stuff #######################

# find the best concentration for the point estimate estimation
# Best is exact match
# Next best is something at a lower concentration
# Worst is at a higher concentration
# Also, prefer later rounds
def find_best_conc(sorts, search_conc):

    for allow_avid in [False, True]:

        for bc in sorted(sorts.index, key=lambda x: -sorts.loc[x]['round']):
            row = sorts.loc[bc]
            if ( row['special'] != "" ):
                continue
            if ( (row['avidity'] == "avid") != allow_avid ):
                continue
            if ( np.isnan(row['collection_fraction'])):
                continue

            conc = row['concentration']

            if ( np.isclose(conc, search_conc, rtol=0.05) ):
                return bc
        best_bc = None
        closest = 0

        #better to go under than over
        for bc in sorted(sorts.index, key=lambda x: -sorts.loc[x]['round']):
            row = sorts.loc[bc]
            if ( row['special'] != "" ):
                continue
            if ( (row['avidity'] == "avid") != allow_avid ):
                continue
            if ( np.isnan(row['collection_fraction'])):
                continue

            conc = row['concentration']

            if ( conc > search_conc ):
                continue

            if ( conc*4 < search_conc ):
                continue

            if (conc > closest):
                closest = conc
                best_bc = bc

        if (not best_bc is None):
            return best_bc


        closest = 100000000000

        #forced to go over
        for bc in sorted(sorts.index, key=lambda x: -sorts.loc[x]['round']):
            row = sorts.loc[bc]
            if ( row['special'] != "" ):
                continue
            if ( (row['avidity'] == "avid") != allow_avid ):
                continue
            if ( np.isnan(row['collection_fraction'])):
                continue

            conc = row['concentration']

            if ( conc < search_conc ):
                continue

            if ( conc/4 > search_conc ):
                continue

            if (conc < closest):
                closest = conc
                best_bc = bc


        if (not best_bc is None):
            return best_bc

    return None




for kd_est in point_estimates:
    use_kd_est = kd_est

    # Ok, we want to find the sorting concentration where someone with kd_est
    # will have a 0.20 collection rate

    # collection_rate = conc / (kd + conc)
    # (kd + conc) * collection_rate = conc
    # kd * collection_rate = conc * ( 1 - collection_rate )
    # conc = kd * collection_rate / ( 1 - collection_rate )

    search_conc = use_kd_est * args.collection_rate_saturation / ( 1 - args.collection_rate_saturation )


    best_bc = find_best_conc(sorts, search_conc)

    if ( best_bc is None ):
        print("Warning! You don't have any sorts where you can "+
                        "point-estimate a KD of %.1f. Was looking for a sort around %.1f nM"%(kd_est, search_conc))
        continue

    row = sorts.loc[best_bc]

    label = "%i"%(kd_est)


    conc = row['concentration']

    if ( conc > search_conc * 1.05 ):

        use_kd_est = use_kd_est * conc / search_conc
        label = "%i"%(round(use_kd_est))


        print("Warning! Could not find a good sort for point-estimating a KD of %i. Was looking for a "%(kd_est)+
                "sort around %.1f nM or lower. Will use %s. However, this requires changing the point"%(search_conc, row['name'])+
                " estimate to %s nM."%(label))

    if ( row['avidity'] == "avid" ):
        label += "w"

        print("Warning! Could not find a good sort for point-estimating a KD of %i. Was looking for a "%(kd_est)+
                "sort around %.1f nM or lower. Will use %s. However, note that this sort was performed with avidity"%(search_conc, row['name'])+
                " and these results cannot be compared with sorts done without avidity.")


    print("")
    print("Estimating %s nM binders from %s."%(label, row['name']))

    expected_frac = conc / ( use_kd_est + conc )
    expected_frac = min( expected_frac, args.collection_rate_saturation )

    enough_counts = counts[best_bc + "_cells"] > row['max_collection_enrich']
    enough_frac = counts[best_bc + "_frac"] > expected_frac

    its_a_binder = enough_counts & enough_frac
        
    records['binder_%s_nm'%label] = its_a_binder



if ( args.target != "" ):
    records['target'] = args.target


# sketch_kd stuff
def norm(x):
    return x / x.sum()

cool_sorts = sorts[(sorts['special'] == "") & (sorts['avidity'] != "avid")]['pool_name']


records['description'] = counts['description']


df = pd.DataFrame(records)
df = df.merge(counts, 'inner', 'description')



df['sketch_kd'] = False
df['sketch_kd_level'] = 1
for bc in cool_sorts:
    expression_sort = sorts.loc[bc]['expression_parent']
    enrich = norm(df[bc]) / norm(df[expression_sort])
    conc = sorts.loc[bc]['concentration']


    m = (df['kd_ub'] < 10000) & ~(df['low_conf'])
    to_plot = df[m]
    enrich_to_plot = enrich[m]

    to_plot = to_plot.sort_values("kd_ub", ascending=False)
    max_at = np.maximum.accumulate(enrich_to_plot)
    mask = to_plot['kd_ub'] < conc*10
    df.loc[to_plot[mask].index, 'sketch_kd_level'] = np.minimum(to_plot[mask]['sketch_kd_level'], enrich_to_plot[mask]/max_at[mask])
    df.loc[df['sketch_kd_level'] < args.sketch_kd_threshold, 'sketch_kd'] = True





df.to_csv(args.output_file, sep=" ", index=None, na_rep="NaN")







