#!/usr/bin/env python

# Estimate the binding SC50 from yeast-surface-display sorting experiments
#  (the kind performed in Cao et. al 2021).
#
# Brian Coventry 2021, 2025
#   -- Based in concept on the EC50 calculator from Rocklin 2017, Science
#
# This is the clean, single-method predictor: it implements ONLY "Method N20", the
# configuration selected at the end of the 2025 method search. Every modeling choice below
# is fixed; there are no method-selection flags. See `what_the_best_predictor_looks_like.md`
# for the full justification of each choice. In one sentence: a real binomial likelihood
# (peaked, so mid-range binders localize), a posterior-median point estimate (smooth and
# replicate-stable), a Bayesian credible interval tempered by a replicate-fit overdispersion
# factor (honestly calibrated), and an identifiability gate (never report an affinity beyond
# the range you actually tested).
#
# SC50 -- the concentration of target at which 50% of a design's yeast are collected by the
#         cell sorter (lower = tighter binder).


# ============================ What this script computes ============================
#
# For each design we treat the cell sorter as a binomial process: given the parent-pool cell
# count and the per-SC50 expected collection probability, how likely is the observed child
# count?  Summing the binomial log-PMF over all titration sorts gives a log-likelihood profile
# over a grid of candidate SC50 values.  From that profile we report:
#
#   - sc50_est : the posterior MEDIAN under a uniform prior on log(SC50), built from the
#                UNtempered likelihood (full mid-range localization).  This is the headline
#                number.  It is forced to inf (non-binder) when the identifiability gate fires.
#   - sc50_lb / sc50_ub : an equal-tailed Bayesian CREDIBLE interval at --ci_level, built from
#                the likelihood TEMPERED by --likelihood_temp_ci.  The tempering is a
#                replicate-fit overdispersion correction (the raw binomial interval is
#                overconfident ~2.5-6x); it touches ONLY the interval, p_bind and low_conf --
#                never sc50_est.
#   - sc50_unmeasurable : True when sc50_est lands beyond --identifiability_gate x the highest
#                tested (titration) concentration.  Such an estimate is pure counting noise
#                (flips between replicates); we censor sc50_est / lb / ub to inf and flag it.
#                Treat True as "non-binder / off-scale".
#   - sc50_unmeasurable_cutoff : the concentration cutoff used for the gate above
#                (highest titration concentration x --identifiability_gate); inf if disabled.
#   - p_bind : P(SC50 < identifiability_gate × highest_titration_concentration), a soft
#                "is this a binder within the measurable range" probability from the (tempered)
#                posterior over titration sorts.  The threshold equals sc50_unmeasurable_cutoff,
#                so p_bind and the hard gate use the same trust ceiling.
#   - p_bind_threshold : the SC50 threshold used for p_bind (= sc50_unmeasurable_cutoff).
#   - rel_enrich_titr : the passenger-plasmid score.  How far this design's enrichment from the
#                naive expression pool to the highest-concentration titration sort sits BELOW the
#                maximum enrichment among designs with the same apparent SC50 (real/fraction units).
#                Values near 1.0 => enriches as well as same-affinity peers; near 0 => far less
#                enriched => likely a doubly-transformed passenger.  A design alone in its SC50
#                window scores 1.0 (not flagged).
#                This is a SOFT score, not a reliable hard gate -- it does not alter sc50_est.
#   - suspected_passenger : provisional bool, rel_enrich_titr < --suspected_passenger_threshold.
#                Default threshold 0.05 (= 1/20): enriches at less than 5% of the cohort maximum.
#                Reliable only in SSM/designed-pool experiments; prefer the continuous rel_enrich_titr.
#   - low_conf : the data had too few counts to say anything (interval spans the whole grid,
#                or is enormous while pinned at the floor).
#
# NOTE on avidity: avidity sorts are EXCLUDED from the fit.  Many input files use a "fake avid"
# label on prior enrichment rounds to keep them out of the titration; that idiom still works.
# This script does not produce the legacy avid_lb / avid_ub / avid_doesnt_agree columns.


# ================================ Input data format ================================
#
# ./estimate_sc50_from_ngs.py pooled_counts.list sorts.csv
#
# pooled_counts.list -- experimental counts (whitespace separated)
#   -- one column per pool (column titles match pool_name in sorts.csv)
#   -- a "description" column naming each design
#   -- absolute column sums do not matter (sorts.csv collected_cells handles normalization)
#   -- Example header: bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8 description
#
# sorts.csv -- details of the sorting experiment (comma separated). Required columns:
#   -- pool_name        : name of this pool in pooled_counts.list
#   -- parent_pool      : the pool that entered the sorter to produce this one (blank if it is
#                         not in pooled_counts.list / is a root)
#   -- concentration    : target concentration (historically nM; units are arbitrary)
#   -- avidity          : { avid, "" }.  "avid" rows are excluded from the fit (see note above)
#   -- special          : { expression, "" , ... }.  You MUST have exactly the expression sort
#                         marked "expression".  Any other non-empty value means "not a normal
#                         titration sort" and the row is ignored by the fit.
#   -- special_num      : unused; leave blank
#   -- collection_fraction : fraction of cells collected by your final gate (0-1, not percent)
#   -- collected_cells  : how many cells you collected (estimating is fine). "num_cells" is
#                         accepted as a synonym.
#   -- non_standard_gate: leave blank; any value here rejects the row
#   -- notes            : free text
#   -- frac_matching_seqs : optional; fraction of reads matching expected sequences (default 1)
#
# Example header:
#   pool_name,parent_pool,concentration,avidity,special,special_num,collection_fraction,collected_cells,non_standard_gate,notes


# ================================ Output columns ===================================
#
#   sc50_est, sc50_lb, sc50_ub, sc50_unmeasurable, sc50_unmeasurable_cutoff, p_bind,
#   rel_enrich_titr, suspected_passenger, low_conf, lowest_conc, highest_conc,
#   [target], description
#
# With --dump_round_information, the per-pool diagnostic columns are also emitted:
#   {pool}            : counts of this design in this pool
#   {pool}_cells      : cells calculated to have been collected in this pool
#   {pool}_frac       : fraction of this design's cells collected in this pool
#   {pool}_sorted_cells : back-calculated cells sorted for the expression sort


# ================================ Key assumptions ==================================
#
#  1. Each design expresses perfectly such that 100% of cells would be collected on an
#     expression-only sort.
#     -- Why: Avoids co-fitting expression and SC50 (two unknowns per design).
#     -- Effect: Higher-expressing designs get a modest SC50 advantage (~10x for 75% vs 25%
#                expression). Acceptable given the difficulty of separating the two effects.
#
#  2. Collection rates >= --collection_rate_saturation (default 0.20) are all treated as
#     saturating (equivalent to 20%).
#     -- Why: Assumption 1 is not perfectly true; some designs express poorly. Anything over
#             20% is "saturating" so poorly-expressing designs don't confound the fit.
#     -- Effect: SC50s are fit more accurately near the saturation end.
#
#  3/4. Doubly-transformed yeast cells are the dominant noise source, handled by flooring
#     the expected collection rate at the most-copied-cell-since-transformation level
#     (max_collection_enrich) per pool.
#     -- Why: A doubly-transformed yeast carries non-binding DNA through successive sorts.
#             If we assume only one doubly-transformed cell per sequence (the same DNA
#             didn't get lucky twice), the worst-case cell count equals the number of
#             copies of that initial cell -- we take the per-pool maximum as the floor.
#     -- Effect: Doubly-transformed passenger contamination is largely suppressed. Designs
#                with very low transformation efficiency may be called non-binders.
#
#  5/6. The sorter is modeled as a binomial process; non-integer child counts are linearly
#     interpolated into the PMF.
#     -- Why: Allows computing P(SC50 == X | data) directly; sharply peaked so mid-range
#             binders localize well.
#     -- Effect: More accurate point estimates than the legacy two-sided CDF approach;
#                requires the overdispersion temperature (--likelihood_temp_ci) for honest CIs.


# ========================= Can you use this script? ================================
#
# Critical experimental requirements. If these are met, the analysis is likely valid.
# If they are not, results may be misleading:
#
# 1. You are performing a sorting analysis where you collect the BINDING population.
#
# 2. Your gate is set such that nearly NONE of the cells are collected when your
#    binder target is not present (negative control). See Figure 1 below.
#
# 3. Your gate is set to collect EVERYTHING that shows signal above the negative
#    control. See Figure 1 below.
#
# 4. You are fully NGS-sequencing every pool of interest.
#
# 5. Your selection step is NON-COMPETITIVE (strong binders cannot crowd out weak
#    binders). Running experiments in a column, taking the top 1% by FACS, or
#    sorting at 1 pM in a tiny volume are not compatible with this script.
#
# 6. You are recording the FRACTION of cells collected (collection_fraction in sorts.csv).
#
# 7. Each design expresses well enough that ~20% would be collected if you sorted on
#    expression alone (--collection_rate_saturation default). Designs that express
#    much lower will appear as weaker binders than they are.
#
#
# Figure 1. How to set your gate after running a negative control
#           (negative control = fluorescently labeled target is absent)
#
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
#    #ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd/
#                                  Expression Axis (log)
#
# https://www.text-image.com/convert/pic2ascii.cgi


import sys
import argparse
import scipy.stats
import scipy.special
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(
    description="Estimate binding SC50 from YSD sorting NGS data (Method N20). "
                "Open this file in a text editor for the full input/output format.")
parser.add_argument("pooled_counts.list", type=str)
parser.add_argument("sorts.csv", type=str)
parser.add_argument("--output_file", default="affinity_estimates.sc", type=str,
                    help="The output file.")
parser.add_argument("--dump_round_information", action="store_true",
                    help="Also emit the per-pool diagnostic columns (counts, cells, frac, "
                         "sorted_cells).")
parser.add_argument("--target", default="", type=str,
                    help="Add a constant 'target' column to the output.")
parser.add_argument("--cell_copies_before_first_sort", type=float, default=1,
                    help="Between transformation and the first sort with target, how many "
                         "copies of each cell happened? (i.e. 2**cell_divisions)")

# --- Foundational assumptions (defaults are the validated values; change with care) ---
parser.add_argument("--min_selection_rate", default=0.0001, type=float,
                    help="Expected background (non-binding) collection rate.")
parser.add_argument("--collection_rate_saturation", default=0.20, type=float,
                    help="Collection rates at or above this are treated as saturating.")
parser.add_argument("--max_parent_cells", default=1000, type=float,
                    help="Cap on parent-pool cells used per design. Counts above this are "
                         "scaled down so the binomial does not become absurdly sharp for "
                         "high-count designs. Matters MORE under the binomial likelihood and "
                         "is part of why the CI temperature is tractable -- keep it.")
parser.add_argument("--sc50_fold_step", default=1.1, type=float,
                    help="SC50 grid resolution (fold step between grid points).")
parser.add_argument("--p_dna_is_sequenced", default=0.90, type=float,
                    help="Minor parameter; only used to avoid divide-by-zero for 0-count cells.")

# --- The interval / threshold knobs ---
parser.add_argument("--ci_level", default=0.9973, type=float,
                    help="Credible-interval mass (default 0.9973 ~ 3 sigma). Calibrated by "
                         "--likelihood_temp_ci.")
parser.add_argument("--likelihood_temp_ci", default=0.05, type=float,
                    help="Overdispersion temperature applied to the likelihood for the CI / "
                         "p_bind / low_conf ONLY (never the point estimate). The raw binomial "
                         "interval is overconfident ~2.5-6x; std(z) ~ sqrt(temp), so ~0.05 "
                         "calibrates it. This value is FIT FROM REPLICATE DATA -- re-fit it "
                         "from replicate std(z) if the assay's noise changes; do not hand-tune.")
parser.add_argument("--identifiability_gate", default=20.0, type=float,
                    help="Trust ceiling. If sc50_est exceeds (highest tested concentration x "
                         "this), censor it to inf and set sc50_unmeasurable=True. A policy "
                         "choice about how far past the titration range you trust extrapolation "
                         "(default 20x). Set 0 to disable.")
parser.add_argument("--suspected_passenger_threshold", default=0.05, type=float,
                    help="Fractional threshold for suspected_passenger: flags designs where "
                         "exp(rel_enrich_titr) < this value, i.e. the design's enrichment is "
                         "less than (threshold x cohort reference). "
                         "Default 0.05 (= 1/20: enriches at less than 5%% of cohort reference). Meaningful "
                         "only in SSM/designed-pool experiments where the cohort is rich in real "
                         "binders; unreliable in random library screens. PROVISIONAL -- treat "
                         "rel_enrich_titr as the primary output and tune per use case.")

args = parser.parse_args()


min_prob = -1e8


# ---------------------------------------------------------------------------------
# Helpers (ported unchanged from estimate_affinity_from_ngs_sc50er.py)
# ---------------------------------------------------------------------------------

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


# The real binomial log-likelihood log P(x | n, p). It is sharply peaked (unlike a two-sided
# goodness-of-fit p-value), so it localizes mid-range binders. The (non-integer) child count x
# is linearly interpolated into the PMF the same way the legacy CDF helpers did.
def my_binom_logpmf(x, n, p):
    x = x.clip(0, n)
    lo = scipy.stats.binom.logpmf(np.floor(x).astype(int), n=n, p=p)
    hi = scipy.stats.binom.logpmf(np.ceil(x).astype(int), n=n, p=p)
    frac = x - np.floor(x)
    out = frac * hi + (1 - frac) * lo
    # logpmf can be -inf when p is exactly 0/1 and x disagrees; keep it finite so the
    # softmax/median downstream don't choke.
    return np.clip(out, min_prob, None)


# For a given sorting round, find the pool sorted at the highest concentration (used for the
# doubly-transformed back-calculation).
def get_highest_conc_at_round(sorts, roundd, expression_parent, quiet=False):
    highest_conc = 0
    highest_bc = None
    for bc in sorts[(sorts['round'] == roundd) & (sorts['expression_parent'] == expression_parent)].index:
        if ( sorts.loc[bc]['special'] != "" ):
            continue
        concentration = sorts.loc[bc]['concentration']
        if ( sorts.loc[bc]['avidity'] == "avid" ):
            concentration *= 100
        if ( concentration > highest_conc ):
            highest_bc = bc
            highest_conc = concentration

    if ( not highest_bc is None):
        return highest_bc

    for bc in sorts[(sorts['round'] == roundd) & (sorts['expression_parent'] == expression_parent)].index:
        if ( sorts.loc[bc]['special'] != "other_target" ):
            continue
        concentration = sorts.loc[bc]['special_num']
        if ( sorts.loc[bc]['avidity'] == "avid" ):
            concentration *= 100
        if ( concentration > highest_conc ):
            highest_bc = bc
            highest_conc = concentration

    if ( not highest_bc is None):
        if (not quiet):
            print("Warning: Using other_target for highest conc at round", roundd)
        return highest_bc

    return highest_bc


# definition of SC50
def sc50_to_frac(sc50, conc):
    return conc / (sc50 + conc)


# Cohort-relative enrichment (the rel_enrich_titr passenger score). For each design with a
# finite apparent SC50, return log(enrichment) minus the MAXIMUM log(enrichment) among designs
# whose apparent SC50 is within `window` log-units (the "same affinity" cohort).
# Very negative => this design is far less enriched than same-affinity peers => likely a
# doubly-transformed passenger (diluted by its non-binding single-transformant copies).
# A design alone in its SC50 window scores 0 (at the reference) and is not flagged as a
# passenger — this correctly handles isolated real binders in passenger-dominated experiments.
# Vectorized via unique (lo,hi) pairs so the max is computed once per unique window.
def cohort_relative_enrich(log_enrich, log_sc50, valid, window=np.log(3.0)):
    n = len(log_enrich)
    out = np.full(n, np.nan)
    idx = np.where(valid & np.isfinite(log_enrich) & np.isfinite(log_sc50))[0]
    if len(idx) == 0:
        return out
    ls = log_sc50[idx]
    le = log_enrich[idx]
    order = np.argsort(ls)
    ls_s = ls[order]
    le_s = le[order]
    lo = np.searchsorted(ls_s, ls_s - window, side="left")
    hi = np.searchsorted(ls_s, ls_s + window, side="right")
    window_ids, inverse = np.unique(np.stack([lo, hi], axis=1), axis=0, return_inverse=True)
    ref_values = np.array([le_s[l:h].max() for l, h in window_ids])
    out[idx[order]] = le_s - ref_values[inverse]
    return out


def assign_rounds_and_names(sorts):
    # Assign round by starting from expression and moving forwards
    sorts['expression_parent'] = ""
    sorts.loc[sorts['special'] == 'naive', 'round'] = -1
    expression_mask = sorts['special'] == 'expression'
    sorts.loc[expression_mask, 'round'] = 0
    sorts.loc[expression_mask, 'expression_parent'] = sorts[expression_mask]['pool_name']
    missing = True
    iterr = 0
    while ( missing ):
        missing = False
        iterr += 1
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
            sys.exit("Could not assign rounds. Parent structure does not make sense. All pools "
                     "must be derived from expression (with the exception of a naive pool)")

    sorts['name'] = sorts.apply(get_sort_name, axis=1)


# Totally non-scientific function to sanity-check that you formatted sorts.csv correctly
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
            return "%020.2f"%(1e10 + 1e8 - row['special_num'])
        return row['special']
    def build_display_hierarchy(sorts, parent):
        children = sorts[sorts['parent_pool'] == parent ]
        if ( len(children) == 0 ):
            return []
        child_hiers = []
        child_hier_sizes = []
        for child in children['pool_name']:
            hier_child = build_display_hierarchy(sorts, child)
            child_hiers.append([child, hier_child])
            child_hier_sizes.append(get_hier_size(hier_child))
        size_0s = [x for x, y in zip(child_hiers, child_hier_sizes) if y == 0]
        size_0s = sorted(size_0s, key=lambda item: row_sorter(sorts.loc[item[0]]) )
        this_hier = []
        for i in range(len(size_0s)):
            this_hier.append(size_0s[i])
        if ( len(size_0s) < len(child_hiers) ):
            size_other, other_sizes = list(zip(*[(x, y) for x, y in zip(child_hiers, child_hier_sizes) if y != 0]))
            arg_sorted = np.argsort(other_sizes)
            for i in arg_sorted:
                this_hier.append(size_other[i])
        return this_hier
    display_hierarchy = build_display_hierarchy(sorts, "")
    delta_indent = 2
    def display_it(sorts, hier, indent, fmt="%-40s"):
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
pretty_print_sorts(sorts)


# The sorts that go into the SC50 fit: normal titration sorts only -- non-special AND non-avid.
# (Avid rows are excluded; the "fake avid" idiom for ignoring prior rounds still works.)
fit_sort_mask = (sorts['special'] == "") & (sorts['avidity'] != "avid")
if ( fit_sort_mask.sum() == 0 ):
    sys.exit("No titration sorts to fit (every non-special sort is marked avid). This script "
             "excludes avid sorts from the fit and cannot analyze an all-avidity experiment.")


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


# Back-calculate how many cells were sorted for R1 (used for doubly-transformed calculations)
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

# Doubly-transformed floor: the enrichment for the design with the most counts in each sort.
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


# Concentration bounds. lowest is from fit sorts only; highest uses all non-special sorts
# (including avid) so the identifiability gate reflects the full tested range — an avid sort at
# 1000 nM does constitute a "tested at 1000 nM" event even though it is excluded from the fit.
lowest_concentration = sorts[fit_sort_mask]['concentration'].min()
highest_concentration = sorts[sorts['special'] == ""]['concentration'].max()
low_sc50 = lowest_concentration / 10
high_sc50 = highest_concentration * 1e8

# The x-axis for SC50 calculations (log-spaced -> uniform prior on log(SC50))
sc50_x = np.exp(np.arange( np.log(low_sc50), np.log(high_sc50), np.log(args.sc50_fold_step)))

# Expectation value of the number of cells given 0 counts, to prevent divide-by-0.
maybe_the_counts_were = np.arange(0, 1000)
p_the_counts_were = scipy.stats.binom.pmf(0, p=args.p_dna_is_sequenced, n=maybe_the_counts_were)
min_allowable_cells = (p_the_counts_were * maybe_the_counts_were).sum() / p_the_counts_were.sum()


# ---------------------------------------------------------------------------------
# Likelihood: sum the binomial log-PMF over the titration sorts.
# ---------------------------------------------------------------------------------
print("")
print("Calculating binomial log-likelihoods over titration sorts")

sc50_llh = np.zeros((len(counts), len(sc50_x)))
for bc in sorts[fit_sort_mask].index:
    row = sorts.loc[bc]
    print("   %s -- %s"%(row['pool_name'], row['name']))

    parent_bc = row['parent_pool']

    # parent_cells entered the sorter; child_cells came out.
    child_cells = counts[bc + "_cells"].copy()
    total_parent_cells = row['collected_cells'] / row['collection_fraction']
    parent_cells = counts[parent_bc] / counts[parent_bc].sum() * sorts.loc[parent_bc, 'frac_matching_seqs'] * total_parent_cells
    parent_cells = parent_cells.clip(min_allowable_cells, None)
    parent_cells = parent_cells / parent_cells.sum() * total_parent_cells * sorts.loc[parent_bc, 'frac_matching_seqs']

    assert(np.isclose(child_cells.sum() / row['frac_matching_seqs'], row['collected_cells']))

    # Probability of selection per candidate SC50, plus the doubly-transformed floor.
    p_sel_sc50 = sc50_to_frac(sc50_x, row['concentration'])
    p_sel_doubly = ( row['max_collection_enrich'] / parent_cells.values ).clip(0, 1)
    p_sel_not_sc50 = p_sel_doubly.clip(args.min_selection_rate, 1)

    p_sel = np.maximum(p_sel_sc50[None, :], p_sel_doubly[:, None])
    # Clip to [floor, saturation]: below-floor (non-binder) and above-saturation (tight) SC50
    # regions become flat (no information), while an enriched design IS penalized under the
    # floor rate so the likelihood is peaked at the matching SC50.
    p_sel = p_sel.clip(p_sel_not_sc50[:, None].clip(0, args.collection_rate_saturation), args.collection_rate_saturation)

    # Cap parent cells so the binomial doesn't get absurdly sharp for high-count designs.
    to_max_parent_cells = (args.max_parent_cells / parent_cells).clip(0, 1)
    parent_cells *= to_max_parent_cells
    child_cells *= to_max_parent_cells

    n = np.ceil(parent_cells).values.astype(int)[:, None]
    x = child_cells.values[:, None]
    sc50_llh += my_binom_logpmf(x, n, p_sel)


# ---------------------------------------------------------------------------------
# Posterior bookkeeping.
#   - The point estimate (sc50_est) uses the UNtempered likelihood (full localization).
#   - The interval / p_bind / low_conf use the likelihood tempered by --likelihood_temp_ci.
# ---------------------------------------------------------------------------------
temp_ci = args.likelihood_temp_ci
sc50_llh_ci = sc50_llh * temp_ci

low_conf = np.zeros(len(counts), bool)

# p_bind threshold is the same trust ceiling as the identifiability gate so both use a
# consistent definition of "measurable range."
p_bind_threshold = highest_concentration * args.identifiability_gate if args.identifiability_gate > 0 else np.inf

# p_bind = P(SC50 < threshold) from the tempered posterior (uniform prior on log SC50).
_thresh_i = min(int(np.searchsorted(sc50_x, p_bind_threshold)), len(sc50_x))
_log_norm = scipy.special.logsumexp(sc50_llh_ci, axis=-1)
if _thresh_i > 0:
    _log_below = scipy.special.logsumexp(sc50_llh_ci[:, :_thresh_i], axis=-1)
    p_bind = np.exp(_log_below - _log_norm)
else:
    p_bind = np.zeros(len(sc50_llh_ci))

# Equal-tailed Bayesian credible interval from the tempered posterior.
_ci_max = sc50_llh_ci.max(axis=-1)
_post = np.exp(sc50_llh_ci - _ci_max[:, None])
_post /= _post.sum(axis=-1, keepdims=True).clip(1e-300, None)
_post_cdf = np.cumsum(_post, axis=-1)

_alpha = 1.0 - args.ci_level
sc50_lb_i = (_post_cdf >= _alpha / 2.0).argmax(axis=-1)
_ub_reached = _post_cdf >= 1.0 - _alpha / 2.0
sc50_ub_i = _ub_reached.argmax(axis=-1)
sc50_ub_i[~_ub_reached.any(axis=-1)] = len(sc50_x) - 1

_grid_i = np.arange(len(sc50_x))[None, :]
prob_above_0 = (_grid_i >= sc50_lb_i[:, None]) & (_grid_i <= sc50_ub_i[:, None])

# Triggers when starting counts are so low the interval includes ~all SC50 values.
prob_entirely_above_0 = ( prob_above_0.all(axis=-1) | (prob_above_0[:, 0] & prob_above_0[:, -1])
                          | ~prob_above_0.any(axis=-1) )
low_conf[prob_entirely_above_0] = True

sc50_lb = sc50_x[sc50_lb_i]
sc50_ub = sc50_x[sc50_ub_i]
sc50_lb[sc50_lb_i == 0] = 0
sc50_ub[sc50_ub_i == 0] = 0
sc50_lb[sc50_lb_i == len(sc50_x)-1] = np.inf
sc50_ub[sc50_ub_i == len(sc50_x)-1] = np.inf
sc50_lb[prob_entirely_above_0] = np.inf
sc50_ub[prob_entirely_above_0] = np.inf

# Point estimate: posterior MEDIAN from the UNtempered likelihood.
_pp = np.exp(sc50_llh - sc50_llh.max(axis=-1, keepdims=True))
_pp /= _pp.sum(axis=-1, keepdims=True).clip(1e-300, None)
_pp_cdf = np.cumsum(_pp, axis=-1)
sc50_median_i = (_pp_cdf >= 0.5).argmax(axis=-1)
sc50_est = sc50_x[sc50_median_i].astype(float)

# If the CI range is massive while pinned at the floor, throw out the datapoint.
# Use the UN-tempered posterior here: likelihood_temp_ci widens the credible interval
# intentionally (calibrated to replicate noise), so the tempered lb/ub naturally fall
# below lowest_concentration even for real in-range binders.  The un-tempered posterior
# correctly reflects information content: for truly uninformative designs it also spans
# the full grid, but for real binders the un-tempered lb stays within the titration range.
_untempered_lb_i = (_pp_cdf >= _alpha / 2.0).argmax(axis=-1)
_untempered_ub_reached = _pp_cdf >= 1.0 - _alpha / 2.0
_untempered_ub_i = _untempered_ub_reached.argmax(axis=-1)
_untempered_ub_i[~_untempered_ub_reached.any(axis=-1)] = len(sc50_x) - 1
_untempered_lb = sc50_x[_untempered_lb_i].astype(float)
_untempered_ub = sc50_x[_untempered_ub_i].astype(float)
_untempered_lb[_untempered_lb_i == 0] = 0.0
_untempered_ub[_untempered_ub_i == len(sc50_x) - 1] = np.inf
low_conf[(_untempered_ub / _untempered_lb.clip(lowest_concentration, None) > 1000) & (_untempered_lb < lowest_concentration)] = True
# Raw (finite, uncensored) median grid value -- the apparent SC50 used to form the passenger
# cohort below. Kept before the non-binder/gate inf-censoring so cohorts span the full
# population (matching the validated rel_enrich_titr, which bins every design by apparent SC50).
sc50_apparent = sc50_est.copy()
sc50_est[sc50_median_i == len(sc50_x) - 1] = np.inf
sc50_est[prob_entirely_above_0] = np.inf


# Identifiability gate: censor estimates beyond the trustable extrapolation range to inf.
sc50_unmeasurable = np.zeros(len(sc50_est), dtype=bool)
gate_cutoff = highest_concentration * args.identifiability_gate if args.identifiability_gate > 0 else np.inf
if args.identifiability_gate > 0:
    sc50_unmeasurable = np.isfinite(sc50_est) & (sc50_est > gate_cutoff)
    print("")
    print("Identifiability gate at est > %.3g (%.3g top sort x %g): %d designs censored to inf "
          "(unmeasurable beyond titration range)"
          % (gate_cutoff, highest_concentration, args.identifiability_gate,
             int(sc50_unmeasurable.sum())))
    sc50_est[sc50_unmeasurable] = np.inf
    sc50_lb[sc50_unmeasurable] = np.inf
    sc50_ub[sc50_unmeasurable] = np.inf


# ---------------------------------------------------------------------------------
# Passenger score: rel_enrich_titr (the validated passenger-plasmid detector).
# Enrichment from the naive expression pool to the highest-concentration titration sort,
# expressed RELATIVE to the same-apparent-SC50 cohort. A doubly-transformed passenger -- diluted
# by its non-binding single-transformant copies -- under-enriches versus a real binder of the
# same apparent SC50, so a very negative score flags it. The apparent SC50 used for the cohort
# is our own sc50_est (non-binders / gated designs have inf SC50 and are excluded). This is a
# SOFT score (not a reliable hard gate); see progress_passenger.md.
n_designs = len(counts)
fit_bcs = list(sorts[fit_sort_mask].index)
fit_concs = sorts.loc[fit_bcs, 'concentration'].values.astype(float)
hi_conc_bc = fit_bcs[int(np.argmax(fit_concs))]
expr_bc = sorts.loc[hi_conc_bc, 'expression_parent']

log_enrich_titr = np.full(n_designs, np.nan)
if expr_bc != "" and expr_bc in counts.columns:
    norm_expr = (counts[expr_bc] / counts[expr_bc].sum()).values
    norm_hi = (counts[hi_conc_bc] / counts[hi_conc_bc].sum()).values
    with np.errstate(divide="ignore", invalid="ignore"):
        log_enrich_titr = np.log(norm_hi) - np.log(norm_expr)
    log_enrich_titr[~np.isfinite(log_enrich_titr)] = np.nan
else:
    print("Warning! Could not find an expression pool for the top titration sort %s; "
          "rel_enrich_titr will be NaN." % hi_conc_bc)

# Cohort by the apparent (uncensored) SC50 so the same-affinity reference is built over the
# full population, not just the in-range binders that survive the gate.
log_sc50 = np.log(sc50_apparent)
rel_enrich_titr = cohort_relative_enrich(log_enrich_titr, log_sc50, valid=np.isfinite(log_sc50))

# Threshold in fractional space: exp(rel_enrich_titr) < threshold means this design enriches
# at less than (threshold * cohort reference). Default 1/20 = 0.05.
suspected_passenger = np.isfinite(rel_enrich_titr) & (np.exp(rel_enrich_titr) < args.suspected_passenger_threshold)


# ---------------------------------------------------------------------------------
# Assemble output.
# ---------------------------------------------------------------------------------
records = {
    "sc50_est": sc50_est,
    "sc50_lb": sc50_lb,
    "sc50_ub": sc50_ub,
    "sc50_unmeasurable": sc50_unmeasurable,
    "sc50_unmeasurable_cutoff": np.full(n_designs, gate_cutoff),
    "p_bind": p_bind,
    "p_bind_threshold": np.full(n_designs, p_bind_threshold),
    "rel_enrich_titr": np.exp(rel_enrich_titr),
    "suspected_passenger": suspected_passenger,
    "low_conf": low_conf,
    "lowest_conc": np.full(n_designs, lowest_concentration),
    "highest_conc": np.full(n_designs, highest_concentration),
}
if ( args.target != "" ):
    records['target'] = args.target
records['description'] = counts['description']

df = pd.DataFrame(records)
df = df.merge(counts, 'inner', 'description')


# Column order: analysis columns first, description last. Drop merged per-pool count columns
# unless --dump_round_information was requested.
new_key_order = list(df)
if not args.dump_round_information:
    for key in list(counts):
        if key in new_key_order:
            if key == 'description':
                continue
            new_key_order.remove(key)
new_key_order.remove('description')
new_key_order.append('description')

df[new_key_order].to_csv(args.output_file, sep=" ", index=None, na_rep="NaN")
