#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import sys
import json
import itertools
import re
import random
import numpy as np
import math
import glob
from collections import defaultdict
from scipy.optimize import curve_fit
import argparse
import torch
print("Done importing")


# In[39]:


# args = os


default_contig = "A1-6,A7-15,0 B1-118"
default_pdb_name = "third_go_-_fences_7_v1_2_2_4_2_4_2_5_n402__1_t005_0001.pdb"
default_int_out = "xy_str=[[3.929,4.439],[-3.929,-4.439]];aligned_w_z=[True,False];motif_sizes=[6,9];motif_can_add=[[1,19],[3,19]]"

default_contig = "A1-11,0 B1-146"
default_pdb_name = "asdf"
default_int_out = "xy_str=[[0.000,0.000]];aligned_w_z=[True];motif_sizes=[11];motif_can_add=[[2,33]]"

default_contig = "A1-12,A13-20,0 B1-153"
default_pdb_name = "asdf"
default_int_out = "xy_str=[[-5.014,-1.239],[5.014,1.239]];aligned_w_z=[True,False];motif_sizes=[12,8];motif_can_add=[[33,33],[33,6]]"

# ferredoxin
default_contig = 'A1-6,A7-19,A20-24,A25-29,A30-42,A43-48,0 B1-153'
default_pdb_name = 'asdf'
default_int_out = 'xy_str=[[-3.590,-1.611],[7.702,3.350],[-1.245,7.564],[-2.380,2.956],[4.874,-6.298],[-5.361,-5.961]];aligned_w_z=[True,False,True,False,True,False];motif_sizes=[6,13,5,5,13,6];motif_can_add=[[16,16],[33,33],[16,16],[16,16],[2,0],[16,0]];strand_isheet=[0,np.nan,0,0,np.nan,0];isheet_units=[[0.288,0.958]];sheet_can_extend=[[False, True]]'

# test sheet
default_contig = 'A1-5,A6-10,0 B1-1'
default_int_out = 'xy_str=[[5.117,-0.168],[-5.117,0.168]];aligned_w_z=[True,True];motif_sizes=[5,5];motif_can_add=[[33,33],[33,33]];strand_isheet=[np.nan,0];isheet_units=[[0.342,0.940]];sheet_can_extend=[[False,True]]'

# test 2sheet
default_contig = 'A1-5,A6-10,0_B1-1' 
default_int_out= 'xy_str=[[5.117,-0.168],[-5.117,0.168]];aligned_w_z=[True,True];motif_sizes=[5,5];motif_can_add=[[33,33],[33,33]];strand_isheet=[np.nan,0];isheet_units=[[0.342,0.940]];sheet_can_extend=[[False,True]]'

#
default_contig = 'A1-4,A5-8,0 B1-198'
default_int_out= 'xy_str=[[-4.194,-0.544],[4.194,0.544]];aligned_w_z=[True,False];motif_sizes=[4,4];motif_can_add=[[16,1],[16,3]];strand_isheet=[0,0];isheet_units=[[-0.911,-0.411]];sheet_can_extend=[[True,False]]'

parser = argparse.ArgumentParser()
parser.add_argument("-contig", type=str, default=default_contig)
parser.add_argument("-pdb_name", type=str, default=default_pdb_name)
parser.add_argument("-interface_output", type=str, default=default_int_out)
parser.add_argument("-output_dir", type=str, default="./")
parser.add_argument("-num_to_output", type=int, default=100)
parser.add_argument("-flank_angle", type=float, default=45)
parser.add_argument("-max_delta_helix_in_protein", type=int, default=12)
parser.add_argument("-max_delta_touching_helix", type=int, default=8)
parser.add_argument("-loop_provision", type=int, default=4)
parser.add_argument("-min_helix_size", type=int, default=12)
parser.add_argument("-min_protein_size", type=int, default=0)
parser.add_argument("-max_protein_size", type=int, default=120)
parser.add_argument("-max_helix_size", type=int, default=35)
parser.add_argument("-max_considered_helices", type=int, default=9)
parser.add_argument("-min_considered_helices", type=int, default=3)
parser.add_argument("-max_helices_given_n_input_helices", type=str, default="1:4,2:6,3:7,4:9,5:12")
parser.add_argument("-max_protein_size_helix_per_res", type=int, default=35)
parser.add_argument("-helix_len_step", type=int, default=2)
parser.add_argument("-dont_balance_lengths", action='store_true')
parser.add_argument("-min_sheet_width", type=int, default=2)
parser.add_argument("-max_sheet_width", type=int, default=4)
parser.add_argument("-min_strand_size", type=int, default=3)
parser.add_argument("-max_strand_size", type=int, default=15)
parser.add_argument("-no_motif_actually", action="store_true")
parser.add_argument("-no_helices", action="store_true")

is_cmd = False
if ("ipykern" in sys.argv[0] ):
    args = parser.parse_args([])
else:   
    args = parser.parse_args(sys.argv[1:])
    is_cmd = True

fname_prefix = os.path.join(args.output_dir, os.path.basename(args.pdb_name).replace(".pdb", "") + "__")


# In[40]:



max_helices_given_n_input_helices = {}
if ( args.max_helices_given_n_input_helices != "" ):
    max_helices_given_n_input_helices = eval('{' + args.max_helices_given_n_input_helices + '}')


# In[41]:


if ( is_cmd ):

    def do_nothing(asdf=None, b=None, c=None, **kwargs):
        pass
    def return_something(asdf=None, b=None, c=None, **kwargs):
        return []

    def return_plt():
        return plt

    plt = sys
    plt.show = do_nothing
    plt.figure = do_nothing
    plt.xlim = do_nothing
    plt.ylim = do_nothing
    plt.scatter = do_nothing
    plt.gca = return_plt
    plt.add_patch = do_nothing
    plt.plot = do_nothing
    matplotlib = plt
    matplotlib.patches = plt
    plt.Circle = return_something
else:
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    


# In[42]:


# xy_str=[[-8.935,-0.590],[0.192,3.779],[8.743,-3.189]];aligned_w_z=[True,False,False];motif_sizes=[9,8,12];motif_can_add=[[2,19],[19,19],[19,3]]
# empty_contig = "A1-9,A10-17,A18-29,0 B1-118"

# xy_str=[[8.162,2.968],[-0.033,3.497],[-8.129,-6.465]];aligned_w_z=[True,False,True];motif_sizes=[9,9,8];motif_can_add=[[6,19],[2,19],[19,1]]
# empty_contig = "A1-9,A10-18,A19-26,0 B1-118"

strand_isheet = []
isheet_units = []
sheet_can_extend = []

empty_contig = args.contig
exec(args.interface_output)
assert("xy_str" in locals())
assert("aligned_w_z" in locals())
assert("motif_sizes" in locals())
assert("motif_can_add" in locals())

xy_str = np.array(xy_str)
# z_bounds = np.array(z_bounds)
motif_sizes = np.array(motif_sizes)
motif_can_add = np.array(motif_can_add)
strand_isheet = np.array(strand_isheet)
isheet_units = np.array(isheet_units)

assert(len(xy_str) == len(motif_sizes))
assert(len(xy_str) == len(aligned_w_z))
assert(len(xy_str) == len(motif_can_add))

if len(strand_isheet) == 0:
    strand_isheet = np.full(len(xy_str), np.nan)
else:
    max_helices_given_n_input_helices = {}
    
assert(len(xy_str) == len(strand_isheet))

use_flank_positive_x = args.flank_angle
use_flank_negative_x = args.flank_angle

if ( "flank_positive_x" in args.interface_output ):
    use_flank_positive_x = flank_positive_x
    
if ( "flank_negative_x" in args.interface_output ):
    use_flank_negative_x = flank_negative_x


# In[43]:


def get_sheet_info(xy_str, aligned_w_z, strand_isheet, isheet_units):
    coms = []
    sizes = []
    orderings = []
    for isheet in range(len(isheet_units)):
        unit = isheet_units[isheet]
        mask = strand_isheet == isheet
        sheet_coords = xy_str[mask]
        com = np.mean( sheet_coords , axis=0 )
        coords_from_com = sheet_coords - com
        projected_dist = np.sum( coords_from_com * unit, axis=-1)
        size = np.max( np.abs(projected_dist))
        local_order = np.argsort(projected_dist)
        order = np.where(mask)[0][local_order]
        
        coms.append(com)
        sizes.append(size)
        orderings.append(order)
        
    return coms, sizes, orderings
        
        
        

def plot_circles(xy_str, aligned_w_z, strand_isheet, isheet_units, scale=1.2, looping=None ):
    colors = [
        'blue',
        'red',
    ]


    plt.figure(figsize=(2, 1.5))
    plt.xlim(-20*scale, 20*scale)
    plt.ylim(-15*scale, 15*scale)

    for i in range(len(xy_str)):
        xy = xy_str[i]

        aligned = aligned_w_z[i]

        circle = matplotlib.patches.Circle(xy, radius=4, color=colors[int(aligned)])

        plt.scatter([xy[0]], [xy[1]])
        plt.gca().add_patch(circle)
        
#         assert((circle.center == xy).all())
        
        
    sheet_coms, sheet_sizes, sheet_orderings = get_sheet_info(xy_str, aligned_w_z, strand_isheet, isheet_units)
    for isheet in range(len(isheet_units)):
        start = sheet_coms[isheet] - (sheet_sizes[isheet] + 8) * isheet_units[isheet]
        end = sheet_coms[isheet] + (sheet_sizes[isheet] + 8 ) * isheet_units[isheet]
        
        plt.plot([start[0], end[0]], [start[1], end[1]], color='#00FF00', linewidth='4')
        
        
    if (looping is not None):
        xs = xy_str[looping,0]
        ys = xy_str[looping,1]
        
        plt.plot(xs, ys, color='black')
        
        
plot_circles(xy_str, aligned_w_z, strand_isheet, isheet_units)

    
plt.show()


# In[44]:


def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def lines_intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


# In[45]:




# Function for projecting some vector a onto b
def proj(a, b):
    k = np.sum(a * b) / np.sum(b * b) # use squared norm because we divide by it twice
    return k * b


# Returns the distance from line segment AB to point C
def distanceSegmentToPoint2(A, B, C):
    # Compute vectors AC and AB
    AC = C - A
    AB = B - A

    # Get point D by taking the projection of AC onto AB then adding the offset of A
    D = proj(AC, AB) + A

    AD = D - A
    # D might not be on AB so calculate k of D down AB (aka solve AD = k * AB)
    # We can use either component, but choose larger value to reduce the chance of dividing by zero
    k = AD[0] / AB[0] if np.abs(AB[0]) > np.abs(AB[1]) else AD[1] / AB[1];

    # Check if D is off either end of the line segment
    if (k <= 0.0):
        return np.sum(np.square(C -A))
    elif (k >= 1.0):
        return np.sum(np.square(C - B))
    else:
        return np.sum(np.square(C - D))

    
def line_circle_intersect(pt1, pt2, center, radius):
    distance = distanceSegmentToPoint2(pt1, pt2, center)
    return distance < radius**2

def circle_circle_intersect(center1, center2, radius):
    return np.linalg.norm(center1 - center2) < radius*2

def check_circle_vs_lines(center, radius, lines):
    any_intersection = False
    for start, end in lines:
        intersect = line_circle_intersect(start, end, center, radius)
        if ( intersect ):
            any_intersection = True
            break
    return any_intersection

def check_circle_vs_circles(center, radius, centers):
    any_intersection = False
    for other_center in centers:
        intersect = circle_circle_intersect(center, other_center, radius)
        if ( intersect ):
            any_intersection = True
            break
    return any_intersection

def whose_clashing_circles(i, radius, centers):
    whose_clashing_circle = np.zeros(len(centers), bool)
    for other_i in range(len(centers)):
        if ( other_i == i ):
            continue
        intersect = circle_circle_intersect(centers[i], centers[other_i], radius)
        whose_clashing_circle[other_i] = intersect
    return whose_clashing_circle


# In[46]:


def redraw_circles(circles, circle_base, color):
    for i in range(len(circles)):
        xy = circles[i]
        base = circle_base[i]
        if ( base < 0 ):
            continue
        
        circle = matplotlib.patches.Circle(xy, radius=4, color=color)
        plt.gca().add_patch(circle)


# In[47]:



def bump_circle(circles, circle_base, i, amount):
    circles = circles.copy()
    base_xy = circles[circle_base[i]]
    
    center = circles[i]
    
    delta = center - base_xy
    
    angle = np.degrees(np.arctan2(delta[1], delta[0]))
    dist = np.linalg.norm(delta)
    
    new_angle = np.radians(angle + amount)
    
    new_center = base_xy + np.array([dist * np.cos(new_angle), dist * np.sin(new_angle)])

    circles[i] = new_center
    
    return circles

        
def generate_connecting_lines(xy_str, close_dist=1100):

    dists = np.linalg.norm( xy_str[:,None] - xy_str[None,:], axis=-1 )

    close_enough = (dists < close_dist) & (dists > 0)

    lines = []

    for i, j in zip(*np.where(close_enough)):
        if ( i < j ):
            continue
        xs = xy_str[[i, j], 0]
        ys = xy_str[[i, j], 1]

        lines.append([xy_str[i], xy_str[j]])
    return lines


def find_lower_edge(xy_str, lines):
    vertical_lines = []
    for i in range(len(xy_str)):

        pt = xy_str[i]

        start_pt = pt.copy()
        start_pt[1] -= 1
        end_pt = pt.copy()
        end_pt[1] = -30

        any_intersection = False
        for iline in range(len(lines)):

            seg_start, seg_end = lines[iline]

            intersect = lines_intersect(start_pt, end_pt, seg_start, seg_end)
            if (intersect ):
                any_intersection = True
                break

        if ( any_intersection ):
            continue
        vertical_lines.append(i)

#         plt.plot([start_pt[0], end_pt[0]], [start_pt[1], end_pt[1]], color='grey')

    vertical_lines = np.array(vertical_lines)
    
    return vertical_lines

def generate_flanks(xy_str, vertical_lines, negative_flank_angle=45, positive_flank_angle=45):
    
    lines = []
    x_sorted = list(sorted(vertical_lines, key=lambda x: xy_str[x][0]))

#     angle = 20
    left_frac = np.tan(np.radians(negative_flank_angle))
    left_edge = xy_str[x_sorted[0]]
    left_pt = left_edge.copy()
    left_pt += np.array([-10, -10*left_frac])

    lines.append([left_edge, left_pt])


    right_frac = np.tan(np.radians(positive_flank_angle))
    right_edge = xy_str[x_sorted[-1]]
    right_pt = right_edge.copy()
    right_pt += np.array([10, -10*right_frac])

    lines.append([right_edge, right_pt])

    return lines
def generate_lower_guard_lines(lines):

    new_lines = []
    for start_pt, end_pt in lines:
        xs = np.linspace(start_pt[0], end_pt[0], 5)
        ys = np.linspace(start_pt[1], end_pt[1], 5)

        for x, y in zip(list(xs), list(ys)):
            pt = np.array([x, y])
            lower = pt.copy()
            lower[1] = -30
            new_lines.append([pt, lower])
            
    return new_lines
        
    
def plot_connecting_lines(lines):
    for pt1, pt2 in lines:
        xs = [pt1[0], pt2[0]]
        ys = [pt1[1], pt2[1]]
        plt.plot(xs, ys, color='black')
        

def generate_starting_circles(xy_str, lines, is_new, theta_step=-5):

    x_sorted = list(sorted(np.arange(len(xy_str)), key=lambda x: xy_str[x][0]))

    circles = list(xy_str)
    circle_base = []
    for i in range(len(circles)):
        if is_new[i]:
            circle_base.append(-2)
        else:
            circle_base.append(-1)
#     theta_step = -10

    for i_x, i_basis in enumerate(x_sorted):

        cur_theta = 270

        basis_xy = xy_str[i_basis]

        while (cur_theta > -90):
            cur_theta += theta_step
            trial_x = basis_xy[0] + 8*np.cos(np.radians(cur_theta))
            trial_y = basis_xy[1] + 8*np.sin(np.radians(cur_theta))

            trial_pt = np.array([trial_x, trial_y])

            intersect_line = check_circle_vs_lines(trial_pt, 3.99, np.array(lines))
            if ( intersect_line ):
                continue
            intersect_circle = check_circle_vs_circles(trial_pt, 3.99, circles)
            if ( intersect_circle ):
                continue

#             circle = matplotlib.patches.Circle(trial_pt, radius=4, color='green')
#             plt.gca().add_patch(circle)
            circles.append(trial_pt)
            circle_base.append(i_basis)

    circles = np.array(circles)
    circle_base = np.array(circle_base)
    
    return circles, circle_base


    
def shove_towards_the_center(circles, circle_base, start_left=True):
    
    shove_mask = circle_base >= 0
    circles = circles.copy()
    movable = circle_base >= 0
    movable_idx = np.where(movable)[0]
    
    if ( len(movable_idx) == 0 ):
        return circles
    move_order = []
    for i in range(len(movable_idx)//2):
        
        if ( start_left ):
            move_order.append(-movable_idx[i])
            move_order.append(movable_idx[-i-1])
        else:
            move_order.append(movable_idx[-i-1])
            move_order.append(-movable_idx[i])
    if ( len(movable_idx) % 2 == 1):
        move_order.append((-1 if start_left else 1) * movable_idx[len(movable_idx)//2])
    
    assert(len(move_order) == len(movable_idx))
    assert(len(np.unique(np.abs(move_order))) == len(movable_idx))
    
    for imove in move_order:
        
        move_amount = np.sign(imove) * 5
        imove = abs(imove)
        while True:
            circles, result = shove_circle_recursive(circles, circle_base, imove, move_amount, shove_mask)

            if ( not result ):
                break
    
        shove_mask[imove] = False
    return circles


# In[48]:



# returns true on success
def shove_circle_recursive(circles, circle_base, i, amount, shove_mask, trial_bump=5 ):
#     print("Shoving %i by %i"%(i, amount))
    is_template = (circle_base < 0) | ~shove_mask
    
    test_circles = bump_circle(circles, circle_base, i, amount)

    intersect_line = check_circle_vs_lines(test_circles[i], 3.99, np.array(lines))
    if ( intersect_line ):
        return circles, False
    
    
    while True:
        whose_clashing = whose_clashing_circles(i, 3.99, test_circles)
#         print(whose_clashing)
        
        if ( not whose_clashing.any() ):
            return test_circles, True

        clash_is_template = (whose_clashing & is_template).any()
        if ( clash_is_template ):
            return circles, False

        expected_clash = i + 1 if amount < 0 else i - 1
        unexpected_clash = i - 1 if amount < 0 else i + 1
        
        # ok, there's a clash that's unexpected, it's probably a backwards clash
        # with a circle above us
        if ( expected_clash == len(circles) or is_template[expected_clash] or not whose_clashing[expected_clash] ):
            try:
                if not (whose_clashing[unexpected_clash]):
                    print("Truly unexpected clash")
            except:
                print("Truly unexpected clash")
            return test_circles, True
    
        test_circles, success = shove_circle_recursive(test_circles, circle_base, expected_clash,
                                                                   np.sign(amount)*trial_bump, shove_mask)
        
        if ( not success ):
            return circles, False
        


def hamiltonian_recursive(edge_mask, path):
    if ( len(path) == len(edge_mask) ):
        return [path]
    our_i = path[-1]
    
    potential_neighbors = np.where( edge_mask[our_i]  )[0]
    
    all_paths = []
    for neighbor in potential_neighbors:
        if ( neighbor in path ):
            continue
        all_paths += hamiltonian_recursive(edge_mask, path + [neighbor])
    
    return all_paths
    

def find_all_loopings(edge_mask):
    all_paths = []
    for i in range(len(edge_mask)):
        all_paths += hamiltonian_recursive(edge_mask, [i])
    return all_paths


# In[49]:



if ( len(xy_str) in max_helices_given_n_input_helices):
    args.max_considered_helices = max_helices_given_n_input_helices[len(xy_str)]





# how many Å wide is a strand?
strand_width = 4.7

# extend sheets
# 

# digit 0 is just a placeholder so there can be no sheets
# digit 0 is always 0
# digit 1/2 are extensions of the negative and positive sides of the first sheet
# digit 3/4 are ...

digits = 1 + len(sheet_can_extend)*2
digit_base = args.max_sheet_width

sheet_scenario_ub = digit_base ** digits

ideas_by_size_by_scenario = {}

for sheet_scenario in range(sheet_scenario_ub):
    if sheet_scenario % digit_base != 0:
        continue
    
    # figure out how many to add to each sheet based on scenario number or reject
    acceptable = True
    sheet_additions = []
    for isheet in range(len(sheet_can_extend)):
        negative_digit = isheet*2 + 1
        positive_digit = isheet*2 + 2
        
        sheet_cur_size = (strand_isheet == isheet).sum()
        
        add_negative = ( sheet_scenario // (digit_base ** negative_digit ) ) % digit_base
        add_positive = ( sheet_scenario // (digit_base ** positive_digit ) ) % digit_base
        
        if add_negative > 0 and not sheet_can_extend[isheet][0]:
            acceptable = False
            break
        if add_positive > 0 and not sheet_can_extend[isheet][1]:
            acceptable = False
            break
        
        new_size = add_negative + sheet_cur_size + add_positive
        if new_size < args.min_sheet_width:
            acceptable = False
            break
        if new_size > args.max_sheet_width:
            acceptable = False
            break
            
        sheet_additions.append([add_negative, add_positive])

    if not acceptable:
        continue
    
    print(sheet_additions)
    
    
    # actually perform the additions to the sheets
    sheet_coms, sheet_sizes, sheet_orderings = get_sheet_info(xy_str, aligned_w_z, strand_isheet, isheet_units)
#     print(len(aligned_w_z))
    
    local_xy_str = list(xy_str)
    local_strand_isheet = list(strand_isheet)
    local_aligned_w_z = list(aligned_w_z)
    is_new = [False]*len(xy_str)
    
    for isheet in range(len(isheet_units)):
        add_negative, add_positive = sheet_additions[isheet]
        unit = isheet_units[isheet]
        
        negative_term = sheet_orderings[isheet][0]
        positive_term = sheet_orderings[isheet][-1]
        
        negative_aligned = aligned_w_z[negative_term]
        positive_aligned = aligned_w_z[positive_term]
        
        negative_xy = xy_str[negative_term]
        positive_xy = xy_str[positive_term]
        
        for iadd in range(1, add_negative+1):
            xy = negative_xy + unit * iadd * strand_width * -1
            with_z = negative_aligned
            for i in range(iadd):
                with_z = not with_z
            local_xy_str.append(xy)
            local_strand_isheet.append(isheet)
            local_aligned_w_z.append(with_z)
            is_new.append(True)
            
        for iadd in range(1, add_positive+1):
            xy = positive_xy + unit * iadd * strand_width
            with_z = positive_aligned
            for i in range(iadd):
                with_z = not with_z
            local_xy_str.append(xy)
            local_strand_isheet.append(isheet)
            local_aligned_w_z.append(with_z)
            is_new.append(True)
    
    local_xy_str = np.array(local_xy_str)
    local_strand_isheet = np.array(local_strand_isheet)
    local_aligned_w_z = np.array(local_aligned_w_z)
#     print(len(local_aligned_w_z))
    
            
    lines = generate_connecting_lines(local_xy_str)

    vertical_lines = find_lower_edge(local_xy_str, lines)


    lines += generate_flanks(local_xy_str, vertical_lines, 
                             negative_flank_angle=use_flank_negative_x, 
                             positive_flank_angle=use_flank_positive_x)

    lines += generate_lower_guard_lines(lines)



    starting_circles, starting_circle_base = generate_starting_circles(local_xy_str, lines, is_new)


    num_test_circles = (starting_circle_base >= 0).sum()
    num_templ = (starting_circle_base < 0).sum()

    circle_ideas = []
    circle_idea_bases = []
    
    
    motif_accounted_size = ( np.isnan(local_strand_isheet).sum() * (args.min_helix_size + args.loop_provision)
                             + (~np.isnan(local_strand_isheet)).sum() * (args.min_strand_size + args.loop_provision)
                           )
    max_realistic_helices = (args.max_protein_size - motif_accounted_size) // args.min_helix_size + 1 + num_templ


    def circle_patterns_are_equal(circles1, circles2, tol=1):
        if ( len(circles1) != len(circles2 )):
            return False
        dists2 = np.sum(np.square(circles1-circles2), axis=-1)
        return not np.any( dists2 > tol**2 )

    def circle_pattern_exists(new_circles, all_circles, tol=1):
        for other_circles in all_circles:
    #         print("here")
            if ( circle_patterns_are_equal(new_circles, other_circles, tol=tol)):
                return True
        return False

    for i_idea in range(2**num_test_circles):

        fmt = "0%ib"%num_test_circles
        binary = format(i_idea, fmt).replace("b", "")
        mask = np.array(list(binary)) == "1"

        use_circle_mask = np.zeros(len(starting_circles), bool)
        use_circle_mask[:num_templ] = True
        use_circle_mask[num_templ:] = mask
        
        
        if args.no_helices:
            if mask.sum() > 0:
                continue

#         print(use_circle_mask.sum(), args.min_considered_helices, max_realistic_helices )
        if ( use_circle_mask.sum() > max_realistic_helices): #args.max_considered_helices ):
            continue
        if ( use_circle_mask.sum() < args.min_considered_helices ):
            continue
            
#         print("here?")

        for start_left in [True, False]:
            circles = starting_circles[use_circle_mask]
            circle_base = starting_circle_base[use_circle_mask]

            circles = shove_towards_the_center(circles, circle_base, start_left)
            if ( not circle_pattern_exists( circles, circle_ideas ) ):
                circle_ideas.append(circles)
                circle_idea_bases.append(circle_base)
    #         break


    #     if ( i_idea > 100 ):
    #         break
    # circles = bump_circle(circles, circle_base, 4, -10)

    print("Ideas:", len(circle_ideas))
    

    ideas_by_size = defaultdict(list)
    for idea, bases in zip(circle_ideas, circle_idea_bases):
        size = len(idea)
        ideas_by_size[size].append((idea, bases, local_xy_str, local_strand_isheet, local_aligned_w_z))




    for circles, circle_base in zip(circle_ideas, circle_idea_bases):

        plot_circles(local_xy_str, local_aligned_w_z, local_strand_isheet, isheet_units)
        plot_connecting_lines(lines)
        redraw_circles(starting_circles, starting_circle_base, 'green')
        redraw_circles(circles, circle_base, 'cyan')
        plt.show()
        
        print("debug on printing")
        break
        
        
    ideas_by_size_by_scenario[sheet_scenario] = ideas_by_size

#     break


# In[51]:


max_realistic_helices


# In[52]:


final_ideas = []
final_bases = []
final_locals = []

for scenario in ideas_by_size_by_scenario:
    ideas_by_size = ideas_by_size_by_scenario[scenario]

    for size in ideas_by_size:
        scores = []
        ideas = []
        bases = []
        for idea, base, local_xy_str, local_strand_isheet, local_aligned_w_z in ideas_by_size[size]:
            ideas.append(idea)
            bases.append(base)

            localer_strand_isheet = np.full(len(idea), np.nan)
            localer_strand_isheet[:len(local_strand_isheet)] = local_strand_isheet
            
            is_templ = (base <= -1) & np.isnan(localer_strand_isheet)
            is_strand = (base <= -1) & ~np.isnan(localer_strand_isheet)

            templs = idea[is_templ]
            non_templs = idea[~is_templ & ~is_strand]
            strands = idea[is_strand]
            strand_belongs_to = localer_strand_isheet[is_strand]

            if ( len(non_templs) == 0 ):
                scores.append(0)
                continue

            touching_templ = (np.linalg.norm( non_templs[:,None] - templs[None,:], axis=-1 ) < 9).sum()
            non_temple_touch = (np.linalg.norm( non_templs[:,None] - non_templs[None,:], axis=-1) < 9).sum()/2
            
            # you only get one point for touching a sheet
            strand_touch_matrix = ( np.linalg.norm( non_templs[:,None] - strands[None,:], axis=-1) < 9)
            sheet_touches = 0
            for isheet in range(len(isheet_units)):
                mask = strand_belongs_to == isheet
                sheet_touches += strand_touch_matrix[:,mask].any(axis=-1).sum(axis=0)
                

            score = touching_templ + sheet_touches + non_temple_touch * 0.1

            scores.append(score)

        scores = np.array(scores)
        ideas = np.array(ideas)
        bases = np.array(bases)

        best_to_worst = np.argsort(-scores)
        best_idx = best_to_worst[0]

        best_score = scores[best_idx]
        if len(isheet_units) > 0:
            best_score -= 1

        print(size)
        for i, idx in enumerate(np.where(scores >= best_score)[0]):
    #         print(idx)
            circles = ideas[idx]
            circle_base = bases[idx]

            if i == 0:
                plot_circles(local_xy_str, local_aligned_w_z, local_strand_isheet, isheet_units)
                redraw_circles(circles, circle_base, 'green')
                plt.show()
                print("debug on priting")

            final_ideas.append(circles)
            final_bases.append(circle_base)
            final_locals.append((local_xy_str, local_strand_isheet, local_aligned_w_z))

    print("Final ideas:", len(final_ideas))


# In[53]:



#https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
def line_segment_intersection(start1, end1, start2, end2):
    p = start1
    r = end1 - start1
    q = start2
    s = end2 - start2
    
    # collinear or parallel 
    if np.cross(r, s) == 0:
        # parallel
        if np.cross(q - p, r) != 0:
            return False
        t0 = np.sum( (q - p) * r ) / np.sum( r * r )
        t1 = t0 + np.sum( s * r ) / np.sum( r * r )
        if t1 < t0:
            t0,t1 = t1,t0
        # check if t0,t1 intersects 0,1
        if t0 > 1:
            return False
        if t1 < 0:
            return False
        return True
    
    t0 = np.cross(q-p, s) / np.cross(r, s)
    t1 = np.cross(q-p, r) / np.cross(r, s)
    
    if np.abs(t0-0.5) > 0.5 or np.abs(t1-0.5) > 0.5:
        return False
    
    return True


# In[ ]:





# In[58]:



final_looped_dreams = []

super_break = False
for idea, base, (local_xy_str, local_strand_isheet, local_aligned_w_z) in zip(
                                                                    final_ideas, final_bases, final_locals):

#     print("Debug")
#     if not (base == -2).any():
#         continue
# idea = final_ideas[9]
# base = final_bases[9]

    localer_strand_isheet = np.full(len(idea), np.nan)
    localer_strand_isheet[:len(local_strand_isheet)] = local_strand_isheet

    loop_dist = 12
    if len(isheet_units) > 0:
        loop_dist = 13 # the ferredoxin needed 13

    all_by_dist = np.linalg.norm(idea[:,None] - idea[None,:], axis=-1)
    can_loop_grid = (all_by_dist < loop_dist) & (all_by_dist > 0)

    num_test_circles = (base >= 0).sum()
    num_templ = (base < 0).sum()
    is_templ = (base < 0) & (np.isnan(localer_strand_isheet))
    is_strand = (base < 0) & (~np.isnan(localer_strand_isheet))
    
    same_sheet = localer_strand_isheet[:,None] == localer_strand_isheet[None,:]
    assert(np.nan != np.nan)

    looped_dreams = []

    total = 0
    for idream in range(2**num_test_circles):
    #     idream = 7

        fmt = "0%ib"%num_test_circles
        binary = format(idream, fmt).replace("b", "")
        mask = np.array(list(binary)) == "1"

        dream_with_z = np.zeros(len(idea), bool)
        dream_with_z[:num_templ] = local_aligned_w_z
        dream_with_z[num_templ:] = mask

        with_z_can_loop = dream_with_z[:,None] ^ dream_with_z[None,:]

        dream_can_loop = can_loop_grid & with_z_can_loop

        all_loopings = find_all_loopings(dream_can_loop)
        if ( len(all_loopings) == 0 ):
            continue
        all_loopings = np.array(all_loopings)

        total += len(all_loopings)

        for looping in all_loopings:
            looped_dreams.append((dream_with_z, looping))

#     break
    
#     print(total)
    
    cos = []
    for with_z, looping in looped_dreams:
        list_looping = list(looping)
        order = []
        for i in range(len(looping)):
            order.append(list_looping.index(i))
        order = np.array(order)
        assert((order >= 0).all())
        
        reordered_idea = idea[looping]
        reordered_can_loop = can_loop_grid[looping][:,looping]
        reordered_is_templ = is_templ[looping]
        reordered_is_strand = is_strand[looping]
        reordered_same_sheet = same_sheet[looping][:,looping]
        
        reordered_is_fixed = reordered_is_templ | reordered_is_strand

        # bundle is part of a triangle
        #  considered part of a triangle if its contacts can bridge each other
        # we don't care what its contacts are
        #  anything that bridges two fixeds counts
        is_bridged = np.zeros(len(reordered_can_loop), int)
        for i in range(1, len(order)-1):
            is_bridge = reordered_can_loop[i-1, i+1]
            if ( reordered_is_fixed[i-1] and reordered_is_fixed[i+1]  ):
                is_bridge = True
            if ( is_bridge ):
                is_bridged[i-1:i+1+1] = True
        
        # figure out how far away nearest primary-adjacent contact is for non-bridged helices
        # if it's part of a bridge, then 0
        # if it's touching one-past a neighbor that's 2 (two past is 3)
        #
        # if it's not touching anyone besides the primary-adjecent, it's unsupported
        #
        # code below tosses out anything that's not bridged
        #
        # min_delta_is is a histogram
        min_delta_is = np.zeros(len(reordered_can_loop), int)
        unsupported_helix = False
        for i_elem in range(len(order)):
            if ( is_bridged[i_elem] ):
                min_delta_is[0] += 1
                continue
            contacts = np.where(reordered_can_loop[i_elem])[0]
            delta = np.abs(i_elem - contacts)
            if ( i_elem == 0 or i_elem == len(reordered_can_loop)-1):
                delta = np.r_[delta, np.array([1])]
                
            if ( len(delta) == 2 ):
                unsupported_helix = True
                continue
            
            best_idx = np.argsort(delta)[2]
            best_delta = delta[best_idx]
            min_delta_is[best_delta] += 1

            
        # strand rules
        # you're only allowed to form hairpins, no triple strands
        strand_rule_broken = False
        for i_elem in range(1, len(order)-1):
            if reordered_same_sheet[i_elem, i_elem-1] and reordered_same_sheet[i_elem, i_elem+1]:
                strand_rule_broken = True
                break
            
            
        # no crossing loops
        even_istarts = np.arange(0, len(order)-1, 2)
        odd_istarts = np.arange(1, len(order)-1, 2)
        
        even_iends = even_istarts + 1
        odd_iends = odd_istarts + 1
        
        even_starts = reordered_idea[even_istarts]
        odd_starts = reordered_idea[odd_istarts]
        even_ends = reordered_idea[even_iends]
        odd_ends = reordered_idea[odd_iends]
        
        crossing = False
        for i in range(len(even_starts)-1):
            for j in range(i+1, len(even_starts)):
                if line_segment_intersection( even_starts[i], even_ends[i], even_starts[j], even_ends[j] ):
                    crossing = True
                    break
        
        for i in range(len(odd_starts)-1):
            for j in range(i+1, len(odd_starts)):
                if line_segment_intersection( odd_starts[i], odd_ends[i], odd_starts[j], odd_ends[j] ):
                    crossing = True
                    break
                    
        
        
            
#         try:
#             if ( np.allclose(looping, np.array([5, 4, 3, 0, 1, 2])) ):
#                 print(is_bridged)
#                 print(min_delta_is)
#                 plot_circles(idea, with_z, localer_strand_isheet, isheet_units, looping=looping)
#                 plt.show()
#                 super_break = True
#         except:
#             pass


        if min_delta_is.sum() == 0:
            co1 = 1000
        else:
            co1 = np.where(min_delta_is)[0][-1]
            
        if ( super_break ):
            break
        # they get duplicated but we're taking a mean anyways
        contacts1, contacts2 = np.where(reordered_can_loop)

        order_delta = np.abs(contacts1 - contacts2)
#         order_delta[order_delta <= 2] = 0

        co = np.mean(order_delta)

        if ( unsupported_helix ):
            co = np.nan
        if ( co1 > 0 ):
            co = np.nan
        if strand_rule_broken:
            co = np.nan
        if crossing:
            co = np.nan
        cos.append(co)
        
    if ( super_break ):
        break
    if (len(looped_dreams) == 0):
        continue
    if ( np.isnan(cos).all() ):
        continue
    cos = np.array(cos)
    lowest = np.nanmin(cos)
#     lowest = 0

    
    good_enough = cos <= lowest + 0.15
    argsort = np.argsort(cos)
    
#     print("%3i -- %.1f -- %.1f"%(total, cos.min(), cos.max()))
    for i, iarg in enumerate(argsort):
        if ( not good_enough[iarg] ):
            continue
#         if ( i % 2 != 0 ):
#             continue
#         print("%.2f"%cos[iarg])

        with_z = looped_dreams[iarg][0]
        looping = looped_dreams[iarg][1]
        print(looped_dreams[iarg][1])
        if i == 0:
            plot_circles(idea, with_z, localer_strand_isheet, isheet_units, looping=looping)
            plt.show()
            print("debug on printing")

        final_looped_dreams.append((idea, base, with_z, looping, local_xy_str, local_strand_isheet, local_aligned_w_z ))
       
#     print("Debuggin!")
#     break
print("Looped dreams:", len(final_looped_dreams))


# In[57]:


unsupported_helix


# In[21]:


empty_binder_contig, target_contig = empty_contig.split()
motif_contig_parts = empty_binder_contig.split(",")[:-1]


# In[23]:


max_delta_helix_in_protein = args.max_delta_helix_in_protein
max_delta_touching_helix = args.max_delta_touching_helix
loop_provision = args.loop_provision
min_helix_size = args.min_helix_size
max_protein_size = args.max_protein_size
max_helix_size = args.max_helix_size
min_strand_size = args.min_strand_size
max_strand_size = args.max_strand_size

helix_step = args.helix_len_step
strand_step = 2

# ang_per_helix_res = 1.56
# ang_per_strand_res = 3.45
# helix_res_per_sheet_res = ang_per_strand_res / ang_per_helix_res
helix_res_per_sheet_res = 2

final_looped_gap_sizes = []

overkill_dreams = len(final_looped_dreams) / args.num_to_output

# print("debug")
most_overkill = 5
for idea, base, with_z, looping, local_xy_str, local_strand_isheet, local_aligned_w_z in final_looped_dreams:
    
#     print("debug")
#     print(len(idea))
#     if len(idea) != 8:
#         continue
        
        
    localer_strand_isheet = np.full(len(idea), np.nan)
    localer_strand_isheet[:len(local_strand_isheet)] = local_strand_isheet
    
    is_temple = None
    is_og_templ = base == -1
#     in_z_bounds = np.zeros((len(idea), 2))
#     in_z_bounds[is_templ] = z_bounds
    in_motif_sizes = np.zeros(len(idea), int)
    in_motif_sizes[is_og_templ] = motif_sizes
    in_motif_contig_parts = np.full(len(idea), None, dtype=object)
    in_motif_contig_parts[is_og_templ] = np.array(motif_contig_parts)
    in_motif_can_add = np.full((len(idea), 2), 1000, dtype=int)
    in_motif_can_add[is_og_templ] = motif_can_add
    
    arr_is_strand = ~np.isnan(localer_strand_isheet)
    arr_is_helix = np.isnan(localer_strand_isheet)
    arr_min_size = np.full(len(idea), min_helix_size)
    arr_min_size[arr_is_strand] = min_strand_size
    
    # allow helices to be longer than the motif extractor said if they will be smaller
    #  than the minimum helix size
    for i in range(1000):
        max_sizes = in_motif_sizes + in_motif_can_add.sum(axis=-1)
        if ( (max_sizes >= arr_min_size).all() ):
            break
        needs_fixing = np.where(max_sizes < arr_min_size)[0]
        for ifix in needs_fixing:
            in_motif_can_add[ifix,:] += 1
#         print("Fixing")
    
    
    
#     local_xy_str = local_xy_str[looping]
    localer_strand_isheet = localer_strand_isheet[looping]
    idea = idea[looping]
    with_z = with_z[looping]
    base = base[looping]
    is_og_templ = is_og_templ[looping]
    arr_is_strand = arr_is_strand[looping]
    arr_is_helix = arr_is_helix[looping]
    arr_min_size = arr_min_size[looping]
#     in_z_bounds = in_z_bounds[looping]
    in_motif_sizes =in_motif_sizes[looping]
    in_motif_contig_parts = in_motif_contig_parts[looping]
    in_motif_can_add = in_motif_can_add[looping]
    looping =  np.arange(len(looping), dtype=int) #looping[looping.copy()]
    
    ordered_motif_contig_parts = [x for x in in_motif_contig_parts if x is not None]
    
    all_by_dist = np.linalg.norm(idea[:,None] - idea[None,:], axis=-1)
    
    is_touching = (all_by_dist < 12) & (all_by_dist > 0)
    is_far = all_by_dist > 16
    
    this_max_protein_size = min(max_protein_size, len(idea)*args.max_protein_size_helix_per_res )
    
    num_loops = len(idea)-1
    res_in_loops = num_loops * loop_provision
    min_res_in_helices = min_helix_size * arr_is_helix.sum() + min_strand_size * arr_is_strand.sum()
    min_protein_size = res_in_loops + min_res_in_helices
    
    max_res_to_provision = max_protein_size - min_protein_size
    max_res_in_helices = max_protein_size - res_in_loops
    
    this_max_helix_size = min(max_helix_size, max_res_to_provision + min_helix_size )
    
    if ( max_res_to_provision < 0 ):
        continue
        
    arr_max_size = np.full(len(idea), this_max_helix_size)
    arr_max_size[arr_is_strand] = min(max_strand_size, this_max_helix_size)
    
#     print(max_res_to_provision)
    each_helix_potential_sizes = []
    for i in range(len(idea)):
        min_size = max( arr_min_size[i], in_motif_sizes[i] )
        max_size = max( arr_max_size[i], min_size )
        step = helix_step if arr_is_helix[i] else strand_step
        each_helix_potential_sizes.append(list(np.arange(min_size, max_size+0.01, step, dtype=int)))
#     plot_circles(idea, with_z, looping=looping)
#     plt.show()
    
    min_remaining_size = np.zeros(len(idea), int)
    for i in range(len(idea)):
        for j in range(i, len(idea)):
            min_remaining_size[i] += np.min(each_helix_potential_sizes[i])
            
    to_helix_len = np.ones(len(idea), float)
    to_helix_len[arr_is_strand] = helix_res_per_sheet_res
    
    def recursive_helix_lengths(
                    sizes,
                    is_touching,
                    to_helix_len,
                    each_helix_potential_sizes,
                    max_res_in_helices,
                    max_delta_touching_helix,
                    min_remaining_size
                    ):
        
        i_helix = len(sizes)
        cur_size = np.sum(sizes)
        if ( i_helix == len(is_touching) ):
            if ( cur_size + res_in_loops < args.min_protein_size ):
                return []
            return [sizes]
        
        minimum_to_add = min_remaining_size[i_helix]
        
        provisional_left = max_res_in_helices - cur_size - minimum_to_add
#         print(provisional_left)
        if ( provisional_left < 0 ):
            return []
        
        
#         print(cur_size, provisional_left)
        # deal with is touching
        # get previous sizes with is_touching
        # convert everything to be a helix for this section
        is_touching_prev = is_touching[i_helix,:i_helix]
        lb = 0
        ub = 1000
        if ( is_touching_prev.sum() > 0 ):
            prev_touching_sizes = np.array(sizes)[is_touching_prev] * to_helix_len[:len(sizes)][is_touching_prev]
            lb = prev_touching_sizes.max() - max_delta_touching_helix
            ub = prev_touching_sizes.min() + max_delta_touching_helix
#             print(lb, ub)
        if ( len(sizes) > 0 ):
            prev_sizes = np.array(sizes) * to_helix_len[:len(sizes)]
            lb_whole = prev_sizes.max() - max_delta_helix_in_protein
            ub_whole = prev_sizes.min() + max_delta_helix_in_protein
            lb = max(lb, lb_whole)
            ub = min(ub, ub_whole)
#             print(lb, ub)
        # convert helical lengths back to sheet lengths if it's a sheet
        lb = int(lb / to_helix_len[len(sizes)])
        ub = int(ub / to_helix_len[len(sizes)])
        
        results = []
        
        these_sizes = each_helix_potential_sizes[i_helix]
        min_size = these_sizes[0]
        for size in these_sizes:
            this_adds = size - min_size
            if ( this_adds > provisional_left ):
                continue
            if ( size < lb ):
                continue
            if ( size > ub ):
                continue
            
            results += recursive_helix_lengths(
                                    sizes + [size],
                                    is_touching,
                                    to_helix_len,
                                    each_helix_potential_sizes,
                                    max_res_in_helices,
                                    max_delta_touching_helix,
                                    min_remaining_size
                                    )
        return results
    
    results = recursive_helix_lengths(
                    [],
                    is_touching,
                    to_helix_len,
                    each_helix_potential_sizes,
                    max_res_in_helices,
                    max_delta_touching_helix,
                    min_remaining_size
                    )
    results = np.array(results)
#     print((results.sum(axis=-1) > 120).any())
    
#     print("Here", len(results))
    ###########
    
    # figure out who goes to which gaps
    num_gaps = is_og_templ.sum()+1
    gap_comprised_of = [[] for x in range(num_gaps)]
    full_diff_elems_in_gap = ['']*num_gaps
    helix_belongs_to = []
    loops_in_gap = np.zeros(num_gaps, int)
    cur_gap = 0
    for i in range(len(is_og_templ)):
        we_belong_to = [cur_gap]
        gap_comprised_of[cur_gap].append(i)
        if ( is_og_templ[i] ):
            cur_gap += 1
            gap_comprised_of[cur_gap].append(i)
            we_belong_to.append(cur_gap)
        else:
            ss = "H" if arr_is_helix[i] else "E"
            full_diff_elems_in_gap[cur_gap] += ss
#             full_diff_helices_in_gap[cur_gap] += 1
        if ( i < len(is_templ)-1 ):
            loops_in_gap[cur_gap] += 1
        helix_belongs_to.append(we_belong_to)
    
    extra_because_loops = np.zeros(num_gaps, int)
    extra_because_loops += loops_in_gap * loop_provision
            
    # decide how big each gap is going to be in terms of residues
    # does unique
    all_gap_sizes = []
    results = np.array(results)
    
    # there's an itertools product coming after this
    # there's no reason to have more than 3x as many as you need here
    overkill_results = len(results) * overkill_dreams
    if overkill_results > most_overkill: 
        print("Overkill results", overkill_results)
    if overkill_results > 3:
        new_size = int(max( 3 / overkill_dreams, 30))
        sizes = results.sum(axis=-1)
        size_bin = sizes//5
        from_each_bin = int(np.ceil( new_size / len(size_bin) ))
        arange = np.arange(len(results))
        keep_mask = np.zeros(len(results), bool)
        for binn in np.unique(size_bin):
            this_size = min(from_each_bin, (size_bin == binn).sum())
            keep_mask[ np.random.choice(arange[size_bin == binn], this_size, replace=False)] = True

#         indices = np.random.choice(np.arange(len(results)), new_size, replace=False)
        results = results[keep_mask]
    overkill_results = len(results) * overkill_dreams
    print("New Overkill results", overkill_results)
    
    for result in results:
        helix_to_provision = result - in_motif_sizes

        a_helix_was_too_long = False
        
        # helix provision options ends up being a list of vectors of length num_gaps
        # each entry says how that helix is going to divy up the residues it needs to add outside of it's input length
        helix_provision_options = []
        for i_helix in range(len(is_og_templ)):
            we_belong_to = helix_belongs_to[i_helix]
            this_to_prevision = helix_to_provision[i_helix]
            if (len(we_belong_to) == 1):
                these_provisions = np.zeros(num_gaps, int)
                these_provisions[we_belong_to[0]] = this_to_prevision
                helix_provision_options.append([these_provisions])
            else:
                can_add_l, can_add_r = in_motif_can_add[i_helix]
                left = np.arange(0, this_to_prevision+0.01, 1, dtype=int)
                right = this_to_prevision - left
                these_options = []
                for l, r in zip(left, right):
                    if ( l > can_add_l ):
                        continue
                    if ( r > can_add_r ):
                        continue
                    this_option = np.zeros(num_gaps, int)
                    this_option[we_belong_to[0]] = l
                    this_option[we_belong_to[1]] = r
                    these_options.append(this_option)
                if ( len(these_options) == 0 ):
                    a_helix_was_too_long = True
                helix_provision_options.append(these_options)
                
                
#         print('here')
        if ( a_helix_was_too_long ):
            continue

        # final size of each gap using all different combinations of helix provisions
        gap_sizes = np.array(list(np.sum(x, axis=0) for x in itertools.product(*helix_provision_options)))
        gap_sizes += extra_because_loops
        
#         search = np.array([32, 36, 26,  3])
#         if ( gap_sizes.shape[-1] == len(search)):
#             if ( (gap_sizes == search).all(axis=-1).any() ):
#                 print(result)
        
        gap_sizes = np.unique(gap_sizes, axis=0)
        all_gap_sizes.append(gap_sizes)
        
    if ( len(all_gap_sizes) == 0 ):
        print(0)
        continue
    all_gap_sizes = np.concatenate( all_gap_sizes )
    
    all_gap_sizes = np.unique(all_gap_sizes, axis=0)

    # drop gaps where elements differ by less than helix_step
    for i in range(len(all_gap_sizes)):
        if ( i >= len(all_gap_sizes)-1 ):
            break
        this_gap = all_gap_sizes[i]

        deltas = np.abs(all_gap_sizes - this_gap)
        all_close = (deltas < helix_step).all(axis=-1)

        all_gap_sizes = all_gap_sizes[~all_close]

    print(len(all_gap_sizes), "debug on printint")
    
#     plot_circles(idea, with_z, localer_strand_isheet, isheet_units, looping=looping)
#     plt.show()


    
    final_looped_gap_sizes.append((idea, base, with_z, looping, full_diff_elems_in_gap, 
                                                       ordered_motif_contig_parts, all_gap_sizes,
                                   localer_strand_isheet))
    
    
#     print("debug")
#     break
    ###########
#     break

    


# In[61]:





# In[24]:



def generate_ss(binder_parts, helix_size=8, terminal_helices=True, split_motif_size=0):
    full_size = binder_parts[-1][-1]
    
    ss = torch.ones(full_size)*3
    
    helix_bounds = {}
    
    for name, start, end in binder_parts:
        
        tp = name.split("_")[0]
        
        bound_lb = None
        bound_ub = None
        
        letter = None
        
        if ( tp[0] == "m" ):
            bound_lb = start
            bound_ub = end
            letter = tp[1]
        elif ( tp == "f" ):
            pass
        elif ( tp == "fq" or tp == 'fe' or tp == 'fh' ):
            
            if tp == 'fq':
                letter = name.split("_")[-1][0]
            else:
                letter = tp[1]
            
            size = end-start+1
            if ( size < helix_size and letter == 'h' ):
                bound_lb = start
                bound_ub = end
                print("Warning: %s is a small helix"%name)
            else:
                center = (start + end) / 2
                half_helix = (helix_size-1) / 2

                # always ends up centered, errors to earlier on odd-splits
                bound_lb = int(center-half_helix)
                bound_ub = int(center+half_helix)
            
        if ( terminal_helices and tp != 'f' ):
            if ( start == 1 ):
                if ( bound_lb is None ):
                    bound_lb = 1
                    bound_ub = end
#                 assert( not bound_lb is None )
                bound_lb = 1
            if ( end == full_size ):
                if ( bound_lb is None ):
                    bound_lb = start
                    bound_ub = full_size
#                 assert( not bound_lb is None )
                bound_ub = full_size
            
        if ( not bound_lb is None ):
            assert(letter in 'eh')
            ss_value = 0 if letter == 'h' else 1
            ss[bound_lb-1:bound_ub] = ss_value
            
        
        if ( split_motif_size > 0  and tp[0] == "m"):
            center = (bound_lb + bound_ub) / 2
            half_split = (split_motif_size-1) / 2
            split_start = int(center-half_split)
            split_end = int(center+half_split)
            print(split_start, split_end)
            ss[split_start-1:split_end] = 3
            
    
        helix_bounds[name] = (bound_lb, bound_ub)

    return ss, helix_bounds

def generate_adj(ss, helix_bounds, contacts, anti_contacts):
    
    full_size = len(ss)
    adj = torch.ones((full_size, full_size))*2
    
    for is_anti in [False, True]:
        cts = anti_contacts if is_anti else contacts
        store_val = 0 if is_anti else 1
        
        for contact in cts:
            name1, name2 = contact.split(":")

            lb1,ub1 = helix_bounds[name1]
            lb2,ub2 = helix_bounds[name2]

            assert(not lb1 is None)
            assert(not lb2 is None)

            adj[lb1-1:ub1,lb2-1:ub2] = store_val
            adj[lb2-1:ub2,lb1-1:ub1] = store_val
        
    is_mask = ss == 3
    adj[is_mask,:] = 0
    adj[:,is_mask] = 0
    
    return adj

# we're losing some of our helix-length variation because this evenly splits
def parse_contig_into_parts(contig, motif_is_strand, full_diff_elems_in_gap, sheet_factor=0.5):

    binder_contig = contig.split(" ")[0]

    contig_parts = binder_contig.split(",")[:-1]


    free_parts = 0
    motif_parts = 0

    seqpos = 1

    binder_parts = []

    for part in contig_parts:
        if ( part.startswith("A") ):
            motif_parts += 1
            
            letter = 'e' if motif_is_strand[motif_parts-1] else 'h'

            name = "m%s_%i"%(letter, motif_parts)

            m_start, m_end = part[1:].split("-")
            m_start = int(m_start)
            m_end = int(m_end)
            size = m_end-m_start+1


            start = seqpos
            end = start + size-1
            binder_parts.append((name, start, end))

        else:
            this_helices = full_diff_helices_in_gap[free_parts]
            free_parts += 1

            h = ""
            if len(this_helices) > 0:
                h = this_helices.lower() if len(this_helices) == 1 else "q"

            name = "f%s_%i"%(h, free_parts)

            size = int(part.split("-")[0])

            start = seqpos
            end = start + size-1

            if ( len(this_helices) <= 1 ):
                binder_parts.append((name, start, end))

            else:
                num_elems = len(this_helices)
                factors = np.ones(num_elems, float)
                for i, elem in enumerate(list(this_helices)):
                    if elem == 'E':
                        factors[i] = sheet_factor
                total_factor = factors.sum()
                size_per_factor = size / total_factor
                my_size = size_per_factor * factors
                my_start_offset = np.cumsum(my_size).astype(int)
                
                starts = np.array([start] + list(my_start_offset + start)[:-1])
                ends = (list(starts-1) + [end])[-len(this_helices):]
                
                for i, (this_start, this_end) in enumerate(zip(starts, ends)):
                    letter = this_helices[i].lower()
                    this_name = name + "_%s%i"%(letter, i+1)
                    binder_parts.append((this_name, this_start, this_end))



        seqpos += size
        
    assert(motif_parts == len(motif_is_strand))
    
    return binder_parts


# In[25]:


name_ss_adj_contig_groups = []

all_names = []
for i_idea, (idea, base, with_z, looping, full_diff_helices_in_gap, 
                             ordered_motif_contig_parts, all_gap_sizes,
            localer_strand_isheet) in enumerate(final_looped_gap_sizes):

    input_full_diff_helices_in_gap = list(full_diff_helices_in_gap)
    
    arr_is_strand = ~np.isnan(localer_strand_isheet)
    is_motif = base == -1
    motif_is_strand = arr_is_strand[is_motif]
    
    prefix="h%i_id%i"%(len(idea), i_idea)
    
    this_output = []
    for i_gap, these_gaps in enumerate(all_gap_sizes):
        full_diff_helices_in_gap = list(input_full_diff_helices_in_gap)
        assert(len(these_gaps)-1 == len(ordered_motif_contig_parts))
        
        
        # start by assembling a really dumb empty contig without motifs
        building_contig_parts = []
        needs_pop = []

        for i_part in range(len(these_gaps)):
            this_gap = these_gaps[i_part]
            if ( this_gap != 0 ):
                building_contig_parts.append("%i-%i"%(this_gap, this_gap))
            else:
                assert(i_part == 0 or i_part == len(these_gaps)-1)
                needs_pop.append(i_part)
            if ( i_part < len(ordered_motif_contig_parts) ):
                building_contig_parts.append(ordered_motif_contig_parts[i_part])
                
        for i_part in reversed(needs_pop):
            full_diff_helices_in_gap.pop(i_part)

        building_contig_parts.append("0")
        
        contig = ",".join(building_contig_parts) + " " + target_contig

        # next, list of parts that will form final protein
        binder_parts = parse_contig_into_parts(contig, motif_is_strand, full_diff_helices_in_gap)
#         print(binder_parts)
        
#         plot_circles(idea, with_z, localer_strand_isheet, isheet_units, looping=looping)
#         plt.show()
        
#         break
        
        part_is_circle = []
        i_circle = 0
        for i in range(len(binder_parts)):
            part_name = binder_parts[i][0]
            if ( "m" in part_name or "h" in part_name or 'e' in part_name):
                part_is_circle.append(i_circle)
                i_circle += 1
            else:
                part_is_circle.append(-1)
        assert(i_circle == len(idea))
        part_is_circle = np.array(part_is_circle)
                
        circle_name = []
        for i_circle in range(len(idea)):
            idx = list(part_is_circle).index(i_circle)
            circle_name.append(binder_parts[idx][0])
        
        all_by_dist = np.linalg.norm(idea[:,None] - idea[None,:], axis=-1)
        
        needs_close_cst = (all_by_dist < 9) & (all_by_dist > 0)
        needs_far_cst = (all_by_dist > 16)
        
        contacts = []
        anti_contacts = []
        
        for i, j in zip(*np.where(needs_close_cst)):
            if ( i > j ):
                continue
            if ( circle_name[i].startswith("m") and circle_name[j].startswith("m") ):
                continue
            contacts.append("%s:%s"%(circle_name[i], circle_name[j]))
            
        for i, j in zip(*np.where(needs_far_cst)):
            if ( i > j ):
                continue
            if ( circle_name[i].startswith("m") and circle_name[j].startswith("m") ):
                continue
            anti_contacts.append("%s:%s"%(circle_name[i], circle_name[j]))
        
        
        
        ss, helix_bounds = generate_ss(binder_parts)
        adj = generate_adj(ss, helix_bounds, contacts, anti_contacts)
        
        
#         plot_circles(idea, with_z, localer_strand_isheet, isheet_units, looping=looping)
#         plt.show()
#         print(ss)
        
#         sns.heatmap(adj)
#         break

        if args.no_motif_actually:
            contig = "%i-%i,0 "%(len(ss), len(ss)) + target_contig
        
        
        name = prefix + "_gap%i_s%03i"%(i_gap, binder_parts[-1][-1])
        
        all_names.append(name)
        this_output.append((name, ss, adj, contig))

        
    name_ss_adj_contig_groups.append(this_output)
    
    
#     has_q = False
#     for part in binder_parts:
#         if "q" in part[0]:
#             has_q = True
#             break
#     if has_q:
#         break


    


# In[38]:



num_to_output = args.num_to_output
num_per_idea = int(np.ceil(num_to_output / len(name_ss_adj_contig_groups)))

sizes = np.array(list(len(group) for group in name_ss_adj_contig_groups))

if ( args.dont_balance_lengths ):
    final_name_ss_adj_contig = []
    for group, size in zip(name_ss_adj_contig_groups, sizes):

        size = num_per_idea
        repeat_number = int(np.ceil(size / len(group)))

        keep_idx = np.random.choice(np.arange(len(group), dtype=int).repeat(repeat_number), size, replace=False)

        this_group = list(np.array(group, dtype=object)[keep_idx])

        final_name_ss_adj_contig += this_group
        
else:
    mega_list = []
    mega_lengths = []
    for group, size in zip(name_ss_adj_contig_groups, sizes):

        if len(group) == 0:
            continue
        
        lengths = np.array([len(x[1]) for x in group])
        size = num_per_idea
        repeat_number = int(np.ceil(size / len(group)))
        mega_list += list(np.array(group, dtype=object)[np.arange(len(group)).repeat(repeat_number)])
        mega_lengths += list(lengths.repeat(repeat_number))
    
    mega_list = np.array(mega_list, dtype=object)
    mega_lengths = np.array(mega_lengths)
    weights = np.zeros(len(mega_lengths))
    smallest = np.min(mega_lengths)
    biggest = np.max(mega_lengths)
    lb_sizes = np.arange(smallest, biggest, 5)
    accounted = 0
    for lb_size in lb_sizes:
        ub_size = lb_size + 5-1
        mask = (mega_lengths >= lb_size) & (mega_lengths <= ub_size)
        if mask.sum() > 0:
            accounted += 1
        
    for lb_size in lb_sizes:
        goal_frac = 1/accounted
        ub_size = lb_size + 5-1
        mask = (mega_lengths >= lb_size) & (mega_lengths <= ub_size)
        weight = goal_frac / mask.sum()
        weights[mask] = weight
#         print(mask.sum())
    
    assert(np.isclose(weights.sum(), 1))
    
    use_indices = torch.multinomial(torch.tensor(weights), num_to_output).numpy()
    
    
    final_name_ss_adj_contig = list(mega_list[use_indices])
        
    
    
    
import random
random.shuffle(final_name_ss_adj_contig)


# In[233]:



# fname_prefix = "/home/bcov/sc/design/tie2/attempt3/loop_it/go/scaffolds/my_test_scaff_arc"

output_dict = {}
output_list = []
for name, ss, adj, contig in final_name_ss_adj_contig:
    output_dict[name + "_ss"] = ss
    output_dict[name + "_adj"] = adj
    output_list.append("%s %s"%(name, contig))
    
with open(fname_prefix + ".txt", "w") as f:
    f.write("\n".join(output_list) + "\n")
    
with open(fname_prefix + ".pt", "wb") as f:
    torch.save(output_dict, f)


# In[234]:


fname_prefix
# print(contig.replace(" ", "\\ "))

# name = "my_test_scaff"

# with open("/home/bcov/sc/design/tie2/attempt3/loop_it/go/scaffolds/my_test_scaff/my_test_scaff_ss.pt", "wb") as f:
#     torch.save(ss, f)
    
# with open("/home/bcov/sc/design/tie2/attempt3/loop_it/go/scaffolds/my_test_scaff/my_test_scaff_adj.pt", "wb") as f:
#     torch.save(adj, f)


# In[22]:





# In[23]:


fname_prefix


# In[ ]:




