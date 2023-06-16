
# This Python code is for evaluation of Relative Closest Branch Length (RCBL).
# The first argument should be original tree data in Newick format.
# The second should be original tree data in Newick format.
# The third should be output file name.
# [e.g.] python3 calc_rcbl.py original_tree.newick prune_tree.newick output.tsv

import sys
import math
import phylotreelib as pt

# Detect closest leaf
def detect_closest_leaf( target_leaf, leaves_list, tree ):
    closest_dist = math.inf
    closest_leaf = ''
    for leaf in leaves_list:
        if ( leaf != target_leaf ):
            dist = tree.nodedist( leaf, target_leaf )
            if ( dist < closest_dist ):
                closest_dist = dist
                closest_leaf = leaf

    return( closest_leaf, closest_dist )

print( 'Reading original tree ...' )
with pt.Newicktreefile( sys.argv[ 1 ] ) as treefile:
    original_tree = treefile.readtree()

print( 'Reading pruned tree ...' )
with pt.Newicktreefile( sys.argv[ 2 ] ) as treefile:
    pruned_tree = treefile.readtree()

print( 'Creating output file ...' )
fout = open( sys.argv[ 3 ], 'w' )
fout.write(
    'TargetLeaf'   + '\t' +
    'ClosestLeaf'  + '\t' +
    'BranchLength' + '\t' +
    'State'        + '\n'
)

print( 'Calculating closest branchs lengths of original tree ...\n' )
leaves_list       = original_tree.leaflist()
closest_dist_list = list()

iter = 0
for target_leaf in leaves_list:
    iter += 1
    #print( target_leaf )
    closest_leaf, closest_dist = detect_closest_leaf( target_leaf, leaves_list, original_tree )
    closest_dist_str = '%.5f' % round( closest_dist, 5 )
    print(
        str( iter )      + '\t' +
        target_leaf      + '\t' +
        closest_leaf     + '\t' +
        closest_dist_str
    )
    fout.write(
        target_leaf      + '\t' +
        closest_leaf     + '\t' +
        closest_dist_str + '\t' +
        'original'       + '\n'
    )
    closest_dist_list.append( closest_dist )

# Calculate average closest branch length of original tree
acbl_original = 0.0
sum           = 0.0
for closest_dist in closest_dist_list: sum += closest_dist
acbl_original = sum / len( closest_dist_list )

print( 'Calculating closest branchs lengths of pruned tree ...\n' )
# Clear and reset valiables
leaves_list.clear()
closest_dist_list.clear()
leaves_list       = pruned_tree.leaflist()
closest_dist_list = list()

iter = 0
for target_leaf in leaves_list:
    iter += 1
    #print( target_leaf )
    closest_leaf, closest_dist = detect_closest_leaf( target_leaf, leaves_list, pruned_tree )
    closest_dist_str = '%.5f' % round( closest_dist, 5 )
    print(
        str( iter )      + '\t' +
        target_leaf      + '\t' +
        closest_leaf     + '\t' +
        closest_dist_str
    )
    fout.write(
        target_leaf      + '\t' +
        closest_leaf     + '\t' +
        closest_dist_str + '\t' +
        'pruned'         + '\n'
    )
    closest_dist_list.append( closest_dist )

# Calculate average closest branch length of pruned tree
acbl_pruned = 0.0
sum         = 0.0
for closest_dist in closest_dist_list: sum += closest_dist
acbl_pruned = sum / len( closest_dist_list )

# Show result
print( '\n' )
print( 'ACBL of original tree : ' + str( acbl_original               ) )
print( 'ACBL of pruned   tree : ' + str( acbl_pruned                 ) )
print( 'RCBL                  : ' + str( acbl_pruned / acbl_original ) )

fout.close()

'''
print( 'Reading tree ...' )
with pt.Newicktreefile( sys.argv[ 1 ] ) as treefile:
    tree = treefile.readtree()

print( 'Calculating closest branch lengths using nearest_n_leaves...' )
leaves_list       = tree.leaflist()
closest_leaf      = ''
closest_dist      = 0.0
closest_leaf_list = list()
closest_dist_list = list()

iter = 0
for leaf in leaves_list:
    iter += 1
    closest_leaf = tree.nearest_n_leaves( leaf , 1 )
    closest_leaf = closest_leaf.pop()
    closest_dist = tree.nodedist( leaf, closest_leaf )
    print(
        str( iter )  + '\t' +
        leaf         + '\t' +
        closest_leaf + '\t' +
        str( round( closest_dist, 5 ) )
    )
    closest_leaf_list.append( closest_leaf )
    closest_dist_list.append( closest_dist )
'''