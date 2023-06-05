
import os                 # For file operation
import sys                # For command-line args
import time               # For elapse time reporting
from   enum import Enum   # For Enum error handling
import phylotreelib as pt # Package for tree handling

print( 'Reading original tree ...' )
with pt.Newicktreefile( sys.argv[ 1 ] ) as treefile:
    original_tree = treefile.readtree()

print( 'Creating output file ...' )
fout = open( sys.argv[ 1 ] + '_output.tsv', 'w' )
fout.write(
    'TargetLeaf'   + '\t' +
    'ClosestLeaf'  + '\t' +
    'BranchLength' + '\t' +
    'State'        + '\n'
)

print( 'Calculating closest branch lengths ...' )
leaves_list       = original_tree.leaflist()
closest_leaf      = ''
closest_dist      = 0.0
closest_leaf_list = list()
closest_dist_list = list()

for leaf in leaves_list:
    closest_leaf = original_tree.nearest_n_leaves( leaf , 1 )
    closest_leaf = closest_leaf.pop()
    closest_dist = original_tree.nodedist( leaf, closest_leaf )

    closest_leaf_list.append( closest_leaf )
    closest_dist_list.append( closest_dist )

print( '\ndistance of leaves to their closest ones:' )
for i in range( 0, len( leaves_list ) ):
    closest_dist = round( closest_dist_list[ i ], 5 )
    print( leaves_list[ i ], '\t:', closest_dist, '(', closest_leaf_list[ i ], ')' )
    fout.write(
        leaves_list[ i ]       + '\t' +
        closest_leaf_list[ i ] + '\t' +
        str( closest_dist )    + '\t' +
        'original'             + '\n'
    )

print( '\nCalculating average ...' )
sum = 0.0
for elem in closest_dist_list: sum += elem
original_average = sum / len( closest_dist_list )
print( '=>', round( original_average, 5 ) )

print( '\nReading pruned tree ... ' )
with pt.Newicktreefile( sys.argv[ 2 ] ) as treefile:
    pruned_tree = treefile.readtree()

print( 'Calculating closest branch lengths ...' )
leaves_list       = pruned_tree.leaflist()
closest_leaf      = ''
closest_dist      = 0.0
closest_leaf_list.clear()
closest_dist_list.clear()
for leaf in leaves_list:
    closest_leaf = pruned_tree.nearest_n_leaves( leaf , 1 )
    closest_leaf = closest_leaf.pop()
    closest_dist = pruned_tree.nodedist( leaf, closest_leaf )

    closest_leaf_list.append( closest_leaf )
    closest_dist_list.append( closest_dist )

print( '\ndistance of leaves to their closest ones:' )
for i in range( 0, len( leaves_list ) ):
    closest_dist = round( closest_dist_list[ i ], 5 )
    print( leaves_list[ i ], '\t:', closest_dist, '(', closest_leaf_list[ i ], ')' )
    fout.write(
        leaves_list[ i ]       + '\t' +
        closest_leaf_list[ i ] + '\t' +
        str( closest_dist )    + '\t' +
        'pruned'               + '\n'
    )

print( '\nCalculating average ...' )
sum = 0.0
for elem in closest_dist_list: sum += elem
pruned_average = sum / len( closest_dist_list )
print( '=>', round( pruned_average, 5 ) )

print( '\nDintance score =', round( pruned_average / original_average, 5 ) )

fout.close()