
import time               # For elapse time reporting
import phylotreelib as pt # Package for tree handling

DIGIT   = 6            # Global variable to round float values
RED     = '\033[1;31m' # Escape sequence colour code : RED
GREEN   = '\033[1;32m' # Escape sequence colour code : GREEN
RESET   = '\033[0m'    # Reset escape sequence colour code

# --------------------------------------------------------------------------- #
#                  GENERAL FUNCTIONS (FINISH PROGRAM REPORT)                  #
# --------------------------------------------------------------------------- #

# Print 'PROGRAM FINISHED !!!'
def finish_program():
    print( '\n' + GREEN + 'PROGRAM FINISHED !!!' + RESET )

# --------------------------------------------------------------------------- #
#             GENERAL CLASSES (OPTION SETTING AND TIME HANDLING)              #
# --------------------------------------------------------------------------- #

class Options:
    def __init__(
        self,
        nl     = 3,        # Number of leaves remain, default 3
        rtl    = 0.95,     # Threshold of RTL, default 0.95
        res    = 1,        # Resolution of leaf pruning, default 1
        use_nl = False,    # Whether NL are used as stop option, default False
        #format = 'newick', # Input file format, default newick
        #outfmt = 'newick', # Output file format, default newick
    ):
        self.nl       = nl
        self.rtl      = rtl
        self.res      = res
        self.use_nl   = use_nl
        #self.format   = format
        #self.outfmt   = outfmt

    def setOptions(
        self,
        stop_option,
        threshold,
        resolution
    ):
        # Set stop option and its threshold
        if ( stop_option == '# of leaves remain' ):
            self.use_nl = True
            self.nl     = int( threshold )
        elif ( stop_option == 'Relative tree length' ):
            self.use_nl = False
            self.rtl    = float( threshold )

        # Set resolution
        self.res = int( resolution )

    def showOptions( self ):
        # Make message for stop option
        stop_option = ''
        if   ( self.use_nl == True  ): stop_option = ' * Stop option (NL)  : ' + str( self.nl  )
        elif ( self.use_nl == False ): stop_option = ' * Stop option (RTL) : ' + str( self.rtl )

        # Make message for resolution
        resolution  = ' * Resolution        : ' + str( self.res )

        # Show options
        print( '\nParameter set :'                                               )
        print( '===============================================================' )
        print(  stop_option                                                      )
        print(  resolution                                                       )
        print( '===============================================================' )

# Elapsed time handling
class Time:
    def __init__(
        self,
        start_time   = None, # Start time
        end_time     = None, # End time
        elapsed_time = None  # Total time
    ):
        self.start_time   = start_time
        self.end_time     = end_time
        self.elapsed_time = elapsed_time

    def startTime( self ):
        self.start_time = time.time()

    def endTime( self ):
        self.end_time = time.time()

    def calcElapsedTime( self ):
        self.elapsed_time = self.end_time - self.start_time

    def showElapsedTime( self ):
        self.elapsed_time = '%.3f' % round( self.elapsed_time, 3 )
        print( '\nTotal elapsed time : ' + self.elapsed_time + ' sec.' )

# --------------------------------------------------------------------------- #
#                    FUNCTIONS FOR PHYLOGENETIC TREE PRUNING                  #
# --------------------------------------------------------------------------- #

# Detect one leaf pruned
def detect_pruned_leaf( tree, leaves, root ):
    # initialise target leaf
    target_leaf = leaves[ 0 ]
    # initialise target leaf branch length
    target_leaf_bl = tree.nodedist( leaves[ 0 ], root )

    for leaf in leaves:
        path    = tree.nodepath( leaf, root )
        leaf_bl = tree.nodedist( path[ 0 ], path[ 1 ] )
        if ( leaf_bl < target_leaf_bl ):
            target_leaf    = leaf
            target_leaf_bl = leaf_bl

    return target_leaf

def detect_pruned_leaves( tree, leaves, root, resolution ):
    # get leaves branch length
    leaves_bl = list()
    for leaf in leaves:
        path    = tree.nodepath( leaf, root )
        leaf_bl = tree.nodedist( path[ 0 ], path[ 1 ] )
        leaves_bl.append( leaf_bl )

    # sort branch lengths
    for i in range( 0, len( leaves ) - 1 ):
        for j in range( i + 1, len( leaves ) ):
            if ( leaves_bl[ j ] < leaves_bl[ i ] ):
                tmp            = leaves_bl[ i ]
                leaves_bl[ i ] = leaves_bl[ j ]
                leaves_bl[ j ] = tmp

                tmp         = leaves[ i ]
                leaves[ i ] = leaves[ j ]
                leaves[ j ] = tmp

    # Return first [resolution] elements
    return leaves[ : resolution ]

# --------------------------------------------------------------------------- #
#                    CLASSES FOR PHYLOGENETIC TREE PRUNING                    #
# --------------------------------------------------------------------------- #

class Tree:
    def __init__(
        self,
        original_tree      = None,   # Original phylogenetic tree
        all_leaves_list    = list(), # Leaves list in original tree
        num_all_leaves     = None,   # Number of all leaves
        rtl_denominator    = None,   # Total branch length of original tree
        root_node          = None,   # Root of the tree
        diameter           = None,   # Diameter of original tree
        height             = None    # Height of tree (longest root-to-tip distance)
    ):
        self.original_tree   = original_tree
        self.all_leaves_list = all_leaves_list
        self.num_all_leaves  = num_all_leaves
        self.rtl_denominator = rtl_denominator
        self.root_node       = root_node
        self.diameter        = diameter
        self.height          = height

    # Read phylogenetic tree and get information
    def readTree( self, input_tree ):
        print( '\nReading tree data ...' )
        self.original_tree = pt.Tree.from_string( input_tree )
        print( '=> DONE' )

    # Get list of all leaves in original tree
    def getLeavesList( self ):
        print( '\nGetting leaf names ...' )
        self.all_leaves_list = ( self.original_tree ).leaflist()
        #print( self.all_leaves_list )
        print( '=> DONE' )

    # Mid-point rooting
    def midPointRooting( self ):
        print( '\nMid-point rooting ...' )
        ( self.original_tree ).rootmid()
        print( '=> DONE' )

    # Detect root node
    def detectRootNode( self ):
        print( '\nDetecting root node ...' )
        self.root_node = ( self.original_tree ).root
        print( '=> DONE' )

    # Calculate Statistics of input tree
    def calculateTreeStatistics( self ):
        print( '\nCalculating basic Statistics of the input tree ...' )
        # Get number of leaves
        self.num_all_leaves = len( self.all_leaves_list )
        #print( self.num_all_leaves )

        # Get total branch length as denominator of RTL
        self.rtl_denominator = ( self.original_tree ).length()
        #print( self.rtl_denominator )

        # Get original tree diameter
        self.diameter = ( self.original_tree ).diameter()
        #print( self.diameter )

        # Get original tree height
        self.height = ( self.original_tree ).height()
        #print( self.height )

        #for leaf in self.all_leaves_list: print( leaf )
        print( '=> DONE' )

    # Show Statistics of input tree
    def showTreeStatistics( self ):
        # Make some variables String
        global DIGIT
        num_all_leaves   = str(        self.num_all_leaves           )
        rtl_denominator  = str( round( self.rtl_denominator, DIGIT ) )
        diameter         = str( round( self.diameter,        DIGIT ) )
        height           = str( round( self.height,          DIGIT ) )
        print( '\nInput tree Statistics :' )
        print( '==============================================================='   )
        print( ' * Number of the leaves in the tree          : ' + num_all_leaves  )
        print( ' * Total branch length of the tree           : ' + rtl_denominator )
        print( ' * Diameter (the longest leaf-leaf distance) : ' + diameter        )
        print( ' * Height (the longest root-to-tip distance) : ' + height          )
        print( '==============================================================='   )

class Pruner( Tree ):
    def __init__(
        self,
        pruned_tree       = None,   # Pruned Tree
        pruned_leaves     = list(), # Pruned leaves list
        remain_leaves     = list(), # Remaining leaves list
        stop_at_nl        = 3,      # Stop option of remaining leaves number, default 3
        stop_at_rtl       = 0.95,   # Stop option of RTL criterion, default 0.95
        current_rtl       = None,   # Current RTL in each iteration
        pruned_leaf       = None,   # Pruned one leaf in each iteration
        pruned_leaf_list  = list(), # Pruned [resolution] leavs in each iteration
        num_remain_leaves = None,   # Number of remaining leaves
        iteration_time    = 1       # Iteration time of pruning process
    ):
        super().__init__()
        self.pruned_tree       = pruned_tree
        self.pruned_leaves     = pruned_leaves
        self.remain_leaves     = remain_leaves
        self.stop_at_nl        = stop_at_nl
        self.stop_at_rtl       = stop_at_rtl
        self.current_rtl       = current_rtl
        self.pruned_leaf       = pruned_leaf
        self.pruned_leaf_list  = pruned_leaf_list
        self.num_remain_leaves = num_remain_leaves
        self.iteration_time    = iteration_time

    # Copy 'self.original_tree' into 'self.pruned_tree'
    def setInitialTreeInfo( self ):
        self.pruned_tree   = self.original_tree
        self.remain_leaves = self.all_leaves_list

    # Get stop option information from 'Options' class
    def setStopOptions( self, stop_at_nl, stop_at_rtl ):
        self.stop_at_nl  = stop_at_nl
        self.stop_at_rtl = stop_at_rtl

    # Pruning leaves at resolution = 1
    def pruneWithDefaultMode( self, use_nl ):
        # Get number of all leaves
        self.num_remain_leaves = self.num_all_leaves
        # Get Initial RTL
        self.current_rtl = 1.0

        print( '\nSTART PRUNING !\n' )
        # If stop option = NL
        if ( use_nl == True ):
            while ( self.num_remain_leaves > self.stop_at_nl ):
                self.pruned_leaf = detect_pruned_leaf( self.pruned_tree, self.remain_leaves, self.root_node )
                ( self.pruned_tree   ).remove_leaf(  self.pruned_leaf ) # Prune one leaf
                ( self.remain_leaves ).remove(       self.pruned_leaf ) # Pop target pruned leaf
                ( self.pruned_leaves ).append(       self.pruned_leaf ) # Add target pruned leaf
                # Calculate current RTL
                self.current_rtl = ( self.pruned_tree ).length() / self.rtl_denominator
                # Show number of remaining leaves
                self.num_remain_leaves -= 1
                current_rtl_round = '%.12f' % round( self.current_rtl, 12 )
                print(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                        'RTL : '                + current_rtl_round             + '\t' + \
                        '# of leaves remain : ' + str( self.num_remain_leaves ) + '\t' + \
                        'Pruned leaf : '        + str( self.pruned_leaf       ) )
                # Inclement iteration time
                self.iteration_time += 1

        # If stop option = RTL
        elif ( use_nl == False ):
                while ( self.current_rtl > self.stop_at_rtl ):
                    self.pruned_leaf = detect_pruned_leaf( self.pruned_tree, self.remain_leaves, self.root_node )
                    ( self.pruned_tree   ).remove_leaf(  self.pruned_leaf ) # Prune one leaf
                    ( self.remain_leaves ).remove(       self.pruned_leaf ) # Pop target pruned leaf
                    ( self.pruned_leaves ).append(       self.pruned_leaf ) # Add target pruned leaf
                    # Calculate current RTL
                    self.current_rtl = ( self.pruned_tree ).length() / self.rtl_denominator
                    # Show number of remaining leaves
                    self.num_remain_leaves -= 1
                    current_rtl_round = '%.12f' % round( self.current_rtl, 12 )
                    print(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                            'RTL : '                + current_rtl_round             + '\t' + \
                            '# of leaves remain : ' + str( self.num_remain_leaves ) + '\t' + \
                            'Pruned leaf : '        + str( self.pruned_leaf       ) )
                    # Stop loop if 3 leaves remain and program still try to run
                    if ( self.num_remain_leaves == 3 ):
                        print( '\nNOTE : Iteration stopped since there are only 3 leaves in the tree now.' )
                        break
                    # Inclement iteration time
                    self.iteration_time += 1

    # Pruning leaves at resolution = n
    def pruneWithFastMode( self, use_nl, resolution ):
        # Get number of all leaves
        self.num_remain_leaves = self.num_all_leaves
        # Get Initial RTL
        self.current_rtl = 1.0

        print( '\nSTART PRUNING !\n' )
        # If stop option = NL
        if ( use_nl == True ):
            while ( self.num_remain_leaves > self.stop_at_nl ):
                self.pruned_leaf_list = detect_pruned_leaves( self.pruned_tree,
                                                              self.remain_leaves,
                                                              self.root_node,
                                                              resolution )
                for leaf in self.pruned_leaf_list: ( self.pruned_tree   ).remove_leaf( leaf ) # Prune [resolution] leaves
                for leaf in self.pruned_leaf_list: ( self.remain_leaves ).remove(      leaf ) # Pop target pruned leaf
                for leaf in self.pruned_leaf_list: ( self.pruned_leaves ).append(      leaf ) # Add target pruned leaves
                # Calculate current RTL
                self.current_rtl = ( self.pruned_tree ).length() / self.rtl_denominator
                # Show number of remaining leaves
                self.num_remain_leaves -= resolution
                current_rtl_round = '%.12f' % round( self.current_rtl, 12 )
                print(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                        'RTL : '                + current_rtl_round             + '\t' + \
                        '# of leaves remain : ' + str( self.num_remain_leaves ) )
                # Stop loop if [# leaves remain] < [resolution] and program still try to run
                if ( ( self.num_remain_leaves - resolution ) <= 3 ):
                    print( '\nNOTE : Iteration stopped since no leaves can be pruned with the resolution.' )
                    break
                # Inclement iteration time
                self.iteration_time += 1

        # If stop option = RTL
        if ( use_nl == False ):
            while ( self.current_rtl > self.stop_at_rtl ):
                self.pruned_leaf_list = detect_pruned_leaves( self.pruned_tree,
                                                              self.remain_leaves,
                                                              self.root_node,
                                                              resolution )
                for leaf in self.pruned_leaf_list: ( self.pruned_tree   ).remove_leaf( leaf ) # Prune [resolution] leaves
                for leaf in self.pruned_leaf_list: ( self.remain_leaves ).remove(      leaf ) # Pop target pruned leaf
                for leaf in self.pruned_leaf_list: ( self.pruned_leaves ).append(      leaf ) # Add target pruned leaves
                # Calculate current RTL
                self.current_rtl = ( self.pruned_tree ).length() / self.rtl_denominator
                # Show number of remaining leaves
                self.num_remain_leaves -= resolution
                current_rtl_round = '%.12f' % round( self.current_rtl, 12 )
                print(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                        'RTL : '                + current_rtl_round             + '\t' + \
                        '# of leaves remain : ' + str( self.num_remain_leaves ) )
                # Stop loop if [# leaves remain] < [resolution] and program still try to run
                if ( ( self.num_remain_leaves - resolution ) <= 3 ):
                    print( '\nNOTE : Iteration stopped since no leaves can be pruned with the resolution.' )
                    break
                # Inclement iteration time
                self.iteration_time += 1

# --------------------------------------------------------------------------- #
#                      CLASSES FOR OUTPUT FILES HANDING                       #
# --------------------------------------------------------------------------- #

class Output():
    def __init__(
        self,
        out_remain_leaves_str = '',      # List of remaining leaves as string
        out_tree_str          = '',      # Output file for pruned tree
        #out_tree_format       = 'newick' # Output file for pruned tree format, default newick
    ):

        self.out_remain_leaves_str = out_remain_leaves_str
        self.out_tree_str          = out_tree_str
        #self.out_tree_format       = out_tree_format

    def LeavesListToString( self, remain_Leaevs ):
        # Concat all leaves into string with '\n'
        for leaf in remain_Leaevs:
            self.out_remain_leaves_str += leaf + '\n'
        return( self.out_remain_leaves_str )

    def PrunedTreeToString( self, pruned_tree ):
        return( pruned_tree.newick() )

# --------------------------------------------------------------------------- #
#                                CORE FUNCTION                                #
# --------------------------------------------------------------------------- #

def run_pruner(
    input_tree,
    stop_option,
    threshold,
    resolution
):
    # Elapsed time : START
    time = Time()
    time.startTime()

    # Set options
    options = Options()
    options.setOptions(
        stop_option,
        threshold,
        resolution
    )
    #options.checkOptions()
    options.showOptions()

    # Get phylogenetic tree
    tree = Pruner()
    tree.readTree( input_tree )
    tree.getLeavesList()
    tree.midPointRooting()
    tree.detectRootNode()
    tree.calculateTreeStatistics()
    tree.showTreeStatistics()

    # Now pruning time! 😁
    tree.setInitialTreeInfo()
    tree.setStopOptions( options.nl, options.rtl )
    if   ( options.res  > 1 ): tree.pruneWithFastMode(    options.use_nl, options.res )
    elif ( options.res == 1 ): tree.pruneWithDefaultMode( options.use_nl )

    output = Output()
    out_remain_leaves_str = output.LeavesListToString( tree.remain_leaves )
    #print( out_remain_leaves_str )
    out_pruned_tree_str   = output.PrunedTreeToString( tree.pruned_tree   )
    #print( out_pruned_tree_str )

    # Finish program log
    finish_program()

    # Elapsed time : END
    time.endTime()
    time.calcElapsedTime()
    time.showElapsedTime()

    return( out_remain_leaves_str, out_pruned_tree_str )
