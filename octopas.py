
'''
    OCTOPAS: Tool to reduce redundancy of phylogenetic trees by referencing parent nodes
    Copyright (C) 2023  Shintaro Kinoshita

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import os                 # For file operation
import sys                # For command-line args
import time               # For elapse time reporting
from   enum import Enum   # For Enum error handling
import phylotreelib as pt # Package for tree handling

VERSION = '0.2.4'      # Current version of program
ARG_S   = False        # Global variable of argument '--silent', True or False, default False
EXE_LOG = ''           # Global variable of execution log
TSV_LOG = ''           # Global variable of output TSV file
DIGIT   = 6            # Global variable to round float values

''' These global variables are not used anymore !!!
    RED     = '\033[1;31m' # Escape sequence colour code : RED
    GREEN   = '\033[1;32m' # Escape sequence colour code : GREEN
    RESET   = '\033[0m'    # Reset escape sequence colour code
    OUTPUT  = None         # Global variable of output files name prefix
    OUT_DIR = ''           # Glocal variable of output directory
    YELLOW  = '\033[1;33m' # Escape sequence colour code : YELLOW
    BLUE    = '\033[1;34m' # Escape sequence colour code : BLUE
    PURPLE  = '\033[1;35m' # Escape sequence colour code : PURPLE
    CYAN    = '\033[1;36m' # Escape sequence colour code : CYAN
'''

# --------------------------------------------------------------------------- #
#  GENERAL FUNCTIONS (EXECUTION LOGS, OPTIONS MENU AND ERROR HANDLING ETC.)   #
# --------------------------------------------------------------------------- #

# Report log and print it if ARG_S = False
def print_log( content ):
    global EXE_LOG, ARG_S
    EXE_LOG += content + '\n'
    if ( ARG_S == False ): print( content )

# Report log of TSV output
def report_tsv( content ):
    global TSV_LOG
    TSV_LOG += content + '\n'

# Print version
def print_version():
    global VERSION
    print( VERSION )
    sys.exit( 0 )

# Print help
def print_help( script_name ):
    print( '\nUsage: %', script_name, '[Options] \n\n\nOptions :\n' )
    print( '-i | --input      <String>     Input file name, REQUIRED' )
    print( '-o | --output     <String>     Output file name prefix, REQUIRED' )
    print( '-O | --outdir     <String>     Output directory name, if it is empty, all output files are saved in current directory' )
    print( '-f | --format     <String>     Format of input file:\n'           \
           '                                   * newick = Newick (default)\n' \
           '                                   * nexus  = NEXUS' )
    print( '-n | --NL         <Integer>    Number of leaves remain, 3 or more, default 3' )
    print( '-r | --RTL        <Real>       Threshold of relative tree length, at range of [0, 1], default 0.95' )
    print( '-F | --outfmt     <String>     Format of output file:\n'          \
           '                                   * newick = Newick (default)\n' \
           '                                   * nexus  = NEXUS' )
    print( '-R | --resolution <Integer>    Resolution of leaf pruning, 1 or more, default 1' )
    print( '-s | --silent     <Flag>       If true, program does not show any logs on terminal during calculation' )
    print( '-v | --version    <Flag>       Print current version, ignore all other arguments' )
    print( '-h | --help       <Flag>       Print this help, ignore all other arguments' )
    sys.exit( 0 )

# Error handling
# Enum of error cases
class Error( Enum ):
    NO_INPUT_FILE    = 'Input file name is required.'
    NO_OUTPUT_FILE   = 'Output file name is required.'
    NL_NOT_INT       = 'Option -n or --NL must be integer number.'
    RTL_NOT_REAL     = 'Option -r or --RTL must be real number.'
    RES_NOT_INT      = 'Option -R or --resolution must be integer number.'
    NL_TOO_SMALL     = 'Option -n or -NL must be 3 or more'
    RTL_OUT_RANGE    = 'Option -r or --RTL must be at range of [0, 1].'
    INFMT_INCORRECT  = 'Option -f or --format must be \'newick\' or \'nexus\'.'
    OUTFMT_INCORRECT = 'Option -F or --outfmt must be \'newick\' or \'nexus\'.'

# Error bombing function using Enum 'Error'
def error_bomb( MESSAGE ):
    # Colorised message is shown on terminal
    # And normal message is saved to log
    global ARG_S
    ARG_S = True

    print(     '\n\nERROR !!!' )
    print_log( '\n\nERROR !!!' )

    print(     MESSAGE.value )
    print_log( MESSAGE.value )

    print(     '\nPROGRAM HALTED !!!' )
    print_log( '\nRPOGRAM HALTED !!!' )

    sys.exit( 1 )

# Print 'PROGRAM FINISHED !!!'
def finish_program():
    global EXE_LOG, ARG_S
    EXE_LOG += '\n' + 'PROGRAM FINISHED !!!' + '\n'
    if ( ARG_S == False ): print( '\n' + 'PROGRAM FINISHED !!!' )

# --------------------------------------------------------------------------- #
#   GENERAL CLASSES (OPTION SETTING, FILE READING AND ERROR HANDLING ETC.)    #
# --------------------------------------------------------------------------- #

class Options:
    def __init__(
        self,
        args     = list(),   # List of command-line arguments
        input    = None,     # Input file name
        output   = None,     # Output file name
        outdir   = '',       # Output directory name
        format   = 'newick', # Input file format, default newick
        nl       = '3',      # Number of leaes remain, default 3
        rtl      = '0.95',   # Threshold of RTL, default 0.95
        outfmt   = 'newick', # Output file format, default newick
        res      = '1',      # Resolution of lraf pruning, default 1
        nl_flag  = False,    # Flag if '-n' appears in arguments
        rtl_flag = False,    # Flag if '-r' appears in arguments
        res_flag = False,    # Flag if '-R' appears in arguments
        use_nl   = False,    # Whether NL are used as stop option, default False
        silent   = False     # If program does not show logs
    ):
        self.args     = args
        self.input    = input
        self.output   = output
        self.outdir   = outdir
        self.format   = format
        self.nl       = nl
        self.rtl      = rtl
        self.outfmt   = outfmt
        self.res      = res
        self.nl_flag  = nl_flag
        self.rtl_flag = rtl_flag
        self.res_flag = res_flag
        self.use_nl   = use_nl
        self.silent   = silent

    def setOptions( self, argv ):
        # If no arguments, print help
        if ( len( argv ) == 1 ): print_help( argv[ 0 ] )
        # Get command-line arguments
        self.args = argv

        i = 1
        while ( i < len( argv ) ):
            if   ( argv[ i ] == '-i' or argv[ i ] == '--input'      ): self.input  = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-o' or argv[ i ] == '--output'     ): self.output = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-O' or argv[ i ] == '--outdir'     ): self.outdir = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-f' or argv[ i ] == '--format'     ): self.format = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-n' or argv[ i ] == '--NL'         ): self.nl     = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-r' or argv[ i ] == '--RTL'        ): self.rtl    = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-F' or argv[ i ] == '--outfmt'     ): self.outfmt = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-R' or argv[ i ] == '--resolution' ): self.res    = argv[ i + 1 ]; i += 2
            elif ( argv[ i ] == '-s' or argv[ i ] == '--silent'     ): self.silent = True; i += 1
            elif ( argv[ i ] == '-v' or argv[ i ] == '--version'    ): print_version()
            elif ( argv[ i ] == '-h' or argv[ i ] == '--help'       ): print_help( argv[ 0 ] )
            else                                                     : print_help( argv[ 0 ] )

    def checkOptions( self ):
        # Set output file names
        ''' These global variables are not used anymore !!!
            global OUTPUT, OUT_DIR
            OUTPUT  = self.output
            OUT_DIR = self.outdir
        '''

        # Set global variables
        global ARG_S
        ARG_S = self.silent

        # Check if input and output file names are given
        if ( self.input  is None ): error_bomb( Error.NO_INPUT_FILE  )
        if ( self.output is None ): error_bomb( Error.NO_OUTPUT_FILE )

        # Check if '--NL' is not integer
        if ( ( self.nl ).isdigit() == False ): error_bomb( Error.NL_NOT_INT )

        # Check if '--RTL' is not float
        if ( ( self.rtl ).replace( '.', '', 1 ).isdigit() == False ): error_bomb( Error.RTL_NOT_REAL )

        # Check if '--resolution' is not integer
        if ( ( self.res ).isdigit() == False ): error_bomb( Error.RES_NOT_INT )

        # Check if '--format' is 'newick' or 'nexus'
        if ( self.format != 'newick' and self.format != 'nexus' ): error_bomb( Error.INFMT_INCORRECT )

        # Check if '--outfmt' is 'newick' or 'nexus'
        if ( self.outfmt != 'newick' and self.outfmt != 'nexus' ): error_bomb( Error.OUTFMT_INCORRECT )

        # Check if '-n', '-r' or '-R' exist in arguments
        if ( '-n' in self.args or '--NL'         in self.args ): self.nl_flag  = True
        if ( '-r' in self.args or '--RTL'        in self.args ): self.rtl_flag = True
        if ( '-R' in self.args or '--resolution' in self.args ): self.res_flag = True
        #print( self.args )

        # Check if '--NL' is 3 or more
        if ( int( self.nl ) < 3 ): error_bomb( Error.NL_TOO_SMALL )

        # Check if '--RTL' is at range of [0, 1].
        if ( float( self.rtl ) <= 0.0 or float( self.rtl ) >= 1.0 ): error_bomb( Error.RTL_OUT_RANGE )

        # If '-R' is 1, use pruneWithDefaultMode()
        if ( self.res == '1' ): self.res_flag = False

        # If both of stop options appear in arguments, RTL is prioritised
        if   ( self.nl_flag == True  and self.rtl_flag == True  ): self.use_nl = False
        elif ( self.nl_flag == True  and self.rtl_flag == False ): self.use_nl = True
        elif ( self.nl_flag == False and self.rtl_flag == True  ): self.use_nl = False
        elif ( self.nl_flag == False and self.rtl_flag == False ): self.use_nl = False
        #print( 'self.use_nl =', self.use_nl )

    def showOptions( self ):
        # Make message for stop option
        stop_option = ''
        if   ( self.use_nl == True  ): stop_option = ' * Stop option (NL)       : ' + str( self.nl  )
        elif ( self.use_nl == False ): stop_option = ' * Stop option (RTL)      : ' + str( self.rtl )

        # Make message for the rest
        input_filename  = ' * Input filename         : ' + self.input
        output_filename = ' * Output filename prefix : ' + self.output
        output_dirname  = ' * Output dirname         : ' + self.outdir
        silent          = ' * Silent                 : ' + str( self.silent )

        # Show options
        print_log( '\nParameter set :'                                               )
        print_log( '===============================================================' )
        print_log(  input_filename                                                   )
        print_log(  output_filename                                                  )
        print_log(  output_dirname                                                   )
        print_log(  stop_option                                                      )
        print_log(  silent                                                           )
        print_log( '===============================================================' )

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
        print_log( '\nTotal elapsed time : ' + self.elapsed_time + ' sec.' )

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

# Detect pruned leaves at resolution = n
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
    def readTree( self, input_tree, file_format ):
        print_log( '\nReading tree data ...' )
        if ( file_format == 'newick'):
            with pt.Newicktreefile( input_tree ) as treefile:
                self.original_tree = treefile.readtree()
        elif ( file_format == 'nexus' ):
            with pt.Nexustreefile( input_tree ) as treefile:
                self.original_tree = treefile.readtree()
        #print( self.original_tree )
        print_log( '=> DONE' )

    # Get list of all leaves in original tree
    def getLeavesList( self ):
        print_log( '\nGetting leaf names ...' )
        self.all_leaves_list = ( self.original_tree ).leaflist()
        #print( self.all_leaves_list )
        print_log( '=> DONE' )

    # Mid-point rooting
    def midPointRooting( self ):
        print_log( '\nMid-point rooting ...' )
        ( self.original_tree ).rootmid()
        print_log( '=> DONE' )

    # Detect root node
    def detectRootNode( self ):
        print_log( '\nDetecting root node ...' )
        self.root_node = ( self.original_tree ).root
        print_log( '=> DONE' )

    # Calculate Statistics of input tree
    def calculateTreeStatistics( self ):
        print_log( '\nCalculating basic Statistics of the input tree ...' )
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
        print_log( '=> DONE' )

    # Show Statistics of input tree
    def showTreeStatistics( self ):
        # Make some variables String
        global DIGIT
        num_all_leaves   = str(        self.num_all_leaves           )
        rtl_denominator  = str( round( self.rtl_denominator, DIGIT ) )
        diameter         = str( round( self.diameter,        DIGIT ) )
        height           = str( round( self.height,          DIGIT ) )
        print_log( '\nInput tree Statistics :' )
        print_log( '==============================================================='   )
        print_log( ' * Number of the leaves in the tree          : ' + num_all_leaves  )
        print_log( ' * Total branch length of the tree           : ' + rtl_denominator )
        print_log( ' * Diameter (the longest leaf-leaf distance) : ' + diameter        )
        print_log( ' * Height (the longest root-to-tip distance) : ' + height          )
        print_log( '==============================================================='   )

class Pruner( Tree ):
    def __init__(
        self,
        pruned_tree       = None,   # Pruned Tree
        pruned_leaves     = list(), # Pruned leaves list
        remain_leaves     = list(), # Remaining leaves list
        stop_at_nl        = 3,      # Stop option of remaining leaves number, default 3
        stop_at_rtl       = 0.95,   # Stop option of RTL threshold, default 0.95
        current_rtl       = None,   # Current RTL in each iteration
        pruned_leaf       = None,   # Pruned one leaf in each iteration
        pruned_leaf_list  = list(), # Pruned [resolution] leaves in each iteration
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
        # Report TSV, header name
        report_tsv( 'Iteration\tLeavesNumber\tRTL\tPrunedLeaf' )
        # Get initial RTL
        self.current_rtl = 1.0

        print_log( '\nSTART PRUNING !\n' )
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
                print_log(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                            'RTL : '                + current_rtl_round             + '\t' + \
                            'Remaining leaves # : ' + str( self.num_remain_leaves ) + '\t' + \
                            'Pruned leaf : '        + str( self.pruned_leaf       ) )
                # Report TSV log
                report_tsv( str( self.iteration_time    ) + '\t' + \
                            str( self.num_remain_leaves ) + '\t' + \
                            current_rtl_round             + '\t' + \
                            self.pruned_leaf )
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
                    print_log(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                                'RTL : '                + current_rtl_round             + '\t' + \
                                'Remaining leaves # : ' + str( self.num_remain_leaves ) + '\t' + \
                                'Pruned leaf : '        + str( self.pruned_leaf       ) )
                    # Report TSV log
                    report_tsv( str( self.iteration_time    ) + '\t' + \
                                str( self.num_remain_leaves ) + '\t' + \
                                current_rtl_round             + '\t' + \
                                self.pruned_leaf )
                    # Stop loop if 3 leaves remain and program still try to run
                    if ( self.num_remain_leaves == 3 ):
                        print_log( '\nNOTE : Iteration stopped since there are only 3 leaves in the tree now.' )
                        break
                    # Inclement iteration time
                    self.iteration_time += 1

    # Pruning leaves at resolution = n
    def pruneWithFastMode( self, use_nl, resolution ):
        # Get number of all leaves
        self.num_remain_leaves = self.num_all_leaves
        # Report TSV, header name
        report_tsv( 'Iteration\tLeavesNumber\tRTL' )
        # Get initial RTL
        self.current_rtl = 1.0

        print_log( '\nSTART PRUNING !\n' )
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
                print_log(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                            'RTL : '                + current_rtl_round             + '\t' + \
                            'Remaining leaves # : ' + str( self.num_remain_leaves ) )
                # Stop loop if [# leaves remain] < [resolution] and program still try to run
                if ( ( self.num_remain_leaves - resolution ) <= 3 ):
                    print_log( '\nNOTE : Iteration stopped since no leaves can be pruned with the resolution.' )
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
                print_log(  'Iteration : '          + str( self.iteration_time    ) + '\t' + \
                            'RTL : '                + current_rtl_round             + '\t' + \
                            'Remaining leaves # : ' + str( self.num_remain_leaves ) )
                # Stop loop if [# leaves remain] < [resolution] and program still try to run
                if ( ( self.num_remain_leaves - resolution ) <= 3 ):
                    print_log( '\nNOTE : Iteration stopped since no leaves can be pruned with the resolution.' )
                    break
                # Inclement iteration time
                self.iteration_time += 1

# --------------------------------------------------------------------------- #
#                      CLASSES FOR OUTPUT FILES HANDLING                      #
# --------------------------------------------------------------------------- #

class Output():
    def __init__(
        self,
        out_prefix        = None,    # Prefix of output files
        out_dir_name      = '',      # Output directory name
        out_execute_log   = '',      # Execute log file name
        out_tsv_log       = '',      # TSV data file name
        out_remain_leaves = '',      # List of remaining leaves
        out_pruned_leaves = '',      # List of pruned leaves
        out_tree          = None,    # Output pruned tree
        fout_log          = None,    # Output file for exe log
        fout_tsv          = None,    # Output file for STV
        fout_remain       = None,    # Output file for remaining list
        fout_pruned       = None,    # Output file for pruned leaves list
        fout_tree         = None,    # Output file for pruned tree
        fout_tree_format  = 'newick' # Output file for pruned tree format, default newick
    ):
        self.out_prefix        = out_prefix
        self.out_dir_name      = out_dir_name
        self.out_execute_log   = out_execute_log
        self.out_tsv_log       = out_tsv_log
        self.out_remain_leaves = out_remain_leaves
        self.out_pruned_leaves = out_pruned_leaves
        self.out_tree          = out_tree
        self.fout_log          = fout_log
        self.fout_tsv          = fout_tsv
        self.fout_remain       = fout_remain
        self.fout_pruned       = fout_pruned
        self.fout_tree         = fout_tree
        self.fout_tree_format  = fout_tree_format

    def setOutputFileNames( self, output, outdir, outfmt ):
        #global OUTPUT, OUT_DIR
        self.out_prefix   = output
        self.out_dir_name = outdir

        # Create output directory if necessary
        if ( self.out_dir_name != '' and os.path.exists( self.out_dir_name ) == False ):
            if ( self.out_dir_name[ -1 ] != '/' ): self.out_dir_name += '/'
            print_log( '\nOutput directory ' + self.out_dir_name + ' was created.\n' )
            os.mkdir( self.out_dir_name )
        elif ( self.out_dir_name != '' and os.path.exists( self.out_dir_name ) == True ):
            if ( self.out_dir_name[ -1 ] != '/' ): self.out_dir_name += '/'

        # Record output file format and extension
        self.fout_tree_format = outfmt
        extension             = outfmt

        # Set output files name
        self.out_execute_log   = self.out_dir_name + self.out_prefix + '.log'
        self.out_tsv_log       = self.out_dir_name + self.out_prefix + '.tsv'
        self.out_remain_leaves = self.out_dir_name + self.out_prefix + '_remaining_leaves.txt'
        self.out_pruned_leaves = self.out_dir_name + self.out_prefix + '_pruned_leaves.txt'
        self.out_pruned_tree   = self.out_dir_name + self.out_prefix + '_tree.' + extension

    # Save execution log
    def saveExeLog( self ):
        global EXE_LOG
        self.fout_log = open( self.out_execute_log, 'w' )
        ( self.fout_log ).write( EXE_LOG )

    # Save TSV file
    def saveTsv( self ):
        global TSV_LOG
        self.fout_tsv = open( self.out_tsv_log, 'w' )
        ( self.fout_tsv ).write( TSV_LOG )
        print_log( '\n' + self.out_tsv_log + ' was saved.' )

    # Save remaining leaves list
    def saveRemainLeavesList( self, remain_leaves ):
        self.fout_remain = open( self.out_remain_leaves, 'w' )
        for leaf in remain_leaves: ( self.fout_remain ).write( leaf + '\n' )
        print_log( self.out_remain_leaves + ' was saved.' )

    # Save pruned leaves list
    def savePrunedLeavesList( self, pruned_leaves ):
        self.fout_pruned = open( self.out_pruned_leaves, 'w' )
        for leaf in pruned_leaves: ( self.fout_pruned ).write( leaf + '\n' )
        print_log( self.out_pruned_leaves + ' was saved.' )

    # Save output tree
    def saveTree( self, pruned_tree ):
        if ( self.fout_tree_format == 'newick' ):
            with open( self.out_pruned_tree, 'w' ) as outfile:
                outfile.write( pruned_tree.newick() )
        elif ( self.fout_tree_format == 'nexus' ):
            with open( self.out_pruned_tree, 'w' ) as outfile:
                outfile.write( pruned_tree.nexus() )
        print_log( self.out_pruned_tree + ' was saved.' )

    # Close output files
    def closeFiles( self ):
        ( self.fout_log    ).close()
        ( self.fout_tsv    ).close()
        ( self.fout_remain ).close()
        ( self.fout_pruned ).close()

    # Show output file names
    def showOutputFileNames( self ):
        print_log( '\nOutput file name :' )
        print_log( '===============================================================' )
        print_log( ' * Execution log            : ' + self.out_execute_log           )
        print_log( ' * Output TSV               : ' + self.out_tsv_log               )
        print_log( ' * List of remaining leaves : ' + self.out_remain_leaves         )
        print_log( ' * List of pruned leaves    : ' + self.out_pruned_leaves         )
        print_log( ' * Pruned tree              : ' + self.out_pruned_tree           )
        print_log( '===============================================================' )

# --------------------------------------------------------------------------- #
#                                MAIN FUNCTION                                #
# --------------------------------------------------------------------------- #

def main():
    # Elapsed time : START
    time = Time()
    time.startTime()

    # Set options
    options = Options()
    options.setOptions( sys.argv )
    options.checkOptions()
    options.showOptions()

    # Set output files info
    output = Output()
    output.setOutputFileNames( options.output, options.outdir, options.outfmt )

    # Get phylogenetic tree
    tree = Pruner()
    tree.readTree( options.input, options.format )
    tree.getLeavesList()
    tree.midPointRooting()
    tree.detectRootNode()
    tree.calculateTreeStatistics()
    tree.showTreeStatistics()

    # Now pruning time! 😁
    tree.setInitialTreeInfo()
    tree.setStopOptions( int( options.nl ), float( options.rtl ) )
    if   ( options.res_flag == True  ): tree.pruneWithFastMode(    options.use_nl, int( options.res ) )
    elif ( options.res_flag == False ): tree.pruneWithDefaultMode( options.use_nl )

    # Save output files
    output.saveTsv()
    output.saveRemainLeavesList( tree.remain_leaves )
    output.savePrunedLeavesList( tree.pruned_leaves )
    output.saveTree( tree.pruned_tree )

    # Show output file names
    output.showOutputFileNames()

    finish_program()

    # Elapsed time : END
    time.endTime()
    time.calcElapsedTime()
    time.showElapsedTime()

    # Save execution logs and close files
    output.saveExeLog()
    output.closeFiles()

if __name__ == '__main__':
    main()

''' This code block is not used anymore !!!
    print( 'Reading tree data ...' )
    with pt.Newicktreefile( sys.argv[ 1 ] ) as treefile:
        tree = treefile.readtree()

    print( 'Mid-point rooting ...' )
    tree.rootmid()

    print( 'Getting initial info ...' )
    rootnode         = tree.root
    total_branch_len = tree.length()
    rtl              = 1.0
    threshold        = float( sys.argv[ 2 ] )

    leaves = list()

    print( 'Getting leaf names ...' )
    leaves = tree.leaflist()
    #print( leaves )

    def get_target_leaf( leaves ):
        target_leaf    = leaves[ 0 ] # initialise
        target_leaf_bl = tree.nodedist( leaves[ 0 ], rootnode ) # initialise
        for leaf in leaves:
            path    = tree.nodepath( leaf, rootnode )
            leaf_bl = tree.nodedist( path[ 0 ], path[ 1 ] )
            if ( leaf_bl < target_leaf_bl ):
                target_leaf    = leaf
                target_leaf_bl = leaf_bl

        return target_leaf

    iter = 1
    while( rtl > threshold ):
        target_leaf = get_target_leaf( leaves )
        tree.remove_leaf( target_leaf )
        leaves.remove( target_leaf )
        rtl = tree.length() / total_branch_len
        print( 'iter:', iter, '\tRTL:', rtl, '\tpruned:', target_leaf )
        iter += 1
'''
