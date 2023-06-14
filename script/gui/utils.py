
# Input read newick file
def read_newick( file_path ):
    fin = open( file_path, 'r' )
    file_content = fin.readline()
    return file_content

# Check options
def check_options(
    stop_option,
    threshold,
    resolution
):
    #print( stop_option )
    #print( threshold   )
    #print( resolution  )

    if ( stop_option == 'Relative tree length' and threshold.replace( '.', '', 1 ).isdigit() == False ):
        return( 'error', 'Threshold of RTL must be real number.' )

    if ( stop_option == 'Relative tree length' and threshold.replace( '.', '', 1 ).isdigit() == True ):
        if ( float( threshold ) < 0.0 or float( threshold ) > 1.0 ):
            return( 'error', 'Threshold of RTL must be at range of [0, 1].' )

    if ( stop_option == '# of leaves remain' and threshold.isdigit() == False ):
        return( 'error', 'Threshold of remaining leaves number must be integer.' )

    if ( stop_option == '# of leaves remain' and threshold.isdigit() == True ):
        if ( int( threshold ) < 3 ):
            return( 'error', 'Threshold of leaves number must be 3 or more.' )

    if ( resolution.isdigit() == False ):
        return( 'error', 'Resolution must be integer.' )

    if ( int( resolution ) < 1 ):
        return( 'error', 'Resolution must be 1 or more.' )

    return( 'ok', '' )