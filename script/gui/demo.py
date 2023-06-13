
import sys
from PyQt6.QtWidgets import (
    QApplication,
    QWidget,
    QLineEdit,
    QPushButton,
    QTextEdit,
    QFileDialog,
    QLabel,
    QGridLayout,
    QCheckBox,
    QComboBox,
    QMessageBox
)
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui  import QIcon
from pathlib      import Path

import stylesheet_for_qt as st
import pruner
import utils

# -------------------------- WIDGETS STYLE PARAMETERS -------------------------- #
WINDOW_WIDTH      = 900
WINDOW_HEIGHT     = 700
FIN_BUTTON_WIDTH  = 200
FIN_BUTTON_HEIGHT = 30
OPT_BUTTON_WIDTH  = 230
OPT_BUTTON_HEIGHT = 30
RUN_BUTTON_WIDTH  = 180
RUN_BUTTON_HEIGHT = 30
DL_BUTTON_WIDTH   = 300
DL_BUTTON_HEIGHT  = 30
LINK_TEMPLATE     = '<a href={0}>{1}</a>'
PYTHON_SRC_LINK   = 'https://www.cranfield.ac.uk/people/dr-zoltan-kevei-337915'
# ------------------------------------------------------------------------------ #

# Class for title label
class titleLabel( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setText( 'PROTOTYPE OF GUI FOR NEW TREE PRUNER' )
        self.setStyleSheet( st.title )
        self.setAlignment( Qt.AlignmentFlag.AlignCenter )

# Class for python source link
class pythonSourceLinkLabel( QLabel ):
    def __init__( self, parent = None ):
        super().__init__()
        global LINK_TEMPLATE, PYTHON_SRC_LINK
        self.setStyleSheet( st.source )
        self.setOpenExternalLinks( True )
        self.setParent( parent )
        self.setText( LINK_TEMPLATE.format( PYTHON_SRC_LINK, 'Python source' ) )
        self.setAlignment( Qt.AlignmentFlag.AlignCenter )

# Class for input text field label
class inputContentLabel( QTextEdit ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.input_txt_field )

# Class for checkbox message
class checkboxMessageLabel( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setText( 'Use example data? : ' )
        self.setStyleSheet( st.checkbox_msg )

# Class for option title
class optionTitle( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.option_title )
        self.setText( 'Options' )

# Class for oresults title
class resultTitle( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.option_title )
        self.setText( 'Results' )

# Class for combobox of stop option
class comboboxForStopOption( QComboBox ):
    def __init__( self ):
        super().__init__()
        global OPT_BUTTON_WIDTH, FIN_BUTTON_HEIGHT
        self.setStyleSheet( st.combobox )
        self.addItems( [ 'Relative tree length', '# of leaves remain' ] )
        self.setFixedSize( QSize( OPT_BUTTON_WIDTH, FIN_BUTTON_HEIGHT ) )

# Class for options
class textBoxForOptions( QLineEdit ):
    def __init__( self ):
        super().__init__()
        global OPT_BUTTON_WIDTH, FIN_BUTTON_HEIGHT
        self.setStyleSheet( st.linetextbox )
        self.setAlignment( Qt.AlignmentFlag.AlignCenter )
        self.setFixedSize( QSize( OPT_BUTTON_WIDTH, FIN_BUTTON_HEIGHT ) )


# Class for option message style
class optionMessageLabel( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.option_msg )

# Class for result text
class textResult( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setText( 'test' )

# Class for download dialog
class showDownloadDialog( QMessageBox ):
    def __init__( self, file_name ):
        super().__init__()
        self.setWindowTitle( 'DOWNLOADING OUPUT' )
        self.setText( file_name + ' is downloaded!' )
        self.exec()

# Class for option error dialog
class showOptionError( QMessageBox ):
    def __init__( self, message ):
        super().__init__()
        self.setWindowTitle( 'OPTION ERROR' )
        self.setText( 'OPTION ERROR: ' + message )
        self.exec()

# Class of main central part of app
class mainApp( QWidget ):
    def __init__( self ):
        super().__init__()
        global                 \
            WINDOW_WIDTH,      \
            WINDOW_HEIGHT,     \
            FIN_BUTTON_WIDTH,  \
            FIN_BUTTON_HEIGHT, \
            RUN_BUTTON_WIDTH,  \
            RUN_BUTTON_HEIGHT, \
            DL_BUTTON_WIDTH,   \
            DL_BUTTON_HEIGHT

        # Initialise input tree output file contents
        self.input_tree    = ''
        self.output_leaves = ''
        self.output_tree   = ''

        # Create grid layout
        layout = QGridLayout()   # Omajinai
        self.setLayout( layout ) # Omajinai

        # Name of title
        self.setWindowTitle( 'PROTOTYPE OF GUI FOR NEW TREE PRUNER' )

        # Width and height
        self.resize( WINDOW_WIDTH, WINDOW_HEIGHT )

        # Create title layout
        title = titleLabel()

        # Set Python source link label layout
        source = pythonSourceLinkLabel()

        # Set input file button
        button_fin = QPushButton( '&Select Newick file ...', clicked = self.openInputFile )
        button_fin.setStyleSheet( st.button_fin )
        #button_fin.setIcon( QIcon( './icon/upload.icon' ) )
        button_fin.setIconSize( QSize( 24, 24 ) )
        button_fin.setFixedSize( QSize( FIN_BUTTON_WIDTH, FIN_BUTTON_HEIGHT ) )
        button_fin.setToolTip( 'This is a tooltip for the QPushButton widget' )

        # Checkbox of using example data set
        example_cbox = QCheckBox( ' Use example data?' )
        example_cbox.setStyleSheet( st.example_cbox )

        # Set input text field 
        self.input_content = inputContentLabel()

        # Set option title
        option_title = optionTitle()

        # Set stop option label
        option_stop_msg = optionMessageLabel()
        option_stop_msg.setText( 'Stop option :' )

        # Set combobox for stop option
        self.option_stop_opt = comboboxForStopOption()

        # Set stop option threshold label
        option_thresh_msg = optionMessageLabel()
        option_thresh_msg.setText( 'Stop option threshold :' )

        # Set stop option threshold 
        self.option_thresh = textBoxForOptions()
        self.option_thresh.setText( '0.95' )

        # Set resolution label
        option_resol_msg = optionMessageLabel()
        option_resol_msg.setText( 'Resolution :' )

        # Set resolution option threshold 
        self.option_resol = textBoxForOptions()
        self.option_resol.setText( '1' )

        # Set RUN button
        self.button_run = QPushButton( '&RUN', clicked = self.runProgram )
        #self.button_run.setStyleSheet( st.button_run )
        self.button_run.setFixedSize( QSize( RUN_BUTTON_WIDTH, RUN_BUTTON_HEIGHT ) )
        self.button_run.setEnabled( False )
        if ( self.button_run.isEnabled() == False ): self.button_run.setStyleSheet( st.button_run_disable )

        # Set result text box
        self.text_result = textResult()

        # Set result title (style is exactly same as option title)
        result_title = resultTitle()

        # Set results download (remaining leaves) button
        self.button_dl_leaves = QPushButton( '&Download remaining leaves list ...', clicked = self.downloadLeavesList )
        self.button_dl_leaves.setStyleSheet( st.button_download_able )
        self.button_dl_leaves.setFixedSize( QSize( DL_BUTTON_WIDTH, DL_BUTTON_HEIGHT ) )
        self.button_dl_leaves.setEnabled( False )
        if ( self.button_dl_leaves.isEnabled() == False ): self.button_dl_leaves.setStyleSheet( st.button_download_disable )

        # Set results download (pruned tree) button
        self.button_dl_tree = QPushButton( '&Download pruned tree ...', clicked = self.downloadPrunedTree )
        self.button_dl_tree.setStyleSheet( st.button_download_able )
        self.button_dl_tree.setFixedSize( QSize( DL_BUTTON_WIDTH, DL_BUTTON_HEIGHT ) )
        self.button_dl_tree.setEnabled( False )
        if ( self.button_dl_tree.isEnabled() == False ): self.button_dl_tree.setStyleSheet( st.button_download_disable )

        layout.addWidget( title,                  0, 0, 1, 1 )
        layout.addWidget( source,                 1, 0, 1, 1 )
        layout.addWidget( button_fin,             2, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft )
        layout.addWidget( example_cbox,           3, 0, 1, 1 )
        layout.addWidget( self.input_content,     4, 0, 1, 1 )
        layout.addWidget( option_title,           5, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( option_stop_msg,        6, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( self.option_stop_opt,   6, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( option_thresh_msg,      7, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( self.option_thresh,     7, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( option_resol_msg,       8, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( self.option_resol,      8, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( self.button_run,        9, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( QLabel( '' ),          10, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( result_title,          11, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( self.button_dl_leaves, 12, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( self.button_dl_tree,   13, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )

    def openInputFile( self ):
        file_name, _ = QFileDialog.getOpenFileName( self )
        if file_name:
            file_path = Path( file_name )
            # Read file and save content
            file_content = utils.read_newick( file_path )
            ( self.input_content ).setText( file_content )
            # Activate run button
            self.button_run.setEnabled( True )
            if ( self.button_run.isEnabled() == True ): self.button_run.setStyleSheet( st.button_run_able )
            # Disable download buttons
            self.button_dl_leaves.setEnabled( False )
            self.button_dl_tree.setEnabled(   False )
            if ( self.button_dl_leaves.isEnabled() == False ): self.button_dl_leaves.setStyleSheet( st.button_download_disable )
            if ( self.button_dl_tree.isEnabled()   == False ): self.button_dl_tree.setStyleSheet(   st.button_download_disable )

    def runProgram( self ):
        # Check parameters
        check, which = utils.check_options(
            ( self.option_stop_opt ).currentText(),
            ( self.option_thresh   ).text(),
            ( self.option_resol    ).text()
        )
        if   ( check == 'error' ): showOptionError( which )
        elif ( check == 'ok' ): # If 'check' == ok, run program
            # Disable run button immediately not to be pushed again
            self.button_run.setEnabled( False )
            if ( self.button_run.isEnabled() == False ): self.button_run.setStyleSheet( st.button_run_disable )
            # Run pruner
            self.output_leaves, self.output_tree = pruner.run_pruner(
                ( self.input_content ).toPlainText(),
                ( self.option_stop_opt ).currentText(),
                ( self.option_thresh   ).text(),
                ( self.option_resol    ).text()
            )
            # Activate download buttons
            self.button_dl_leaves.setEnabled( True )
            self.button_dl_tree.setEnabled( True )
            if ( self.button_dl_leaves.isEnabled() == True ): self.button_dl_leaves.setStyleSheet( st.button_download_able )
            if ( self.button_dl_tree.isEnabled()   == True ): self.button_dl_tree.setStyleSheet(   st.button_download_able )

    def downloadLeavesList( self ):
        file_path = ''
        content   = self.output_leaves
        dir_name  = QFileDialog.getExistingDirectory( self )
        if dir_name:
            file_path = dir_name + '/result_leaves.txt'
            file_path = Path( file_path )
            fout      = open( file_path, 'w' )
            fout.write( content )
            fout.close()
            showDownloadDialog( 'result_leaves.txt' )

    def downloadPrunedTree( self ):
        file_path = ''
        content   = self.output_tree
        dir_name  = QFileDialog.getExistingDirectory( self )
        if dir_name:
            file_path = dir_name + '/result_tree.newick'
            file_path = Path( file_path )
            fout      = open( file_path, 'w' )
            fout.write( content )
            fout.close()
            showDownloadDialog( 'result_tree.newick' )

app = QApplication( sys.argv )

window = mainApp()
window.setStyleSheet( st.main_flame )
window.show()

app.exec()