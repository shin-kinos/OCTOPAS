
import sys
from PyQt6.QtWidgets import (
    QApplication,
    QWidget,
    QLineEdit,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QFileDialog,
    QLabel,
    QHBoxLayout,
    QGridLayout,
    QCheckBox,
    QComboBox,
    QFormLayout
)
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui  import QIcon
from pathlib      import Path

import stylesheet_for_qt as st

WINDOW_WIDTH      = 800
WINDOW_HEIGHT     = 500
FIN_BUTTON_WIDTH  = 170
FIN_BUTOON_HEIGHT = 30
LINK_TEMPLATE     = '<a href={0}>{1}</a>'
PYTHON_SRC_LINK   = 'https://www.cranfield.ac.uk/people/dr-zoltan-kevei-337915'

# Class for title label
class titleLabel( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setText( 'PROTOPYPE OF GUI FOR NEW TREE PRUNER' )
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

# Class for combobox of stop option
class comboboxForStopOption( QComboBox ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.combobox )
        self.addItems( [ '# of leaves remain', 'Relative tree length' ] )

# Class for option message style
class optionMessageLabel( QLabel ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.option_msg )

class MyApp( QWidget ):
    def __init__( self ):
        super().__init__()
        global                \
            WINDOW_WIDTH,     \
            WINDOW_HEIGHT,    \
            FIN_BUTTON_WIDTH, \
            FIN_BUTOON_HEIGHT

        # Create outer layout
        layout = QGridLayout()   # Omajinai
        self.setLayout( layout ) # Omajinai

        # Name of title
        self.setWindowTitle( 'PROTOPYPE OF GUI FOR NEW TREE PRUNER' )

        # Width and height
        self.resize( WINDOW_WIDTH, WINDOW_HEIGHT )

        # Create title layout
        title = titleLabel()

        # Set Python source link label layout
        source = pythonSourceLinkLabel()

        # Set input file button
        button_fin = QPushButton( '&Select tree file ...', clicked = self.openInputFile  )
        button_fin.setStyleSheet( st.button_fin )
        button_fin.setFixedSize( QSize( FIN_BUTTON_WIDTH, FIN_BUTOON_HEIGHT ) )

        # Checkbox of using example data set
        example_cbox = QCheckBox( ' Use example data?' )
        example_cbox.setStyleSheet( st.example_cbox )

        # Set input text field 
        self.input_content = inputContentLabel()

        # Set stop option label
        option_stop_msg = optionMessageLabel()
        option_stop_msg.setText( 'Stop option :' )

        # Set combobox for stop option
        option_stop_opt = comboboxForStopOption()
        option_stop_opt.setFixedSize( QSize( 230, FIN_BUTOON_HEIGHT ) )

        # Set stop option threshold
        option_thresh_msg = optionMessageLabel()
        option_thresh_msg.setText( 'Stop option threshold :' )

        # Set RUN button
        button_run = QPushButton( '&RUN' )
        button_run.setStyleSheet( st.button_run )
        button_run.setFixedSize( QSize( FIN_BUTTON_WIDTH, FIN_BUTOON_HEIGHT ) )

        layout.addWidget( title,              0, 0, 1, 1 )
        layout.addWidget( source,             1, 0, 1, 1 )
        layout.addWidget( button_fin,         2, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( example_cbox,       3, 0, 1, 1 )
        layout.addWidget( self.input_content, 4, 0, 1, 1 )
        layout.addWidget( QLabel( 'Options' ), 5, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( option_stop_msg,    6, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignLeft   )
        layout.addWidget( option_stop_opt,    6, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )
        layout.addWidget( option_thresh_msg,  7, 0, 1, 1 )
        layout.addWidget( button_run,         8, 0, 1, 1, alignment = Qt.AlignmentFlag.AlignCenter )

    def openInputFile( self ):
        file_name, _ = QFileDialog.getOpenFileName( self )
        if file_name:
            file_path = Path( file_name )
            ( self.input_content ).setText( str( file_path ) )

app = QApplication( sys.argv )

window = MyApp()
window.setStyleSheet( st.main_flame )
window.show()

app.exec()