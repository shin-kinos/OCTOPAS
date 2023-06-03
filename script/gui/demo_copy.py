
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
    QFormLayout,
    QMainWindow
)
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui  import QIcon
from pathlib      import Path
#from PyQt6 import QHBoxLayout, QWindow, QMainWindow, QVBoxLayout

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

class comboboxForStopOption( QComboBox ):
    def __init__( self ):
        super().__init__()
        self.setStyleSheet( st.combobox )
        self.addItems( ['NL', 'RTL' ] )

class Example( QMainWindow ):
    def __init__( self ):
        super().__init__()
        self.initUI()

    def initUI( self ):
        vlayout = QVBoxLayout()
        hlayout = QHBoxLayout()
        widget = QWidget()
        widget.setLayout( hlayout )

        a1 = QLabel( 'label1' )
        a2 = QLabel( 'label2' )
        hlayout.addWidget( a1 )
        hlayout.addWidget( a2 )

        vlayout.addLayout( hlayout )
        self.show()

if __name__ == '__main__':
    app = QApplication( sys.argv )
    window = Example()
    window.show()
    app.exec()