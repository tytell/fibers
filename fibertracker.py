import sys, os
from PyQt5 import QtGui, QtCore, QtWidgets
import logging

from fibertracker.fiberwindow import FiberWindow

def main():
    logging.basicConfig(level=logging.DEBUG)

    app = QtWidgets.QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)

    fw = FiberWindow()
    fw.show()

    return app.exec_()

if __name__ == '__main__':
    sys.exit(main())
