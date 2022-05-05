from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore    import *
from PyQt5.QtGui     import *
from PyQt5.QtWidgets import *
# from PyQt5.QtWidgets import QFileDialog 
import csv

class MyTabs(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)
        layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tabcoord = QWidget()

        # Add tabs
        self.tabs.addTab(self.tabcoord, "Coordinates")

        # Save Button
        self.buttonSave = QtWidgets.QPushButton('Save', self)
        self.buttonSave.clicked.connect(self.handleSave)
        self.buttonSave.clicked.connect(self.file_save)

        # Initiate Tables
        self.createTable()

        # Create tab
        self.tabcoord_layout = QVBoxLayout(self.tabcoord)
        self.tabcoord_layout.addWidget(self.tablewidget)
        self.tabcoord_layout.addWidget(self.buttonSave)

        # Add tabs to widget
        layout.addWidget(self.tabs)

    def createTable(self):
        # Table
        self.tablewidget = QTableWidget()
        self.tablewidget.setColumnCount(3)
        self.tablewidget.setRowCount(10)
        self.tablewidget.setHorizontalHeaderLabels(["Coordinates", "ROIs", "Circle Values"])
        self.tablewidget.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        # for i,row in enumerate(cur): ## <- insert self.ROI info
        #     self.tablewidget.setRowCount(self.tablewidget.rowCount())
        #     for j,val in enumerate(row):
        #         self.tablewidget.setItem(i, j, QtGui.QTableWidgetItem(str(val)))

    # def handleSave(self):
    #     with open('monschedule.csv', 'w') as stream:              
    #         writer = csv.writer(stream, lineterminator='\n')         
    #         for row in range(self.tablewidgetmon.rowCount()):
    #             rowdata = []
    #             for column in range(self.tablewidgetmon.columnCount()):
    #                 item = self.tablewidgetmon.item(row, column)
    #                 if item is not None:
    #                     rowdata.append(item.text())                  
    #                 else:
    #                     rowdata.append('')

    #             writer.writerow(rowdata)
    
    # def file_save(self):
    #     name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')
    #     file = open(name,'w')
    #     text = self.textEdit.toPlainText()
    #     file.write(text)
    #     file.close()



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    main = MyTabs()
    main.show()
    sys.exit(app.exec_())