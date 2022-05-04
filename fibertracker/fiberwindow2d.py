import os, sys
from PyQt5 import QtGui, QtCore, QtWidgets
import numpy as np
import logging
import h5py
import skimage.transform
from copy import copy
import pyqtgraph as pg
import csv

SETTINGS_FILE = "fibertracker.ini"

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class ImageClick(pg.ImageItem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    sigDoubleClick = QtCore.pyqtSignal(object)

    def hoverEvent(self, ev):
        ev.acceptClicks(button=0)

    def mouseClickEvent(self, ev):
        # logging.debug("Mouse click: {}, double {}".format(ev, ev.double()))
        if ev.double():
            self.sigDoubleClick.emit(ev)


class CircleROIwLine(pg.CircleROI):
    def __init__(self, radonWnd=None, **kwargs):
        super().__init__(**kwargs)

        self.radonWnd = radonWnd

        if "angle" in kwargs:
            self._lineAngle = np.deg2rad(kwargs["angle"])
        else:
            self._lineAngle = 0.0

        pos = self.pos()
        r = self.size()[0]/2
        ctrx = pos[0] + r
        ctry = pos[1] + r

        self._dx = r * np.cos(self._lineAngle)
        self._dy = r * np.sin(self._lineAngle)

        self.line = pg.PlotDataItem(x=[ctrx-self._dx, ctrx+self._dx], y=[ctry-self._dy, ctry+self._dy], color='g')

        self.setAcceptedMouseButtons(QtCore.Qt.LeftButton)

        self.sigRegionChanged.connect(self._updateLine)
        self.sigRegionChangeFinished.connect(self.updateRadon)
        self.sigClicked.connect(self.updateRadon)

    def setLineAngle(self, angle):
        self._lineAngle = np.deg2rad(angle)
        self._updateLine(self)

    def addToPlot(self, plot):
        plot.addItem(self)
        plot.addItem(self.line)

    def _updateLine(self, circ):
        pos = self.pos()
        r = self.size()[0]/2
        ctrx = pos[0] + r
        ctry = pos[1] + r

        self._dx = r * np.cos(self._lineAngle)
        self._dy = r * np.sin(self._lineAngle)

        self.line.setData(x=[ctrx-self._dx, ctrx+self._dx], y=[ctry-self._dy, ctry+self._dy])

    def updateRadon(self):
        if self.radonWnd is not None:
            self.radonWnd.updateROI(self)


class RadonWindow(QtWidgets.QWidget):
    def __init__(self, imageData, imageItem):
        super().__init__()

        self.curROI = None
        self.imageData = imageData
        self.imageItem = imageItem

        self.initUI()
        self.readSettings()

    def initUI(self):
        self.setWindowTitle("Radon")

        self.zoomPlot = pg.PlotWidget()
        self.zoomImg = pg.ImageItem()
        self.zoomPlot.addItem(self.zoomImg)

        self.zoomPlot.setAspectLocked()

        self.radonPlot = pg.PlotWidget()
        self.radonImg = pg.ImageItem()
        self.radonPlot.addItem(self.radonImg)

        self.anglePlot = self.radonPlot.plot(x=[], y=[], z=10, pen='r')

        ctr = QtCore.QPointF(32.0, 32.0)
        self.angleMark = [pg.InfiniteLine(pos=0, angle=90, pen='g'),
                          pg.InfiniteLine(pos=ctr, angle=0, pen='g')]
        self.radonPlot.addItem(self.angleMark[0])
        self.zoomPlot.addItem(self.angleMark[1])
        self.angleMark[0].setVisible(False)
        self.angleMark[1].setVisible(False)

        hlayout = QtWidgets.QHBoxLayout()
        hlayout.addWidget(self.zoomPlot)
        hlayout.addWidget(self.radonPlot)

        self.setLayout(hlayout)

    def updateROI(self, roi):
        self.curROI = roi
        self.updateRadon()

    def updateImage(self, imageData):
        self.imageData = imageData
        self.updateRadon()

    def updateRadon(self):
        roi = self.curROI

        if roi is None:
            return

        roiImg = roi.getArrayRegion(self.imageData, self.imageItem)

        radius = roi.size()[0] / 2
        c0, c1 = np.ogrid[0:roiImg.shape[0], 0:roiImg.shape[1]]
        reconstruction_circle = ((c0 - roiImg.shape[0] // 2) ** 2
                                 + (c1 - roiImg.shape[1] // 2) ** 2)

        outside = np.logical_and(roiImg == 0, reconstruction_circle >= (radius-1) ** 2)
        inside = np.logical_not(outside)

        roiImg -= np.mean(roiImg[inside])
        roiImg[outside] = 0.0 ## check

        roiImg /= np.max(roiImg) - np.min(roiImg)

        self.zoomImg.setImage(roiImg, autoLevels=True)

        theta = np.linspace(0, 180, max(roiImg.shape), endpoint=False)
        sinogram = skimage.transform.radon(roiImg, theta, circle=True)

        # if there were no lines, any point on the sinogram would be equal to diameter of the reconstruction circle
        # times the mean intensity

        angleintensity = np.sum(np.abs(sinogram), axis=0)
        angleintensity /= np.max(angleintensity)

        self.radonImg.setImage(sinogram.T)
        self.radonImg.setRect(QtCore.QRectF(0, 0, 180, roiImg.shape[1]))

        self.anglePlot.setData(x=theta, y=angleintensity * sinogram.shape[1])

        angle = theta[np.argmax(angleintensity)]

        roi.setLineAngle(angle)

        self.angleMark[0].setPos(angle)

        ctr = QtCore.QPointF(roiImg.shape[0]/2.0, roiImg.shape[1]/2.0)
        self.angleMark[1].setPos(ctr)
        self.angleMark[1].setAngle(angle)

        self.angleMark[0].setVisible(True)
        self.angleMark[1].setVisible(True)

        self.show()

    def readSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("ZoomWindow")

        self.resize(settings.value("size", QtCore.QSize(200, 200)))
        self.move(settings.value("position", QtCore.QPoint(1000, 200)))

        settings.endGroup()

    def writeSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("ZoomWindow")
        settings.setValue("size", self.size())
        settings.setValue("position", self.pos())
        settings.endGroup()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.writeSettings()
        event.accept()


class FiberWindow2d(QtWidgets.QWidget):
    def __init__(self, datafile, axis='yz'):
        super().__init__()

        self.ROIs = []

        self.datafile = datafile
        self.axis = axis

        self.h5file = h5py.File(datafile, 'r')
        self.data = self.h5file['image']

        self.shape0 = self.data.shape
        self.axisorder = []
        for a1 in axis.lower():
            if a1 == 'x':
                self.axisorder.append(0)
            elif a1 == 'y':
                self.axisorder.append(1)
            elif a1 == 'z':
                self.axisorder.append(2)
            else:
                raise ValueError("Unrecognized axis {}".format(a1))

        lastax = np.setdiff1d([0, 1, 2], self.axisorder)
        self.axisorder = np.concatenate((self.axisorder, lastax))

        # use this to undo the axis ordering
        self.axisunorder = np.argsort(self.axisorder)

        self.shape = np.array(self.shape0)[self.axisorder][:2]
        self.len = self.shape0[self.axisorder[2]]

        self.curFrame = int(self.len/2)

        self.activeROI = None
        self.activeLine = None

        self.initUI()
        self.readSettings()

        self.updateImage()

    def initUI(self):
        self.scrollBar = QtWidgets.QScrollBar(orientation=QtCore.Qt.Horizontal)
        self.scrollBar.setRange(0, self.len-1)
        self.scrollBar.setValue(self.curFrame)
        self.scrollBar.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                                           QtWidgets.QSizePolicy.Fixed))

        self.doneButton = QtWidgets.QPushButton("Done")
        self.saveButton = QtWidgets.QPushButton("Save") # creating saveButton widget

        self.frameNumberBox = QtWidgets.QSpinBox()
        self.frameNumberBox.setRange(0, self.len-1)

        self.frameNumberBox.setValue(self.curFrame)
        lab1 = QtWidgets.QLabel("&Frame number:")
        lab1.setBuddy(self.frameNumberBox)
        lab1.setAlignment(QtCore.Qt.AlignRight)

        hlayout1 = QtWidgets.QHBoxLayout()
        hlayout1.addWidget(self.scrollBar)
        hlayout1.addWidget(lab1)
        hlayout1.addWidget(self.frameNumberBox)
        hlayout1.addWidget(self.doneButton)
        hlayout1.addWidget(self.saveButton) # adding saveButton

        self.plot = pg.PlotWidget()
        self.image = ImageClick(image=self.getImage())
        self.image.sigDoubleClick.connect(self.doubleClick)

        self.plot.addItem(self.image)

        self.plot.setAspectLocked()
        rect = self.plot.viewRect()
        self.plot.setRange(rect=rect, padding=None)

        vlayout1 = QtWidgets.QVBoxLayout()
        vlayout1.addWidget(self.plot)
        vlayout1.addLayout(hlayout1)

        self.setLayout(vlayout1)

        # set up slots
        self.frameNumberBox.valueChanged.connect(self.scrollBar.setValue)
        self.scrollBar.valueChanged.connect(self.frameNumberBox.setValue)

        self.scrollBar.valueChanged.connect(self.setFrame)

        self.saveButton.clicked.connect(self.saveVal) # signal when saveButton clicked

        self.doneButton.clicked.connect(self.close)

        self.radonWnd = RadonWindow(self.getImage(), self.image)

    def getImage(self):
        ind = [slice(None), slice(None), slice(None)]
        ind[self.axisorder[-1]] = self.curFrame
        ind = tuple(ind)

        return self.data[ind]

    def setFrame(self, frame: int) -> None:
        self.curFrame = frame
        self.updateImage()

    def saveVal(self):
         for roi in self.ROIs:
            logging.debug(roi.__dict__)
            pos = roi.pos() # grab the outerleft point
            sz = roi.size()
            r = sz[0]/2

            ctrx = pos[0] + r
            ctry = pos[1] + r
            frame = self.curFrame
            axis = self.axis

            angle = roi._lineAngle

            # self.exportVal(axis, frame, ctrx, ctry, angle, pos, sz, r)
            # logging.debug('ctrxctry: {},{}'.format(pos, ctrx, ctry))
            # logging.debug('angle: {}'.format(angle))
            # logging.debug('')
    
    def exportVal(self, axis, frame, ctrx, ctry, angle, pos, sz, r):
        ## set up points based on axis (x, y, z) ##
        if axis == "yz":
            points = [frame, ctrx, ctry]
        elif axis == "xz":
            points = [ctrx, frame, ctry]
        elif axis == "xy":
            points = [ctrx, ctrx, frame]
        
        entry = [points, angle, frame, pos[0], pos[1], sz, r]
        logging.debug(entry)

        ## export to csv ##
        ## csv files will contain [[frame], [x,y,z coordinates], [ROI angle], [circular coordinates]]
        with open('output.csv', mode='w') as output_file:
            wr = csv.writer(output_file, quoting=csv.QUOTE_MINIMAL)
            for row in entry:
                wr.writerow(row)

    def updateImage(self):
        self.image.setImage(self.getImage())
        self.radonWnd.updateImage(self.getImage())

    def doubleClick(self, ev):
        logging.debug('Click: {}, double {}'.format(ev, ev.double()))
        ev.accept()

        r = 32.0
        ctr = np.array(ev.pos())
        pos = ctr - [r, r]

        roi1 = CircleROIwLine(radonWnd=self.radonWnd,
                              pos=pos, size=[2*r, 2*r], removable=True)

        roi1.addToPlot(self.plot)

        roi1.setZValue(10)      # make sure it's above the image
        self.radonWnd.updateROI(roi1)

        self.ROIs.append(roi1)
    
    def roiClicked(self, roi):
        logging.debug("ROI: {}".format(roi))

        self.radonWnd.updateROI(roi)

    def readSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FiberWindow2D")

        self.resize(settings.value("size", QtCore.QSize(800, 600)))
        self.move(settings.value("position", QtCore.QPoint(200, 200)))

        settings.endGroup()

    def writeSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FiberWindow2D")
        settings.setValue("size", self.size())
        settings.setValue("position", self.pos())
        settings.endGroup()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.writeSettings()

        self.radonWnd.close()

        event.accept()


def main():
    logging.basicConfig(level=logging.DEBUG)

    app = QtWidgets.QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)
    
    fw = FiberWindow2d(datafile='/Users/liu81365/Desktop/Research/fibers/fibertracker/Drerio4.h5', axis='yz')
    fw.show()

    return app.exec_()

if __name__ == '__main__':
    sys.exit(main())
