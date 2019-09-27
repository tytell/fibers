import os, sys
from PyQt5 import QtGui, QtCore, QtWidgets
import numpy as np
import logging
import h5py
import skimage.transform
from copy import copy
import pyqtgraph as pg

SETTINGS_FILE = "fibertracker.ini"

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class ImageROI(pg.ImageItem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    sigDoubleClick = QtCore.pyqtSignal(object)

    def hoverEvent(self, ev):
        ev.acceptClicks(button=0)

    def mouseClickEvent(self, ev):
        logging.debug("Mousce click: {}, double {}".format(ev, ev.double()))
        if ev.double():
            self.sigDoubleClick.emit(ev)

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

        self.plot = pg.PlotWidget()
        self.image = ImageROI(image=self.getImage())
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

        self.doneButton.clicked.connect(self.close)

        # set up extra window for showing zoomed in ROIs and radon results
        self.zoomWnd = QtWidgets.QWidget()
        self.zoomWnd.setWindowTitle("Radon")

        self.zoomPlot = pg.PlotWidget()
        self.zoomImg = pg.ImageItem()
        self.zoomPlot.addItem(self.zoomImg)

        self.zoomPlot.setAspectLocked()

        self.radonPlot = pg.PlotWidget()
        self.radonImg = pg.ImageItem()
        self.radonPlot.addItem(self.radonImg)

        hlayout2 = QtWidgets.QHBoxLayout()
        hlayout2.addWidget(self.zoomPlot)
        hlayout2.addWidget(self.radonPlot)

        self.zoomWnd.setLayout(hlayout2)
        self.zoomWnd.show()

    def getImage(self):
        ind = [slice(None), slice(None), slice(None)]
        ind[self.axisorder[-1]] = self.curFrame
        ind = tuple(ind)

        return self.data[ind]

    def setFrame(self, frame: int) -> None:
        self.curFrame = frame
        self.updateImage()

    def updateImage(self):
        self.image.setImage(self.getImage())

    def doubleClick(self, ev):
        logging.debug('Click: {}, double {}'.format(ev, ev.double()))
        ev.accept()

        sz = np.array([64, 64])
        ctr = np.array(ev.pos())
        pos = ctr - sz/2

        roi1 = pg.RectROI(pos=pos, size=sz, centered=True,
                          removable=True)

        roi1.addRotateHandle([1, 0], [0.5, 0.5])
        roi1.setAcceptedMouseButtons(QtCore.Qt.LeftButton)

        roi1.sigClicked.connect(self.roiClicked)

        self.plot.addItem(roi1)
        roi1.setZValue(10)      # make sure it's above the image
        self.updateZoomImage(roi1)

        self.ROIs.append(roi1)

    def roiClicked(self, roi):
        logging.debug("ROI: {}".format(roi))

        roi.sigRegionChangeFinished.connect(self.updateZoomImage)

        self.updateZoomImage(roi)

    def updateZoomImage(self, roi):
        roiImg = roi.getArrayRegion(self.getImage(), self.image)

        radius = min(roiImg.shape) // 2
        c0, c1 = np.ogrid[0:roiImg.shape[0], 0:roiImg.shape[1]]
        reconstruction_circle = ((c0 - roiImg.shape[0] // 2) ** 2
                                 + (c1 - roiImg.shape[1] // 2) ** 2)
        outside = reconstruction_circle >= radius ** 2
        inside = np.logical_not(outside)

        roiImg -= np.mean(roiImg[inside])
        roiImg[outside] = 0.0

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

        self.radonPlot.plot(x=theta, y=angleintensity * sinogram.shape[1], z=10, pen='r')

        angle = theta[np.argmax(angleintensity)]

        self.radonPlot.addItem(pg.InfiniteLine(pos=angle, angle=90, pen='g'))

        ctr = QtCore.QPointF(roiImg.shape[0]/2.0, roiImg.shape[1]/2.0)
        self.zoomPlot.addItem(pg.InfiniteLine(pos=ctr, angle=angle, pen='g'))

        self.zoomWnd.show()

    def readSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FiberWindow2D")

        self.resize(settings.value("size", QtCore.QSize(800, 600)))
        self.move(settings.value("position", QtCore.QPoint(200, 200)))

        settings.endGroup()

        settings.beginGroup("ZoomWindow")

        self.zoomWnd.resize(settings.value("size", QtCore.QSize(200, 200)))
        self.zoomWnd.move(settings.value("position", QtCore.QPoint(1000, 200)))

        settings.endGroup()

    def writeSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FiberWindow2D")
        settings.setValue("size", self.size())
        settings.setValue("position", self.pos())
        settings.endGroup()

        settings.beginGroup("ZoomWindow")
        settings.setValue("size", self.zoomWnd.size())
        settings.setValue("position", self.zoomWnd.pos())
        settings.endGroup()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.writeSettings()

        self.zoomWnd.close()

        event.accept()


def main():
    logging.basicConfig(level=logging.DEBUG)

    app = QtWidgets.QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)

    fw = FiberWindow2d(datafile='/Users/etytel01/Documents/Fibers/Drerio4.h5', axis='yz')
    fw.show()

    return app.exec_()

if __name__ == '__main__':
    sys.exit(main())