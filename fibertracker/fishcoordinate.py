from ast import And
import os, sys
from turtle import position
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore    import *
from PyQt5.QtGui     import *
from PyQt5.QtWidgets import *
import numpy as np
import logging
import h5py
import skimage.transform
from copy import copy
import pyqtgraph as pg
from scipy import interpolate
from matplotlib import pyplot as plt

SETTINGS_FILE = "fibertracker.ini"

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class ImageCoord(pg.ROI):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        pos = self.pos()
        size = self.size()
        sizex = size[0]
        sizey = size[1]
        ctrx = pos[0]
        ctry = pos[1]
        print("pos: ", pos)

        self.line = pg.PlotDataItem(x=[ctrx], y=[ctry], pen='g')
        self.removeHandle

    def addToPlot(self, plot):
        plot.addItem(self)
        plot.addItem(self.line)

class ImageClick(pg.ImageItem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    sigDoubleClick = QtCore.pyqtSignal(object)

    def hoverEvent(self, ev):
        ev.acceptClicks(button=0)

    def mouseClickEvent(self, ev):
        if ev.double():
            self.sigDoubleClick.emit(ev)

class SubWindow(QtWidgets.QWidget):
    def __init__(self, datafile, coordArr, M, N, axis='yz'):
        super().__init__()

        self.datafile = datafile
        self.axis = axis

        self.coordArr = coordArr
        self.M = M
        self.N = N
        
        self.planeCoords = []

        self.h5file = h5py.File(datafile, 'r')
        self.data = self.h5file['image']

        # setting up image and axis
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
        self.axisunorder = np.argsort(self.axisorder)
        self.shape = np.array(self.shape0)[self.axisorder][:2]
        self.len = self.shape0[self.axisorder[2]]

        self.curPlane = 0

        self.getArbitraryPlanes()

        self.initUI()
        self.readSettings()
        self.updateImage()
    
    def initUI(self):
        self.setWindowTitle("1: Fish Coordiate System")

        self.scrollBar = QtWidgets.QScrollBar(orientation=QtCore.Qt.Horizontal)
        self.scrollBar.setRange(0, self.coordArr.shape[0]-1)
        self.scrollBar.setValue(self.curPlane)
        self.scrollBar.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                                           QtWidgets.QSizePolicy.Fixed))

        self.frameNumberBox = QtWidgets.QSpinBox()
        self.frameNumberBox.setRange(0, self.coordArr.shape[0]-1)

        self.frameNumberBox.setValue(self.curPlane)
        lab1 = QtWidgets.QLabel("&Plane number:")
        lab1.setBuddy(self.frameNumberBox)
        lab1.setAlignment(QtCore.Qt.AlignRight)

        hlayout1 = QtWidgets.QHBoxLayout()
        hlayout1.addWidget(self.scrollBar)
        hlayout1.addWidget(lab1)
        hlayout1.addWidget(self.frameNumberBox)

        self.plot = pg.PlotWidget()
        self.image = ImageClick(image=self.vals[0,:,:])
        self.image.sigDoubleClick.connect(self.doubleClick) 
        self.plot.addItem(self.image)

        self.plot.setAspectLocked()
        rect = self.plot.viewRect()
        self.plot.setRange(rect=rect, padding=None)

        vlayout1 = QtWidgets.QVBoxLayout()
        vlayout1.addWidget(self.plot)
        vlayout1.addLayout(hlayout1)
        self.setLayout(vlayout1)

        # slots
        self.frameNumberBox.valueChanged.connect(self.scrollBar.setValue)
        self.scrollBar.valueChanged.connect(self.frameNumberBox.setValue)
        self.scrollBar.valueChanged.connect(self.setPlane)
    
    def readSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("Window2D")

        self.resize(settings.value("size", QtCore.QSize(800, 600)))
        self.move(settings.value("position", QtCore.QPoint(200, 200)))

        settings.endGroup()
    
    def writeSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FishCoord")
        settings.setValue("size", self.size())
        settings.setValue("position", self.pos())
        settings.endGroup()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.writeSettings()
        # self.sw.close()
        event.accept()

    def getArbitraryPlanes(self, range=[-50, 50]):
        i1 = np.arange(range[0], range[1])
        j1 = np.arange(range[0], range[1])
        [i,j] = np.meshgrid(i1, j1)

        # obtain plane coordinates
        for k, arr in enumerate(self.coordArr):
            ctr = arr
            m = self.M[k]
            n = self.N[k]
            
            self.planeCoords.append(m[np.newaxis, np.newaxis, :] * i[:, :, np.newaxis] + \
                n[np.newaxis, np.newaxis, :] * j[:, :, np.newaxis] + \
                ctr[np.newaxis, np.newaxis, :])
        
        self.planeCoords = np.round(self.planeCoords).astype(int)

        # pull intensity values to display image in transverse plane
        vals = np.zeros(self.planeCoords.shape[:-1])

        for pnum, plane in enumerate(self.planeCoords):
            for i0, row in enumerate(plane):
                for j0, pt in enumerate(row):
                    if (np.all(pt >= 0) and np.all(pt < self.data.shape)):
                        # pull out the value and put it in the right place in the vals matrix
                        vals[pnum, i0, j0] = self.data[pt[0], pt[1], pt[2]]
        print('done')
        ## ERROR HERE. IndexError: tuple index out of range ##
        self.vals = vals

    def setPlane(self, pnum: int) -> None:
        self.curPlane = pnum
        self.updateImage(self.curPlane)
    
    def updateImage(self, pnum = 0):
        self.image.setImage(self.vals[pnum, :, :])
        # self.image.setImage(self.getImage())
        # self.sw.updateImage(self.getImage())

    def doubleClick(self, ev):
        pos = np.array(ev.pos())

        coord = ImageCoord(pos=pos, size=(5.000000, 5.000000), angle=45.0, invertible=True, 
                removable=True, rotatable=False, resizable=False, movable=False)
        
        coord.addToPlot(self.plot)
        self.saveCoord(pos)

class MainWindow(QtWidgets.QWidget):
    def __init__(self, datafile, axis='yz'):
        super().__init__()

        self.datafile = datafile
        self.axis = axis

        self.h5file = h5py.File(datafile, 'r')
        self.data = self.h5file['image']

        # setting up image and axis
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
        self.axisunorder = np.argsort(self.axisorder)
        self.shape = np.array(self.shape0)[self.axisorder][:2]
        self.len = self.shape0[self.axisorder[2]]

        self.curFrame = 0
        self.coordArr = []
        self.M = []
        self.N = []

        self.initUI()
        self.readSettings()
        self.updateImage()
    
    def initUI(self):
        self.setWindowTitle("0: Mark all vertebrae")

        self.scrollBar = QtWidgets.QScrollBar(orientation=QtCore.Qt.Horizontal)
        self.scrollBar.setRange(0, self.len-1)
        self.scrollBar.setValue(self.curFrame)
        self.scrollBar.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                                           QtWidgets.QSizePolicy.Fixed))

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

        self.interpolateButton = QtWidgets.QPushButton("Interpolate")
        hlayout1.addWidget(self.interpolateButton)
        self.interpolateButton.clicked.connect(self.interpolate)
        self.interpolateButton.clicked.connect(self.showSW)

        # slots
        self.frameNumberBox.valueChanged.connect(self.scrollBar.setValue)
        self.scrollBar.valueChanged.connect(self.frameNumberBox.setValue)
        self.scrollBar.valueChanged.connect(self.setFrame)
    
    def readSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("Window2D")

        self.resize(settings.value("size", QtCore.QSize(800, 600)))
        self.move(settings.value("position", QtCore.QPoint(200, 200)))

        settings.endGroup()
    
    def writeSettings(self):
        settings = QtCore.QSettings(SETTINGS_FILE, QtCore.QSettings.IniFormat)

        settings.beginGroup("FishCoord")
        settings.setValue("size", self.size())
        settings.setValue("position", self.pos())
        settings.endGroup()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.writeSettings()
        self.sw.close()
        event.accept()

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
        pos = np.array(ev.pos())

        coord = ImageCoord(pos=pos, size=(5.000000, 5.000000), angle=45.0, invertible=True, 
                removable=True, rotatable=False, resizable=False, movable=False)
        
        coord.addToPlot(self.plot)
        self.saveCoord(pos)

    def saveCoord(self, pos):
        if self.axis.lower() == 'yz':
            self.coordArr.append([self.curFrame, pos[0], pos[1]])
        elif self.axis.lower() == 'xz':
            self.coordArr.append([pos[0], self.curFrame, pos[1]])
        elif self.axis.lower() == 'xy':
            self.coordArr.append([pos[0], pos[1], self.curFrame])
    
    def interpolate(self):
        # prepare coordinate array and sort
        idx = np.argsort(self.sortArr())
        self.coordArr = np.array(self.coordArr)[idx.astype(int)]

        x_arr = []
        y_arr = []
        z_arr = []
        
        for point in self.coordArr:
            x_arr.append(point[0])
            y_arr.append(point[1])
            z_arr.append(point[2])
        
        x_arr = np.array(x_arr)
        y_arr = np.array(y_arr)
        z_arr = np.array(z_arr)

        # interpolation with univariate spline to obtain tangent vector
        dx = x_arr[1:] - x_arr[:-1]
        dy = y_arr[1:] - y_arr[:-1]
        dz = z_arr[1:] - z_arr[:-1]
        ds = np.sqrt(dx**2 + dy**2 + dz**2) 

        s = np.insert(np.cumsum(ds), 0, 0) # arc length
        spx = interpolate.UnivariateSpline(s, x_arr)
        spy = interpolate.UnivariateSpline(s, y_arr)
        spz = interpolate.UnivariateSpline(s, z_arr)

        s1 = np.linspace(s[0], s[-1], 20)
        xs = spx(s1) # function x(s) that gives point 
        ys = spy(s1)
        zs = spz(s1)

        dxs = np.array(spx.derivative()(s1))
        dys = np.array(spy.derivative()(s1))
        dzs = np.array(spz.derivative()(s1))

        # generate perpendicular plane with tangent vector
        t = np.array([dxs, dys, dzs]).T
        pplane = []

        for curr in t:
            U = self.generatePPlane(curr)
            pplane.append(U[:, 1:])

            self.M.append(U[:, 1])
            self.N.append(U[:, 2])

        # replace coordArr with spline coordinates for same indexing
        self.coordArr = np.concatenate((np.vstack(xs), np.vstack(ys), np.vstack(zs)), axis = 1)

        # #%%
        # fig = plt.figure()
        # ax = plt.axes(projection = '3d')
        # ax.scatter(x_arr, y_arr, z_arr)
        # ax.plot(xs, ys, zs)
        # ax.set_xlabel('X-axis', fontweight ='bold')
        # ax.set_ylabel('Y-axis', fontweight ='bold')
        # ax.set_zlabel('Z-axis', fontweight ='bold')
        # plt.show()
        # # %%

        self.sw = SubWindow(self.datafile, self.coordArr, self.M, self.N, self.axis)
        self.sw.show()

    
    def generatePPlane(self, t):
        t = t / np.linalg.norm(t)

        # gram-schmidt orthogonalization
        V = np.vstack((t, [0, 0, 1], [0, 1, 0])).T
        U = np.zeros_like(V)

        U[:, 0] = V[:, 0] / np.linalg.norm(V[:, 0])

        for i in range(1,3):
            U[:, i] = V[:, i]
            for j in range(0,i):
                U[:, i] = U[:, i] - np.dot(U[:, i], U[:, j]) * U[:, j]

            U[:, i] = U[:, i] / np.linalg.norm(U[:, i])

        return U
    
    def sortArr(self, metric='euclidean'):
        aa = np.asarray(self.coordArr)
        bb = np.atleast_2d(np.array([0,0,0]))
        a_dim = aa.ndim
        b_dim = bb.ndim

        if a_dim == 1:
            aa = aa.reshape(1, 1, aa.shape[0])
        if a_dim >= 2:
            aa = aa.reshape(np.prod(aa.shape[:-1]), 1, aa.shape[-1])
        if b_dim > 2:
            bb = bb.reshape(np.prod(bb.shape[:-1]), bb.shape[-1])

        diff = aa - bb
        dist_arr = np.einsum('ijk,ijk->ij', diff, diff)
        if metric[:1] == 'e':
            dist_arr = np.sqrt(dist_arr)
        dist_arr = np.squeeze(dist_arr)
        
        return dist_arr

    def showSW(self):
        self.sw = SubWindow(self.datafile, self.coordArr, self.M, self.N, self.axis)
        self.sw.show()

def main():
    logging.basicConfig(level=logging.DEBUG)

    app = QtWidgets.QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)

    df = '/Users/etytel01/Documents/Fibers/Drerio_6.h5'
    ax = 'xy'
    
    fw = MainWindow(datafile=df, axis=ax)
    fw.show()

    return app.exec_()

if __name__ == '__main__':
    sys.exit(main())