#!/usr/bin/python

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from PyQt5.QtWidgets import QSlider
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QGridLayout
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import random

def get_interface(plot_figure):
    """Return an object that plots stuff"""
    class Window(QDialog):
        def __init__(self, parent=None):
            super(Window, self).__init__(parent)
    
            # a figure instance to plot on
            self.figure = plt.figure()
    
            # this is the Canvas Widget that displays the `figure`
            # it takes the `figure` instance as a parameter to __init__
            self.canvas = FigureCanvas(self.figure)
    
            # this is the Navigation widget
    
            # Just some button connected to `plot` method
            self.button = QPushButton('Plot')
#            self.button.clicked.connect(self.plot)
    
            # set the layout
#            layout = QVBoxLayout()
            layout = QGridLayout()
#            layout.addWidget(self.button)
            self.row = 1
            self.layout = layout
            self.setLayout(layout)
        def plot(self):
            '''Plot the figure'''
            # instead of ax.hold(False)
            import matplotlib.pyplot as plt
            plt.close() # close figure
            self.layout.removeWidget(self.canvas)
            self.figure = plot_figure(self) # get that figure
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)
            # refresh canvas
            self.canvas = FigureCanvas(self.figure)
            # add the new canvas at the position of the old one
#            self.layout.addWidget(self.toolbar,0,0,1,0)
            self.layout.addWidget(self.canvas, 1,0,1,0)
            self.setLayout(self.layout)
            self.canvas.draw()
        def add_combobox(self,cs,name="",label=None):
            """Add a combo box"""
            if label is None: label=name
            combo = QtWidgets.QComboBox(objectName=label)
            combo.addItems(cs)
            combo.currentTextChanged.connect(self.plot)
            self.row += 1 # increase
            if name!="": # empty string
              lb = QtWidgets.QLabel(name)
              self.layout.addWidget(lb,self.row,0)
              self.layout.addWidget(combo,self.row,1)
            else:
              self.layout.addWidget(combo,self.row,1,1,0)
        def get_combobox(self,name):
            """Get the value of a combobox"""
            obj = self.findChild(QtWidgets.QComboBox,name)
            return obj.currentText()
        def add_slider(self,name,label=None,vs=range(100),v0=0):
            """Add a slider"""
            vmin = 0
            vmax = len(vs)-1
            if label is None: label = name
            slider = QSlider(Qt.Horizontal,objectName=label)
            slider.vs = vs # store
            slider.setMinimum(vmin)
            slider.setMaximum(vmax)
            slider.setTickPosition(QSlider.TicksBelow)
            slider.setTickInterval(1)
            if v0 is not None: slider.setValue(v0)
            else: slider.setValue(vmin)
#            slider.sliderMoved.connect(self.plot)
            slider.valueChanged.connect(self.plot)
            label = QtWidgets.QLabel(name)
            self.row += 1 # increase counter
            self.layout.addWidget(label,self.row,0)
            self.layout.addWidget(slider,self.row,1)
            self.setLayout(self.layout)
        def get_slider(self,name):
            """Get the value of a slider"""
            slider = self.findChild(QSlider,name)
            out = slider.value()
            return slider.vs[int(out)]

    app = QApplication(sys.argv)
    main = Window()
    return app,main

if __name__ == '__main__':
    def funfig(obj): # dummy function
        xs = np.linspace(0.,10.,300) # xgrid
        # get the value of the slider
        ys = np.cos(xs*obj.get_slider("k") + obj.get_slider("phi")) 
        fig = plt.figure(0) # initialize figure
        plt.plot(xs,ys,c=obj.get_combobox("c")) # plot data
        return fig # return figure
    app,main = get_interface(funfig) # get the interface
    ks = np.linspace(1.0,3.0,50) # wavevectors
    ps = np.linspace(0.0,2.0,50)*np.pi # phases
    main.add_slider("Wavevector",label="k",vs=ks) # initialize the slider
    main.add_slider("Phase",label="phi",vs=ps) # initialize the slider
    # initialize the combobox
    main.add_combobox(["red","blue","black"],name="Color",label="c") 
    main.plot()
    main.show()
    sys.exit(app.exec_())
