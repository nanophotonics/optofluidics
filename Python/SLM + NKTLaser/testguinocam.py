#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:37:32 2017

@author: Andrei
"""
#python -m PyQt4.uic.pyuic testui.ui -o testui.py
from PyQt4 import QtGui#, QtWidgets # Import the PyQt4 module we'll need
import sys # We need sys so that we can pass argv to QApplication

import testui # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer
import time
from camera import initcam, autoexposure, closecam, takeimagebw, takeimagecolour, exposure
from autofunct import optsetup, compare

import numpy as np
import scipy.misc
from functions import SLMinit, useSLM
from beamshapes import LG, square, gausscomp, nothing, pic, test, angled, displace, rotator
import matplotlib.image as mpimg
intaim=np.load('calib.npy')
defhorpos=5
defhorang=-2.4
defverpos=-25
defverang=8.1
defrotpos=24
defsize=57


class ExampleApp(QtGui.QMainWindow, testui.Ui_MainWindow):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
                            # It sets up layout and widgets that are defined
        
        self.horpos=defhorpos
        self.horv.setText(str(self.horpos))
        self.horang=defhorang
        self.hora.setText(str(self.horang))
        self.verang=defverang
        self.vera.setText(str(self.verang))
        self.verpos=defverpos
        self.verv.setText(str(self.verpos))
        self.rotang=int(defrotpos)
        self.rv.setText(str(self.rotang))        
        self.size=defsize
        self.sv.setText(str(self.size))
        
        self.horposs.valueChanged.connect(self.hp)
        self.horv.editingFinished.connect(self.hpb)    
        self.horangs.valueChanged.connect(self.ha)
        self.hora.editingFinished.connect(self.hab)     
        self.verposs.valueChanged.connect(self.vp)
        self.verv.editingFinished.connect(self.vpb)   
        self.verangs.valueChanged.connect(self.va)
        self.vera.editingFinished.connect(self.vab)      
        self.sizes.valueChanged.connect(self.sizer)
        self.sv.editingFinished.connect(self.sizerb)
        self.dial.valueChanged.connect(self.rot)
        self.rv.editingFinished.connect(self.rotb)
        
        self.hpb()
        self.rotb()
        self.vpb()
        self.vab()
        self.hab()
        self.sizerb()
        self.slmupdate()                     
        
    def hp(self):
        minv=float(self.horposi.text())
        maxv=float(self.horposa.text())
        posn=float(self.horposs.value())/100*(maxv-minv)+minv
        self.horv.editingFinished.disconnect(self.hpb)
        self.horv.setText(str(posn))
        self.horv.editingFinished.connect(self.hpb)
        self.horpos=posn
        self.slmupdate()
        
    def hpb(self):
        minv=float(self.horposi.text())
        maxv=float(self.horposa.text())
        posn=float(self.horv.text())
        self.horposs.valueChanged.disconnect(self.hp)
        self.horposs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horposs.valueChanged.connect(self.hp)
        self.horpos=posn
        self.slmupdate()

    def rot(self):
        posn=float(self.dial.value())
        self.rv.editingFinished.disconnect(self.rotb)
        self.rv.setText(str(posn))
        self.rv.editingFinished.connect(self.rotb)
        self.rotang=posn
        self.slmupdate()
        
    def rotb(self):
        posn=float(self.rv.text())
        self.dial.valueChanged.disconnect(self.rot)
        self.dial.setValue(int(posn))
        self.dial.valueChanged.connect(self.rot)
        self.rotang=posn
        self.slmupdate()

    def ha(self):
        minv=float(self.horangi.text())
        maxv=float(self.horanga.text())
        posn=(float(self.horangs.value())/100*(maxv-minv)+minv)  
        self.hora.editingFinished.disconnect(self.hab)
        self.hora.setText(str(posn))
        self.hora.editingFinished.connect(self.hab)
        self.horang=posn
        self.slmupdate()
        
    def hab(self):
        minv=float(self.horangi.text())
        maxv=float(self.horanga.text())
        posn=float(self.hora.text())
        self.horangs.valueChanged.disconnect(self.ha)
        self.horangs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horangs.valueChanged.connect(self.ha)
        self.horang=posn
        self.slmupdate()
        
    def vp(self):
        minv=float(self.verposi.text())
        maxv=float(self.verposa.text())
        posn=(float(self.verposs.value())/100*(maxv-minv)+minv)
        self.verv.editingFinished.disconnect(self.vpb)
        self.verv.setText(str(posn))
        self.verv.editingFinished.connect(self.vpb)
        self.verpos=posn
        self.slmupdate()
        
    def vpb(self):
        minv=float(self.verposi.text())
        maxv=float(self.verposa.text())
        posn=float(self.verv.text())
        self.verposs.valueChanged.disconnect(self.vp)
        self.verposs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verposs.valueChanged.connect(self.vp)
        self.verpos=posn
        self.slmupdate()        
        
    def va(self):
        minv=float(self.verangi.text())
        maxv=float(self.veranga.text())
        posn=float(self.verangs.value())/100*(maxv-minv)+minv
        self.vera.editingFinished.disconnect(self.vab)
        self.vera.setText(str(posn))
        self.vera.editingFinished.connect(self.vab)
        self.verang=posn
        self.slmupdate()
        
    def vab(self):
        minv=float(self.verangi.text())
        maxv=float(self.veranga.text())
        posn=float(self.vera.text())
        self.verangs.valueChanged.disconnect(self.va)
        self.verangs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verangs.valueChanged.connect(self.va)
        self.verang=posn
        self.slmupdate()
        
    def sizer(self):
        minv=float(self.sizei.text())
        maxv=float(self.sizea.text())
        posn=(float(self.sizes.value())/100*(maxv-minv)+minv)  
        self.sv.editingFinished.disconnect(self.sizerb)
        self.sv.setText(str(posn))
        self.sv.editingFinished.connect(self.sizerb)
        self.size=posn
        self.slmupdate()
        
    def sizerb(self):
        minv=float(self.sizei.text())
        maxv=float(self.sizea.text())
        posn=float(self.sv.text())  
        self.sizes.valueChanged.disconnect(self.sizer)
        self.sizes.setValue(int((posn-minv)/(maxv-minv)*100))
        self.sizes.valueChanged.connect(self.sizer)
        self.size=posn
        self.slmupdate()
    
    def slmupdate(self):
        distribution=np.real((LG(3,2,self.size,slmpix)))
        distribution=displace(angled(rotator(distribution,self.rotang), self.verang, self.horang, slmpix), self.verpos, self.horpos, slmpix)*intaim
        tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
        scipy.misc.imsave('amp.png', np.abs(distribution))
        scipy.misc.imsave('phase.png', np.angle(distribution))
        tosend=tosend.astype('uint16')
        scipy.misc.imsave('hologram.png', tosend[:,:,0]+np.left_shift(tosend[:,:,1],8))
       #def onpress(self):
     #   print('Pressed!')
        
    def closeEvent(self, event):
        # here you can terminate your threads and do other stuff
        closecam(cam)
        # and afterwards call the closeEvent of the super-class
        super(QtGui.QMainWindow, self).closeEvent(event)


                        # and execute the app

#app=QtWidgets.QApplication.instance()  # checks if QApplication already exists
#if app:    # create QApplication if it doesnt exist
#       app.quit()

def main():
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication
    form = ExampleApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app

slmpix=[512,512]
spfrq=0.15
wavelength=675
slm, slmstate =SLMinit(1)

if __name__ == '__main__':              # if we're running file directly and not importing it
    main()                              # run the main function

