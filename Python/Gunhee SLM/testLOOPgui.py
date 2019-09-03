#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:37:32 2017

@author: Andrei
"""
#GUI for SLM control with three different fundamental modes loop
#python -m PyQt4.uic.pyuic testui.ui -o testui.py
from PyQt4 import QtGui, QtCore#, QtWidgets # Import the PyQt4 module we'll need
from PyQt4.QtCore import pyqtSignal
import sys # We need sys so that we can pass argv to QApplication

import testLOOPui # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer
import time
from cameraclass import camera
from autofunct import optsetup, compare, cropcurrent
from matplotlib import cm
import numpy as np
import scipy.misc
from slmclass import SLM
from beamclass2 import beamshapes
import matplotlib.image as mpimg
import threading
import datetime

intaim=np.load('calib.npy')
wavelength=800
slmpix=[512,512]
spfrq=-0.15
slm=SLM(slmpix, wavelength, spfrq, 1)
beamgen1=beamshapes(slmpix)
beamgen2=beamshapes(slmpix)
beamgen3=beamshapes(slmpix)

defhorpos=0
defhorang=0
defverpos=0
defverang=0
defrotpos=0
deffocus=0
defsize=75
defl=0
defg=0

defhorpos_b=0
defhorang_b=0
defverpos_b=0
defverang_b=0
defrotpos_b=0
deffocus_b=0
defsize_b=75
defl_b=1
defg_b=0

defhorpos_c=0
defhorang_c=0
defverpos_c=0
defverang_c=0
defrotpos_c=0
deffocus_c=0
defsize_c=75
defl_c=3
defg_c=0

deftm = 1.0
cam = camera(10)

#thread to save images
#thread was used to control the time flow consistently
def imagethread(num,expo,count,tm):
    now = datetime.datetime.now()
    cam.exposure(expo)
    if num==1:
        scipy.misc.imsave('1stmode_{}.png'.format(count), cam.takeimagebw())
        print('1st mode')
        F = open('1stmode_{}.txt'.format(count),'w')
        F.write('The image was taken at '+str(now)+'\n')
        F.write('1st mode '+str(now-datetime.timedelta(seconds=tm/2))+' ~ '+str(now+datetime.timedelta(seconds=tm/2)))
        F.close()
    if num==2:
        scipy.misc.imsave('2ndmode_{}.png'.format(count), cam.takeimagebw())
        print('2nd mode')
        F = open('2ndmode_{}.txt'.format(count),'w')
        F.write('The image was taken at '+str(now)+'\n')
        F.write('2nd mode '+str(now-datetime.timedelta(seconds=tm/2))+' ~ '+str(now+datetime.timedelta(seconds=tm/2)))
        F.close()
    if num==3:
        scipy.misc.imsave('3rdmode_{}.png'.format(count), cam.takeimagebw())
        print('3rd mode')
        F = open('3rdmode_{}.txt'.format(count),'w')
        F.write('The image was taken at '+str(now)+'\n')
        F.write('3rd mode '+str(now-datetime.timedelta(seconds=tm/2))+' ~ '+str(now+datetime.timedelta(seconds=tm/2)))
        F.close()

#thread to loop three different SLM modes
#thread was used to avoid the loop to be locked in main window
class QThread1(QtCore.QThread):
    sources = {}
    def __init__(self, parent=None):
        QtCore.QThread.__init__(self, parent)
    def on_source(self, dict):
        for key, value in dict.iteritems():
            self.sources[key]=value
    def on_time(self, time):
        self.tm = time
    def run(self):
        self.running = True
        count = 0
        slm.useSLM(self.sources['dist1'])
        expo1=cam.autoexposure()
        slm.useSLM(self.sources['dist2'])
        expo2=cam.autoexposure()
        slm.useSLM(self.sources['dist3'])
        expo3=cam.autoexposure()
        while self.running:
            thrd1 = threading.Thread(target=imagethread,args=(1,expo1,count,self.tm))
            thrd2 = threading.Thread(target=imagethread,args=(2,expo2,count,self.tm))
            thrd3 = threading.Thread(target=imagethread,args=(3,expo3,count,self.tm))
            slm.useSLM(self.sources['dist1'])
            time.sleep(self.tm/2)
            thrd1.start()
            time.sleep(self.tm/2)
            slm.useSLM(self.sources['dist2'])
            time.sleep(self.tm/2)
            thrd2.start()
            time.sleep(self.tm/2)
            slm.useSLM(self.sources['dist3'])
            time.sleep(self.tm/2)
            thrd3.start()
            time.sleep(self.tm/2)
            count = count+1


class ExampleApp(QtGui.QMainWindow, testLOOPui.Ui_MainWindow):
    sig = pyqtSignal(dict)
    sigtm = pyqtSignal(float)
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
    # It sets up layout and widgets that are defined
        self.l=defl
        self.lbox.setValue(self.l)
        self.g=defg
        self.gbox.setValue(self.g)
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
        self.focus = deffocus
        self.fv.setText(str(self.focus))
        self.size=defsize
        self.sv.setText(str(self.size))
        self.tm = deftm
        self.timev.setText(str(self.tm))

        self.horposs.valueChanged.connect(self.hp)
        self.horv.editingFinished.connect(self.hpb)
        self.horangs.valueChanged.connect(self.ha)
        self.hora.editingFinished.connect(self.hab)
        self.verposs.valueChanged.connect(self.vp)
        self.verv.editingFinished.connect(self.vpb)
        self.verangs.valueChanged.connect(self.va)
        self.vera.editingFinished.connect(self.vab)
        self.focuss.valueChanged.connect(self.fo)
        self.fv.editingFinished.connect(self.fob)
        self.sizes.valueChanged.connect(self.sizer)
        self.sv.editingFinished.connect(self.sizerb)
        self.dial.valueChanged.connect(self.rot)
        self.rv.editingFinished.connect(self.rotb)
        self.lbox.valueChanged.connect(self.modechange)
        self.gbox.valueChanged.connect(self.modechange)


        self.l_b=defl_b
        self.lbox_b.setValue(self.l_b)
        self.g_b=defg_b
        self.gbox_b.setValue(self.g_b)
        self.horpos_b=defhorpos_b
        self.horv_b.setText(str(self.horpos_b))
        self.horang_b=defhorang_b
        self.hora_b.setText(str(self.horang_b))
        self.verang_b=defverang_b
        self.vera_b.setText(str(self.verang_b))
        self.verpos_b=defverpos_b
        self.verv_b.setText(str(self.verpos_b))
        self.rotang_b=int(defrotpos_b)
        self.rv_b.setText(str(self.rotang_b))
        self.focus_b = deffocus_b
        self.fv_b.setText(str(self.focus_b))
        self.size_b=defsize_b
        self.sv_b.setText(str(self.size_b))

        self.horpos_b_s.valueChanged.connect(self.hp_b)
        self.horv_b.editingFinished.connect(self.hpb_b)
        self.horang_b_s.valueChanged.connect(self.ha_b)
        self.hora_b.editingFinished.connect(self.hab_b)
        self.verpos_b_s.valueChanged.connect(self.vp_b)
        self.verv_b.editingFinished.connect(self.vpb_b)
        self.verang_b_s.valueChanged.connect(self.va_b)
        self.vera_b.editingFinished.connect(self.vab_b)
        self.focus_b_s.valueChanged.connect(self.fo_b)
        self.fv_b.editingFinished.connect(self.fob_b)
        self.sizes_b.valueChanged.connect(self.sizer_b)
        self.sv_b.editingFinished.connect(self.sizerb_b)
        self.dial_b.valueChanged.connect(self.rot_b)
        self.rv_b.editingFinished.connect(self.rotb_b)
        self.lbox_b.valueChanged.connect(self.modechange_b)
        self.gbox_b.valueChanged.connect(self.modechange_b)


        self.l_c=defl_c
        self.lbox_c.setValue(self.l_c)
        self.g_c=defg_c
        self.gbox_c.setValue(self.g_c)
        self.horpos_c=defhorpos_c
        self.horv_c.setText(str(self.horpos_c))
        self.horang_c=defhorang_c
        self.hora_c.setText(str(self.horang_c))
        self.verang_c=defverang_c
        self.vera_c.setText(str(self.verang_c))
        self.verpos_c=defverpos_c
        self.verv_c.setText(str(self.verpos_c))
        self.rotang_c=int(defrotpos_c)
        self.rv_c.setText(str(self.rotang_c))
        self.focus_c = deffocus_c
        self.fv_c.setText(str(self.focus_c))
        self.size_c=defsize_c
        self.sv_c.setText(str(self.size_c))

        self.horpos_c_s.valueChanged.connect(self.hp_c)
        self.horv_c.editingFinished.connect(self.hpb_c)
        self.horang_c_s.valueChanged.connect(self.ha_c)
        self.hora_c.editingFinished.connect(self.hab_c)
        self.verpos_c_s.valueChanged.connect(self.vp_c)
        self.verv_c.editingFinished.connect(self.vpb_c)
        self.verang_c_s.valueChanged.connect(self.va_c)
        self.vera_c.editingFinished.connect(self.vab_c)
        self.focus_c_s.valueChanged.connect(self.fo_c)
        self.fv_c.editingFinished.connect(self.fob_c)
        self.sizes_c.valueChanged.connect(self.sizer_c)
        self.sv_c.editingFinished.connect(self.sizerb_c)
        self.dial_c.valueChanged.connect(self.rot_c)
        self.rv_c.editingFinished.connect(self.rotb_c)
        self.lbox_c.valueChanged.connect(self.modechange_c)
        self.gbox_c.valueChanged.connect(self.modechange_c)

        self.timev.editingFinished.connect(self.timevchange)

        self.hpb()
        self.rotb()
        self.vpb()
        self.vab()
        self.fob()
        self.hab()
        self.sizerb()
        self.hpb_b()
        self.rotb_b()
        self.vpb_b()
        self.vab_b()
        self.fob_b()
        self.hab_b()
        self.sizerb_b()
        self.hpb_c()
        self.rotb_c()
        self.vpb_c()
        self.vab_c()
        self.fob_c()
        self.hab_c()
        self.sizerb_c()
        self.updateb.clicked.connect(self.paramupdate_a)
        self.updateb_b.clicked.connect(self.paramupdate_b)
        self.updateb_c.clicked.connect(self.paramupdate_c)
        self.sendb.clicked.connect(self.on_send)
        self.pauseb.clicked.connect(self.on_pause)
        #self.slmupdate()

    def modechange(self):
        self.l=self.lbox.value()
        self.g=self.gbox.value()
        #self.slmupdate()


    def hp(self):
        minv=float(self.horposi.text())
        maxv=float(self.horposa.text())
        posn=float(self.horposs.value())/100*(maxv-minv)+minv
        self.horv.editingFinished.disconnect(self.hpb)
        self.horv.setText(str(posn))
        self.horv.editingFinished.connect(self.hpb)
        self.horpos=posn
        #self.slmupdate()

    def hpb(self):
        minv=float(self.horposi.text())
        maxv=float(self.horposa.text())
        posn=float(self.horv.text())
        self.horposs.valueChanged.disconnect(self.hp)
        self.horposs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horposs.valueChanged.connect(self.hp)
        self.horpos=posn
        #self.slmupdate()

    def rot(self):
        posn=float(self.dial.value())
        self.rv.editingFinished.disconnect(self.rotb)
        self.rv.setText(str(posn))
        self.rv.editingFinished.connect(self.rotb)
        self.rotang=posn
        #self.slmupdate()

    def rotb(self):
        posn=float(self.rv.text())
        self.dial.valueChanged.disconnect(self.rot)
        self.dial.setValue(int(posn))
        self.dial.valueChanged.connect(self.rot)
        self.rotang=posn
        #self.slmupdate()

    def ha(self):
        minv=float(self.horangi.text())
        maxv=float(self.horanga.text())
        posn=(float(self.horangs.value())/100*(maxv-minv)+minv)
        self.hora.editingFinished.disconnect(self.hab)
        self.hora.setText(str(posn))
        self.hora.editingFinished.connect(self.hab)
        self.horang=posn
        #self.slmupdate()

    def hab(self):
        minv=float(self.horangi.text())
        maxv=float(self.horanga.text())
        posn=float(self.hora.text())
        self.horangs.valueChanged.disconnect(self.ha)
        self.horangs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horangs.valueChanged.connect(self.ha)
        self.horang=posn
        #self.slmupdate()

    def vp(self):
        minv=float(self.verposi.text())
        maxv=float(self.verposa.text())
        posn=(float(self.verposs.value())/100*(maxv-minv)+minv)
        self.verv.editingFinished.disconnect(self.vpb)
        self.verv.setText(str(posn))
        self.verv.editingFinished.connect(self.vpb)
        self.verpos=posn
        #self.slmupdate()

    def vpb(self):
        minv=float(self.verposi.text())
        maxv=float(self.verposa.text())
        posn=float(self.verv.text())
        self.verposs.valueChanged.disconnect(self.vp)
        self.verposs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verposs.valueChanged.connect(self.vp)
        self.verpos=posn
        #self.slmupdate()

    def va(self):
        minv=float(self.verangi.text())
        maxv=float(self.veranga.text())
        posn=float(self.verangs.value())/100*(maxv-minv)+minv
        self.vera.editingFinished.disconnect(self.vab)
        self.vera.setText(str(posn))
        self.vera.editingFinished.connect(self.vab)
        self.verang=posn
        #self.slmupdate()

    def vab(self):
        minv=float(self.verangi.text())
        maxv=float(self.veranga.text())
        posn=float(self.vera.text())
        self.verangs.valueChanged.disconnect(self.va)
        self.verangs.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verangs.valueChanged.connect(self.va)
        self.verang=posn
        #self.slmupdate()

    def fo(self):
        minv=float(self.focusi.text())
        maxv=float(self.focusa.text())
        posn=float(self.focuss.value())/100*(maxv-minv)+minv
        self.fv.editingFinished.disconnect(self.fob)
        self.fv.setText(str(posn))
        self.fv.editingFinished.connect(self.fob)
        self.focus=posn
        #self.slmupdate()

    def fob(self):
        minv=float(self.focusi.text())
        maxv=float(self.focusa.text())
        posn=float(self.fv.text())
        self.focuss.valueChanged.disconnect(self.fo)
        self.focuss.setValue(int((posn-minv)/(maxv-minv)*100))
        self.focuss.valueChanged.connect(self.fo)
        self.focus=posn
        #self.slmupdate()

    def sizer(self):
        minv=float(self.sizei.text())
        maxv=float(self.sizea.text())
        posn=(float(self.sizes.value())/100*(maxv-minv)+minv)
        self.sv.editingFinished.disconnect(self.sizerb)
        self.sv.setText(str(posn))
        self.sv.editingFinished.connect(self.sizerb)
        self.size=posn
        #self.slmupdate()

    def sizerb(self):
        minv=float(self.sizei.text())
        maxv=float(self.sizea.text())
        posn=float(self.sv.text())
        self.sizes.valueChanged.disconnect(self.sizer)
        self.sizes.setValue(int((posn-minv)/(maxv-minv)*100))
        self.sizes.valueChanged.connect(self.sizer)
        self.size=posn
        #self.slmupdate()

    def modechange_b(self):
        self.l_b=self.lbox_b.value()
        self.g_b=self.gbox_b.value()
        #self.slmupdate()

    def hp_b(self):
        minv=float(self.horpos_b_i.text())
        maxv=float(self.horpos_b_a.text())
        posn=float(self.horpos_b_s.value())/100*(maxv-minv)+minv
        self.horv_b.editingFinished.disconnect(self.hpb_b)
        self.horv_b.setText(str(posn))
        self.horv_b.editingFinished.connect(self.hpb_b)
        self.horpos_b=posn
        #self.slmupdate()

    def hpb_b(self):
        minv=float(self.horpos_b_i.text())
        maxv=float(self.horpos_b_a.text())
        posn=float(self.horv_b.text())
        self.horpos_b_s.valueChanged.disconnect(self.hp_b)
        self.horpos_b_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horpos_b_s.valueChanged.connect(self.hp_b)
        self.horpos_b=posn
        #self.slmupdate()

    def rot_b(self):
        posn=float(self.dial_b.value())
        self.rv_b.editingFinished.disconnect(self.rotb_b)
        self.rv_b.setText(str(posn))
        self.rv_b.editingFinished.connect(self.rotb_b)
        self.rotang_b=posn
        #self.slmupdate()

    def rotb_b(self):
        posn=float(self.rv_b.text())
        self.dial_b.valueChanged.disconnect(self.rot_b)
        self.dial_b.setValue(int(posn))
        self.dial_b.valueChanged.connect(self.rot_b)
        self.rotang_b=posn
        #self.slmupdate()

    def ha_b(self):
        minv=float(self.horang_b_i.text())
        maxv=float(self.horang_b_a.text())
        posn=(float(self.horang_b_s.value())/100*(maxv-minv)+minv)
        self.hora_b.editingFinished.disconnect(self.hab_b)
        self.hora_b.setText(str(posn))
        self.hora_b.editingFinished.connect(self.hab_b)
        self.horang_b=posn
        #self.slmupdate()

    def hab_b(self):
        minv=float(self.horang_b_i.text())
        maxv=float(self.horang_b_a.text())
        posn=float(self.hora_b.text())
        self.horang_b_s.valueChanged.disconnect(self.ha_b)
        self.horang_b_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horang_b_s.valueChanged.connect(self.ha_b)
        self.horang_b=posn
        #self.slmupdate()

    def vp_b(self):
        minv=float(self.verpos_b_i.text())
        maxv=float(self.verpos_b_a.text())
        posn=(float(self.verpos_b_s.value())/100*(maxv-minv)+minv)
        self.verv_b.editingFinished.disconnect(self.vpb_b)
        self.verv_b.setText(str(posn))
        self.verv_b.editingFinished.connect(self.vpb_b)
        self.verpos_b=posn
        #self.slmupdate()

    def vpb_b(self):
        minv=float(self.verpos_b_i.text())
        maxv=float(self.verpos_b_a.text())
        posn=float(self.verv_b.text())
        self.verpos_b_s.valueChanged.disconnect(self.vp_b)
        self.verpos_b_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verpos_b_s.valueChanged.connect(self.vp_b)
        self.verpos_b=posn
        #self.slmupdate()

    def va_b(self):
        minv=float(self.verang_b_i.text())
        maxv=float(self.verang_b_a.text())
        posn=float(self.verang_b_s.value())/100*(maxv-minv)+minv
        self.vera_b.editingFinished.disconnect(self.vab_b)
        self.vera_b.setText(str(posn))
        self.vera_b.editingFinished.connect(self.vab_b)
        self.verang_b=posn
        #self.slmupdate()

    def vab_b(self):
        minv=float(self.verang_b_i.text())
        maxv=float(self.verang_b_a.text())
        posn=float(self.vera_b.text())
        self.verang_b_s.valueChanged.disconnect(self.va_b)
        self.verang_b_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verang_b_s.valueChanged.connect(self.va_b)
        self.verang_b=posn
        #self.slmupdate()

    def fo_b(self):
        minv=float(self.focus_b_i.text())
        maxv=float(self.focus_b_a.text())
        posn=float(self.focus_b_s.value())/100*(maxv-minv)+minv
        self.fv_b.editingFinished.disconnect(self.fob_b)
        self.fv_b.setText(str(posn))
        self.fv_b.editingFinished.connect(self.fob_b)
        self.focus_b=posn
        #self.slmupdate()

    def fob_b(self):
        minv=float(self.focus_b_i.text())
        maxv=float(self.focus_b_a.text())
        posn=float(self.fv_b.text())
        self.focus_b_s.valueChanged.disconnect(self.fo_b)
        self.focus_b_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.focus_b_s.valueChanged.connect(self.fo_b)
        self.focus_b=posn
        #self.slmupdate()

    def sizer_b(self):
        minv=float(self.sizei_b.text())
        maxv=float(self.sizea_b.text())
        posn=(float(self.sizes_b.value())/100*(maxv-minv)+minv)
        self.sv_b.editingFinished.disconnect(self.sizerb_b)
        self.sv_b.setText(str(posn))
        self.sv_b.editingFinished.connect(self.sizerb_b)
        self.size_b=posn
        #self.slmupdate()

    def sizerb_b(self):
        minv=float(self.sizei_b.text())
        maxv=float(self.sizea_b.text())
        posn=float(self.sv_b.text())
        self.sizes_b.valueChanged.disconnect(self.sizer_b)
        self.sizes_b.setValue(int((posn-minv)/(maxv-minv)*100))
        self.sizes_b.valueChanged.connect(self.sizer_b)
        self.size_b=posn
        #self.slmupdate()

    def modechange_c(self):
        self.l_c=self.lbox_c.value()
        self.g_c=self.gbox_c.value()
        #self.slmupdate()


    def hp_c(self):
        minv=float(self.horpos_c_i.text())
        maxv=float(self.horpos_c_a.text())
        posn=float(self.horpos_c_s.value())/100*(maxv-minv)+minv
        self.horv_c.editingFinished.disconnect(self.hpb_c)
        self.horv_c.setText(str(posn))
        self.horv_c.editingFinished.connect(self.hpb_c)
        self.horpos_c=posn
        #self.slmupdate()

    def hpb_c(self):
        minv=float(self.horpos_c_i.text())
        maxv=float(self.horpos_c_a.text())
        posn=float(self.horv_c.text())
        self.horpos_c_s.valueChanged.disconnect(self.hp_c)
        self.horpos_c_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horpos_c_s.valueChanged.connect(self.hp_c)
        self.horpos_c=posn
        #self.slmupdate()

    def rot_c(self):
        posn=float(self.dial_c.value())
        self.rv_c.editingFinished.disconnect(self.rotb_c)
        self.rv_c.setText(str(posn))
        self.rv_c.editingFinished.connect(self.rotb_c)
        self.rotang_c=posn
        #self.slmupdate()

    def rotb_c(self):
        posn=float(self.rv_c.text())
        self.dial_c.valueChanged.disconnect(self.rot_c)
        self.dial_c.setValue(int(posn))
        self.dial_c.valueChanged.connect(self.rot_c)
        self.rotang_c=posn
        #self.slmupdate()

    def ha_c(self):
        minv=float(self.horang_c_i.text())
        maxv=float(self.horang_c_a.text())
        posn=(float(self.horang_c_s.value())/100*(maxv-minv)+minv)
        self.hora_c.editingFinished.disconnect(self.hab_c)
        self.hora_c.setText(str(posn))
        self.hora_c.editingFinished.connect(self.hab_c)
        self.horang_c=posn
        #self.slmupdate()

    def hab_c(self):
        minv=float(self.horang_c_i.text())
        maxv=float(self.horang_c_a.text())
        posn=float(self.hora_c.text())
        self.horang_c_s.valueChanged.disconnect(self.ha_c)
        self.horang_c_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.horang_c_s.valueChanged.connect(self.ha_c)
        self.horang_c=posn
        #self.slmupdate()

    def vp_c(self):
        minv=float(self.verpos_c_i.text())
        maxv=float(self.verpos_c_a.text())
        posn=(float(self.verpos_c_s.value())/100*(maxv-minv)+minv)
        self.verv_c.editingFinished.disconnect(self.vpb_c)
        self.verv_c.setText(str(posn))
        self.verv_c.editingFinished.connect(self.vpb_c)
        self.verpos_c=posn
        #self.slmupdate()

    def vpb_c(self):
        minv=float(self.verpos_c_i.text())
        maxv=float(self.verpos_c_a.text())
        posn=float(self.verv_c.text())
        self.verpos_c_s.valueChanged.disconnect(self.vp_c)
        self.verpos_c_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verpos_c_s.valueChanged.connect(self.vp_c)
        self.verpos_c=posn
        #self.slmupdate()

    def va_c(self):
        minv=float(self.verang_c_i.text())
        maxv=float(self.verang_c_a.text())
        posn=float(self.verang_c_s.value())/100*(maxv-minv)+minv
        self.vera_c.editingFinished.disconnect(self.vab_c)
        self.vera_c.setText(str(posn))
        self.vera_c.editingFinished.connect(self.vab_c)
        self.verang_c=posn
        #self.slmupdate()

    def vab_c(self):
        minv=float(self.verang_c_i.text())
        maxv=float(self.verang_c_a.text())
        posn=float(self.vera_c.text())
        self.verang_c_s.valueChanged.disconnect(self.va_c)
        self.verang_c_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.verang_c_s.valueChanged.connect(self.va_c)
        self.verang_c=posn
        #self.slmupdate()

    def fo_c(self):
        minv=float(self.focus_c_i.text())
        maxv=float(self.focus_c_a.text())
        posn=float(self.focus_c_s.value())/100*(maxv-minv)+minv
        self.fv_c.editingFinished.disconnect(self.fob_c)
        self.fv_c.setText(str(posn))
        self.fv_c.editingFinished.connect(self.fob_c)
        self.focus_c=posn
        #self.slmupdate()

    def fob_c(self):
        minv=float(self.focus_c_i.text())
        maxv=float(self.focus_c_a.text())
        posn=float(self.fv_c.text())
        self.focus_c_s.valueChanged.disconnect(self.fo_c)
        self.focus_c_s.setValue(int((posn-minv)/(maxv-minv)*100))
        self.focus_c_s.valueChanged.connect(self.fo_c)
        self.focus_c=posn
        #self.slmupdate()

    def sizer_c(self):
        minv=float(self.sizei_c.text())
        maxv=float(self.sizea_c.text())
        posn=(float(self.sizes_c.value())/100*(maxv-minv)+minv)
        self.sv_c.editingFinished.disconnect(self.sizerb_c)
        self.sv_c.setText(str(posn))
        self.sv_c.editingFinished.connect(self.sizerb_c)
        self.size_c=posn
        #self.slmupdate()

    def sizerb_c(self):
        minv=float(self.sizei_c.text())
        maxv=float(self.sizea_c.text())
        posn=float(self.sv_c.text())
        self.sizes_c.valueChanged.disconnect(self.sizer_c)
        self.sizes_c.setValue(int((posn-minv)/(maxv-minv)*100))
        self.sizes_c.valueChanged.connect(self.sizer_c)
        self.size_c=posn
        #self.slmupdate()

    def timevchange(self):
        self.tm = float(self.timev.text())


    '''
    def slmupdate(self):
        distribution1=np.real((beamgen1.LG(self.l,self.g,self.size)))
        distribution2=np.real((beamgen2.LG(self.l_b,self.g_b,self.size_b)))
        distribution3=np.real((beamgen3.LG(self.l_c,self.g_c,self.size_c)))
        distribution1=beamgen1.focus(distribution1,self.focus)
        distribution2=beamgen2.focus(distribution2,self.focus_b)
        distribution3=beamgen3.focus(distribution3,self.focus_c)
        distribution1=beamgen1.displace(beamgen1.angled(beamgen1.rotator(distribution1,self.rotang), self.verang, self.horang), self.verpos, self.horpos)*intaim
        distribution2=beamgen2.displace(beamgen2.angled(beamgen2.rotator(distribution2,self.rotang_b), self.verang_b, self.horang_b), self.verpos_b, self.horpos_b)*intaim
        distribution3=beamgen3.displace(beamgen3.angled(beamgen3.rotator(distribution3,self.rotang_c), self.verang_c, self.horang_c), self.verpos_c, self.horpos_c)*intaim
        #while True :
            slm.useSLM(distribution1)
            time.sleep(tm)
            slm.useSLM(distribution2)
            time.sleep(tm)
            slm.useSLM(distribution3)
            time.sleep(tm)
        #tosend=slm.useSLM(distribution)
        #print(self.size)
        #scipy.misc.imsave('amp.png', np.abs(distribution))
        #scipy.misc.imsave('phase.png', np.angle(distribution))
        #tosend=tosend.astype('uint16')
        #scipy.misc.imsave('hologram.png', tosend[:,:,0]+np.left_shift(tosend[:,:,1],8))
        #time.sleep(0.1)
        #scipy.misc.imsave('image3.png',currentimg)
        #scipy.misc.imsave('image3.png',currentimg)
        #val, currentimg=compare(currentimg, self.mask, self.cutmask, self.hexcor)
    '''
    def on_send(self):
        distribution1=np.real((beamgen1.LG(self.l,self.g,self.size)))
        distribution2=np.real((beamgen2.LG(self.l_b,self.g_b,self.size_b)))
        distribution3=np.real((beamgen3.LG(self.l_c,self.g_c,self.size_c)))
        distribution1=beamgen1.focus(distribution1,self.focus)
        distribution2=beamgen2.focus(distribution2,self.focus_b)
        distribution3=beamgen3.focus(distribution3,self.focus_c)
        distribution1=beamgen1.displace(beamgen1.angled(beamgen1.rotator(distribution1,self.rotang), self.verang, self.horang), self.verpos, self.horpos)*intaim
        distribution2=beamgen2.displace(beamgen2.angled(beamgen2.rotator(distribution2,self.rotang_b), self.verang_b, self.horang_b), self.verpos_b, self.horpos_b)*intaim
        distribution3=beamgen3.displace(beamgen3.angled(beamgen3.rotator(distribution3,self.rotang_c), self.verang_c, self.horang_c), self.verpos_c, self.horpos_c)*intaim
        dict = {'dist1':distribution1, 'dist2':distribution2, 'dist3':distribution3}
        self.thread1 = QThread1()
        self.sig.connect(self.thread1.on_source)
        self.sig.emit(dict)
        self.sigtm.connect(self.thread1.on_time)
        self.sigtm.emit(self.tm)
        self.thread1.start()
        self.sendb.setEnabled(False)

    def on_pause(self):
        try:
            self.thread1.running = False
            self.sendb.setEnabled(True)
        except:
            pass

    def paramupdate_a(self):
        try:
            F = open('LG'+str(self.l)+str(self.g)+'.txt','r')
            print(F.readline())
            self.horpos=float(F.readline())
            self.horv.setText(str(self.horpos))
            self.horang=float(F.readline())
            self.hora.setText(str(self.horang))
            self.verpos=float(F.readline())
            self.verv.setText(str(self.verpos))
            self.verang=float(F.readline())
            self.vera.setText(str(self.verang))
            self.size=float(F.readline())
            self.sv.setText(str(self.size))
            self.focus = float(F.readline())
            self.fv.setText(str(self.focus))
            self.rotang= float(F.readline())
            self.rv.setText(str(self.rotang))
            F.close()
            self.hpb()
            self.rotb()
            self.vpb()
            self.vab()
            self.fob()
            self.hab()
            self.sizerb()
        except:
            pass

    def paramupdate_b(self):
        try:
            F = open('LG'+str(self.l_b)+str(self.g_b)+'.txt','r')
            print(F.readline())
            self.horpos_b=float(F.readline())
            self.horv_b.setText(str(self.horpos_b))
            self.horang_b=float(F.readline())
            self.hora_b.setText(str(self.horang_b))
            self.verpos_b=float(F.readline())
            self.verv_b.setText(str(self.verpos_b))
            self.verang_b=float(F.readline())
            self.vera_b.setText(str(self.verang_b))
            self.size_b=float(F.readline())
            self.sv_b.setText(str(self.size_b))
            self.focus_b = float(F.readline())
            self.fv_b.setText(str(self.focus_b))
            self.rotang_b= float(F.readline())
            self.rv_b.setText(str(self.rotang_b))
            F.close()
            self.hpb_b()
            self.rotb_b()
            self.vpb_b()
            self.vab_b()
            self.fob_b()
            self.hab_b()
            self.sizerb_b()
        except:
            pass

    def paramupdate_c(self):
        try:
            F = open('LG'+str(self.l_c)+str(self.g_c)+'.txt','r')
            print(F.readline())
            self.horpos_c=float(F.readline())
            self.horv_c.setText(str(self.horpos_c))
            self.horang_c=float(F.readline())
            self.hora_c.setText(str(self.horang_c))
            self.verpos_c=float(F.readline())
            self.verv_c.setText(str(self.verpos_c))
            self.verang_c=float(F.readline())
            self.vera_c.setText(str(self.verang_c))
            self.size_c=float(F.readline())
            self.sv_c.setText(str(self.size_c))
            self.focus_c = float(F.readline())
            self.fv_c.setText(str(self.focus_c))
            self.rotang_c= float(F.readline())
            self.rv_c.setText(str(self.rotang_c))
            F.close()
            self.hpb_c()
            self.rotb_c()
            self.vpb_c()
            self.vab_c()
            self.fob_c()
            self.hab_c()
            self.sizerb_c()
        except:
            pass

    def closeEvent(self, event):
        # here you can terminate your threads and do other stuff
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



if __name__ == '__main__':              # if we're running file directly and not importing it
    main()                              # run the main function
