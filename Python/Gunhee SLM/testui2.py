# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'testui.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!
#UI for testgui2, which adds the functionality of focus and store

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(467, 660)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QtCore.QSize(100, 100))
        MainWindow.setMaximumSize(QtCore.QSize(2000, 2000))
        MainWindow.setWindowOpacity(1.0)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(30, 10, 421, 581))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setObjectName(_fromUtf8("horizontalLayout_12"))
        self.horpl_12 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_12.sizePolicy().hasHeightForWidth())
        self.horpl_12.setSizePolicy(sizePolicy)
        self.horpl_12.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_12.setObjectName(_fromUtf8("horpl_12"))
        self.horizontalLayout_12.addWidget(self.horpl_12)
        self.horv = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horv.sizePolicy().hasHeightForWidth())
        self.horv.setSizePolicy(sizePolicy)
        self.horv.setMaximumSize(QtCore.QSize(70, 16777215))
        self.horv.setObjectName(_fromUtf8("horv"))
        self.horizontalLayout_12.addWidget(self.horv, QtCore.Qt.AlignHCenter)
        self.horpl_18 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_18.sizePolicy().hasHeightForWidth())
        self.horpl_18.setSizePolicy(sizePolicy)
        self.horpl_18.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_18.setText(_fromUtf8(""))
        self.horpl_18.setObjectName(_fromUtf8("horpl_18"))
        self.horizontalLayout_12.addWidget(self.horpl_18)
        self.verticalLayout.addLayout(self.horizontalLayout_12)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.horposi = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horposi.sizePolicy().hasHeightForWidth())
        self.horposi.setSizePolicy(sizePolicy)
        self.horposi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.horposi.setObjectName(_fromUtf8("horposi"))
        self.horizontalLayout.addWidget(self.horposi)
        self.horposs = QtGui.QSlider(self.verticalLayoutWidget)
        self.horposs.setSliderPosition(50)
        self.horposs.setTracking(False)
        self.horposs.setOrientation(QtCore.Qt.Horizontal)
        self.horposs.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horposs.setTickInterval(1)
        self.horposs.setObjectName(_fromUtf8("horposs"))
        self.horizontalLayout.addWidget(self.horposs)
        self.horposa = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horposa.sizePolicy().hasHeightForWidth())
        self.horposa.setSizePolicy(sizePolicy)
        self.horposa.setMaximumSize(QtCore.QSize(50, 16777215))
        self.horposa.setObjectName(_fromUtf8("horposa"))
        self.horizontalLayout.addWidget(self.horposa)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.horpl_8 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_8.sizePolicy().hasHeightForWidth())
        self.horpl_8.setSizePolicy(sizePolicy)
        self.horpl_8.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_8.setObjectName(_fromUtf8("horpl_8"))
        self.horizontalLayout_7.addWidget(self.horpl_8)
        self.hora = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.hora.sizePolicy().hasHeightForWidth())
        self.hora.setSizePolicy(sizePolicy)
        self.hora.setMaximumSize(QtCore.QSize(70, 16777215))
        self.hora.setObjectName(_fromUtf8("hora"))
        self.horizontalLayout_7.addWidget(self.hora, QtCore.Qt.AlignHCenter)
        self.horpl_13 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_13.sizePolicy().hasHeightForWidth())
        self.horpl_13.setSizePolicy(sizePolicy)
        self.horpl_13.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_13.setText(_fromUtf8(""))
        self.horpl_13.setObjectName(_fromUtf8("horpl_13"))
        self.horizontalLayout_7.addWidget(self.horpl_13)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horangi = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horangi.sizePolicy().hasHeightForWidth())
        self.horangi.setSizePolicy(sizePolicy)
        self.horangi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.horangi.setObjectName(_fromUtf8("horangi"))
        self.horizontalLayout_2.addWidget(self.horangi)
        self.horangs = QtGui.QSlider(self.verticalLayoutWidget)
        self.horangs.setSliderPosition(50)
        self.horangs.setTracking(False)
        self.horangs.setOrientation(QtCore.Qt.Horizontal)
        self.horangs.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horangs.setTickInterval(1)
        self.horangs.setObjectName(_fromUtf8("horangs"))
        self.horizontalLayout_2.addWidget(self.horangs)
        self.horanga = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horanga.sizePolicy().hasHeightForWidth())
        self.horanga.setSizePolicy(sizePolicy)
        self.horanga.setMinimumSize(QtCore.QSize(0, 0))
        self.horanga.setMaximumSize(QtCore.QSize(50, 16777215))
        self.horanga.setObjectName(_fromUtf8("horanga"))
        self.horizontalLayout_2.addWidget(self.horanga)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        spacerItem1 = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.horpl_9 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_9.sizePolicy().hasHeightForWidth())
        self.horpl_9.setSizePolicy(sizePolicy)
        self.horpl_9.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_9.setObjectName(_fromUtf8("horpl_9"))
        self.horizontalLayout_8.addWidget(self.horpl_9)
        self.verv = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.verv.sizePolicy().hasHeightForWidth())
        self.verv.setSizePolicy(sizePolicy)
        self.verv.setMaximumSize(QtCore.QSize(70, 16777215))
        self.verv.setObjectName(_fromUtf8("verv"))
        self.horizontalLayout_8.addWidget(self.verv, QtCore.Qt.AlignHCenter)
        self.horpl_14 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_14.sizePolicy().hasHeightForWidth())
        self.horpl_14.setSizePolicy(sizePolicy)
        self.horpl_14.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_14.setText(_fromUtf8(""))
        self.horpl_14.setObjectName(_fromUtf8("horpl_14"))
        self.horizontalLayout_8.addWidget(self.horpl_14)
        self.verticalLayout.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.verposi = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.verposi.sizePolicy().hasHeightForWidth())
        self.verposi.setSizePolicy(sizePolicy)
        self.verposi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.verposi.setObjectName(_fromUtf8("verposi"))
        self.horizontalLayout_3.addWidget(self.verposi)
        self.verposs = QtGui.QSlider(self.verticalLayoutWidget)
        self.verposs.setSliderPosition(50)
        self.verposs.setTracking(False)
        self.verposs.setOrientation(QtCore.Qt.Horizontal)
        self.verposs.setTickPosition(QtGui.QSlider.TicksBelow)
        self.verposs.setTickInterval(1)
        self.verposs.setObjectName(_fromUtf8("verposs"))
        self.horizontalLayout_3.addWidget(self.verposs)
        self.verposa = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.verposa.sizePolicy().hasHeightForWidth())
        self.verposa.setSizePolicy(sizePolicy)
        self.verposa.setMinimumSize(QtCore.QSize(0, 0))
        self.verposa.setMaximumSize(QtCore.QSize(50, 16777215))
        self.verposa.setObjectName(_fromUtf8("verposa"))
        self.horizontalLayout_3.addWidget(self.verposa)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        spacerItem2 = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.horpl_10 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_10.sizePolicy().hasHeightForWidth())
        self.horpl_10.setSizePolicy(sizePolicy)
        self.horpl_10.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_10.setObjectName(_fromUtf8("horpl_10"))
        self.horizontalLayout_9.addWidget(self.horpl_10)
        self.vera = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.vera.sizePolicy().hasHeightForWidth())
        self.vera.setSizePolicy(sizePolicy)
        self.vera.setMaximumSize(QtCore.QSize(70, 16777215))
        self.vera.setObjectName(_fromUtf8("vera"))
        self.horizontalLayout_9.addWidget(self.vera, QtCore.Qt.AlignHCenter)
        self.horpl_15 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_15.sizePolicy().hasHeightForWidth())
        self.horpl_15.setSizePolicy(sizePolicy)
        self.horpl_15.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_15.setText(_fromUtf8(""))
        self.horpl_15.setObjectName(_fromUtf8("horpl_15"))
        self.horizontalLayout_9.addWidget(self.horpl_15)
        self.verticalLayout.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.verangi = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.verangi.sizePolicy().hasHeightForWidth())
        self.verangi.setSizePolicy(sizePolicy)
        self.verangi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.verangi.setObjectName(_fromUtf8("verangi"))
        self.horizontalLayout_4.addWidget(self.verangi)
        self.verangs = QtGui.QSlider(self.verticalLayoutWidget)
        self.verangs.setSliderPosition(50)
        self.verangs.setTracking(False)
        self.verangs.setOrientation(QtCore.Qt.Horizontal)
        self.verangs.setTickPosition(QtGui.QSlider.TicksBelow)
        self.verangs.setTickInterval(1)
        self.verangs.setObjectName(_fromUtf8("verangs"))
        self.horizontalLayout_4.addWidget(self.verangs)
        self.veranga = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.veranga.sizePolicy().hasHeightForWidth())
        self.veranga.setSizePolicy(sizePolicy)
        self.veranga.setMinimumSize(QtCore.QSize(0, 0))
        self.veranga.setMaximumSize(QtCore.QSize(50, 16777215))
        self.veranga.setObjectName(_fromUtf8("veranga"))
        self.horizontalLayout_4.addWidget(self.veranga)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        spacerItem3 = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem3)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.horpl_11 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_11.sizePolicy().hasHeightForWidth())
        self.horpl_11.setSizePolicy(sizePolicy)
        self.horpl_11.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_11.setObjectName(_fromUtf8("horpl_11"))
        self.horizontalLayout_10.addWidget(self.horpl_11)
        self.sv = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sv.sizePolicy().hasHeightForWidth())
        self.sv.setSizePolicy(sizePolicy)
        self.sv.setMaximumSize(QtCore.QSize(70, 16777215))
        self.sv.setObjectName(_fromUtf8("sv"))
        self.horizontalLayout_10.addWidget(self.sv, QtCore.Qt.AlignHCenter)
        self.horpl_16 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_16.sizePolicy().hasHeightForWidth())
        self.horpl_16.setSizePolicy(sizePolicy)
        self.horpl_16.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_16.setText(_fromUtf8(""))
        self.horpl_16.setObjectName(_fromUtf8("horpl_16"))
        self.horizontalLayout_10.addWidget(self.horpl_16)
        self.verticalLayout.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.sizei = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizei.sizePolicy().hasHeightForWidth())
        self.sizei.setSizePolicy(sizePolicy)
        self.sizei.setMaximumSize(QtCore.QSize(50, 16777215))
        self.sizei.setObjectName(_fromUtf8("sizei"))
        self.horizontalLayout_5.addWidget(self.sizei)
        self.sizes = QtGui.QSlider(self.verticalLayoutWidget)
        self.sizes.setSliderPosition(50)
        self.sizes.setTracking(False)
        self.sizes.setOrientation(QtCore.Qt.Horizontal)
        self.sizes.setTickPosition(QtGui.QSlider.TicksBelow)
        self.sizes.setTickInterval(1)
        self.sizes.setObjectName(_fromUtf8("sizes"))
        self.horizontalLayout_5.addWidget(self.sizes)
        self.sizea = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizea.sizePolicy().hasHeightForWidth())
        self.sizea.setSizePolicy(sizePolicy)
        self.sizea.setMinimumSize(QtCore.QSize(0, 0))
        self.sizea.setMaximumSize(QtCore.QSize(50, 16777215))
        self.sizea.setObjectName(_fromUtf8("sizea"))
        self.horizontalLayout_5.addWidget(self.sizea)
        self.verticalLayout.addLayout(self.horizontalLayout_5)



        spacerItem4 = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem4)
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName(_fromUtf8("horizontalLayout_20"))
        self.horpl_20 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_20.sizePolicy().hasHeightForWidth())
        self.horpl_20.setSizePolicy(sizePolicy)
        self.horpl_20.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_20.setObjectName(_fromUtf8("horpl_20"))
        self.horizontalLayout_20.addWidget(self.horpl_20)
        self.fv = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fv.sizePolicy().hasHeightForWidth())
        self.fv.setSizePolicy(sizePolicy)
        self.fv.setMaximumSize(QtCore.QSize(70, 16777215))
        self.fv.setObjectName(_fromUtf8("fv"))
        self.horizontalLayout_20.addWidget(self.fv, QtCore.Qt.AlignHCenter)
        self.horpl_21 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_21.sizePolicy().hasHeightForWidth())
        self.horpl_21.setSizePolicy(sizePolicy)
        self.horpl_21.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_21.setText(_fromUtf8(""))
        self.horpl_21.setObjectName(_fromUtf8("horpl_21"))
        self.horizontalLayout_20.addWidget(self.horpl_21)
        self.verticalLayout.addLayout(self.horizontalLayout_20)
        self.horizontalLayout_21 = QtGui.QHBoxLayout()
        self.horizontalLayout_21.setObjectName(_fromUtf8("horizontalLayout_21"))
        self.focusi = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.focusi.sizePolicy().hasHeightForWidth())
        self.focusi.setSizePolicy(sizePolicy)
        self.focusi.setMaximumSize(QtCore.QSize(50, 16777215))
        self.focusi.setObjectName(_fromUtf8("focusi"))
        self.horizontalLayout_21.addWidget(self.focusi)
        self.focuss = QtGui.QSlider(self.verticalLayoutWidget)
        self.focuss.setSliderPosition(50)
        self.focuss.setTracking(False)
        self.focuss.setOrientation(QtCore.Qt.Horizontal)
        self.focuss.setTickPosition(QtGui.QSlider.TicksBelow)
        self.focuss.setTickInterval(1)
        self.focuss.setObjectName(_fromUtf8("focuss"))
        self.horizontalLayout_21.addWidget(self.focuss)
        self.focusa = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.focusa.sizePolicy().hasHeightForWidth())
        self.focusa.setSizePolicy(sizePolicy)
        self.focusa.setMinimumSize(QtCore.QSize(0, 0))
        self.focusa.setMaximumSize(QtCore.QSize(50, 16777215))
        self.focusa.setObjectName(_fromUtf8("focusa"))
        self.horizontalLayout_21.addWidget(self.focusa)
        self.verticalLayout.addLayout(self.horizontalLayout_21)




        spacerItem5 = QtGui.QSpacerItem(416, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem5)
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setObjectName(_fromUtf8("horizontalLayout_19"))
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.horizontalLayout_18 = QtGui.QHBoxLayout()
        self.horizontalLayout_18.setObjectName(_fromUtf8("horizontalLayout_18"))
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.label_2 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_7.addWidget(self.label_2)
        self.lbox = QtGui.QSpinBox(self.verticalLayoutWidget)
        self.lbox.setObjectName(_fromUtf8("lbox"))
        self.verticalLayout_7.addWidget(self.lbox)
        spacerItem6 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_7.addItem(spacerItem6)
        self.horizontalLayout_18.addLayout(self.verticalLayout_7)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.label_3 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_8.addWidget(self.label_3)
        self.gbox = QtGui.QSpinBox(self.verticalLayoutWidget)
        self.gbox.setObjectName(_fromUtf8("gbox"))
        self.verticalLayout_8.addWidget(self.gbox)
        spacerItem7 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_8.addItem(spacerItem7)
        self.horizontalLayout_18.addLayout(self.verticalLayout_8)
        self.verticalLayout_6.addLayout(self.horizontalLayout_18)
        self.horizontalLayout_19.addLayout(self.verticalLayout_6)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setObjectName(_fromUtf8("horizontalLayout_11"))
        self.horpl_6 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_6.sizePolicy().hasHeightForWidth())
        self.horpl_6.setSizePolicy(sizePolicy)
        self.horpl_6.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_6.setObjectName(_fromUtf8("horpl_6"))
        self.horizontalLayout_11.addWidget(self.horpl_6)
        self.rv = QtGui.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.rv.sizePolicy().hasHeightForWidth())
        self.rv.setSizePolicy(sizePolicy)
        self.rv.setMaximumSize(QtCore.QSize(70, 16777215))
        self.rv.setObjectName(_fromUtf8("rv"))
        self.horizontalLayout_11.addWidget(self.rv, QtCore.Qt.AlignHCenter)
        self.horpl_17 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.horpl_17.sizePolicy().hasHeightForWidth())
        self.horpl_17.setSizePolicy(sizePolicy)
        self.horpl_17.setMinimumSize(QtCore.QSize(100, 0))
        self.horpl_17.setText(_fromUtf8(""))
        self.horpl_17.setObjectName(_fromUtf8("horpl_17"))
        self.horizontalLayout_11.addWidget(self.horpl_17)
        self.verticalLayout_3.addLayout(self.horizontalLayout_11)
        self.dial = QtGui.QDial(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dial.sizePolicy().hasHeightForWidth())
        self.dial.setSizePolicy(sizePolicy)
        self.dial.setMaximumSize(QtCore.QSize(16777215, 1000))
        self.dial.setMaximum(360)
        self.dial.setProperty("value", 0)
        self.dial.setTracking(False)
        self.dial.setOrientation(QtCore.Qt.Vertical)
        self.dial.setWrapping(True)
        self.dial.setNotchTarget(5.0)
        self.dial.setNotchesVisible(True)
        self.dial.setObjectName(_fromUtf8("dial"))
        self.verticalLayout_3.addWidget(self.dial)
        self.horizontalLayout_19.addLayout(self.verticalLayout_3)
        self.verticalLayout.addLayout(self.horizontalLayout_19)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.verticalLayout.addLayout(self.horizontalLayout_6)

        self.storb = QtGui.QPushButton(self.centralwidget)
        self.storb.setGeometry(QtCore.QRect(30, 620, 100, 30))
        self.storb.setObjectName("storb")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "SLM Control", None))
        self.horpl_12.setText(_translate("MainWindow", "Horizontal Position:", None))
        self.horv.setText(_translate("MainWindow", "0", None))
        self.horposi.setText(_translate("MainWindow", "-50", None))
        self.horposa.setText(_translate("MainWindow", "50", None))
        self.horpl_8.setText(_translate("MainWindow", "Horizontal Angle:", None))
        self.hora.setText(_translate("MainWindow", "0", None))
        self.horangi.setText(_translate("MainWindow", "-15", None))
        self.horanga.setText(_translate("MainWindow", "15", None))
        self.horpl_9.setText(_translate("MainWindow", "Vertical Position:", None))
        self.verv.setText(_translate("MainWindow", "0", None))
        self.verposi.setText(_translate("MainWindow", "-50", None))
        self.verposa.setText(_translate("MainWindow", "50", None))
        self.horpl_10.setText(_translate("MainWindow", "Vertical Angle:", None))
        self.vera.setText(_translate("MainWindow", "0", None))
        self.verangi.setText(_translate("MainWindow", "-15", None))
        self.veranga.setText(_translate("MainWindow", "15", None))
        self.horpl_20.setText(_translate("MainWindow", "Focus:", None))
        self.fv.setText(_translate("MainWindow", "0", None))
        self.focusi.setText(_translate("MainWindow", "-0.1", None))
        self.focusa.setText(_translate("MainWindow", "0.1", None))
        self.horpl_11.setText(_translate("MainWindow", "Size:", None))
        self.sv.setText(_translate("MainWindow", "0", None))
        self.sizei.setText(_translate("MainWindow", "0", None))
        self.sizea.setText(_translate("MainWindow", "300", None))
        self.label_2.setText(_translate("MainWindow", "L:", None))
        self.label_3.setText(_translate("MainWindow", "G:", None))
        self.horpl_6.setText(_translate("MainWindow", "Rotation:", None))
        self.rv.setText(_translate("MainWindow", "0", None))
        self.storb.setText(_translate("MainWindow", "Store",None))
