# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'interface.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1017, 659)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.tabWidget_2 = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.gridLayout_19 = QtWidgets.QGridLayout(self.tab_3)
        self.gridLayout_19.setObjectName("gridLayout_19")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.swave = QtWidgets.QLineEdit(self.tab_3)
        self.swave.setObjectName("swave")
        self.gridLayout.addWidget(self.swave, 13, 1, 1, 1)
        self.Bz = QtWidgets.QLineEdit(self.tab_3)
        self.Bz.setObjectName("Bz")
        self.gridLayout.addWidget(self.Bz, 5, 1, 1, 1)
        self.Bx = QtWidgets.QLineEdit(self.tab_3)
        self.Bx.setObjectName("Bx")
        self.gridLayout.addWidget(self.Bx, 3, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.tab_3)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 5, 0, 1, 1)
        self.mAB = QtWidgets.QLineEdit(self.tab_3)
        self.mAB.setObjectName("mAB")
        self.gridLayout.addWidget(self.mAB, 11, 1, 1, 1)
        self.By = QtWidgets.QLineEdit(self.tab_3)
        self.By.setObjectName("By")
        self.gridLayout.addWidget(self.By, 4, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.tab_3)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.peierls = QtWidgets.QLineEdit(self.tab_3)
        self.peierls.setEnabled(True)
        self.peierls.setObjectName("peierls")
        self.gridLayout.addWidget(self.peierls, 2, 1, 1, 1)
        self.kanemele = QtWidgets.QLineEdit(self.tab_3)
        self.kanemele.setObjectName("kanemele")
        self.gridLayout.addWidget(self.kanemele, 7, 1, 1, 1)
        self.label_24 = QtWidgets.QLabel(self.tab_3)
        self.label_24.setObjectName("label_24")
        self.gridLayout.addWidget(self.label_24, 9, 0, 1, 1)
        self.antihaldane = QtWidgets.QLineEdit(self.tab_3)
        self.antihaldane.setObjectName("antihaldane")
        self.gridLayout.addWidget(self.antihaldane, 10, 1, 1, 1)
        self.label_26 = QtWidgets.QLabel(self.tab_3)
        self.label_26.setObjectName("label_26")
        self.gridLayout.addWidget(self.label_26, 13, 0, 1, 1)
        self.fermi = QtWidgets.QLineEdit(self.tab_3)
        self.fermi.setEnabled(True)
        self.fermi.setObjectName("fermi")
        self.gridLayout.addWidget(self.fermi, 0, 1, 1, 1)
        self.rashba = QtWidgets.QLineEdit(self.tab_3)
        self.rashba.setObjectName("rashba")
        self.gridLayout.addWidget(self.rashba, 6, 1, 1, 1)
        self.haldane = QtWidgets.QLineEdit(self.tab_3)
        self.haldane.setObjectName("haldane")
        self.gridLayout.addWidget(self.haldane, 9, 1, 1, 1)
        self.mAF = QtWidgets.QLineEdit(self.tab_3)
        self.mAF.setObjectName("mAF")
        self.gridLayout.addWidget(self.mAF, 12, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.tab_3)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 6, 0, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.tab_3)
        self.label_11.setObjectName("label_11")
        self.gridLayout.addWidget(self.label_11, 7, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.tab_3)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.tab_3)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.tab_3)
        self.label_12.setObjectName("label_12")
        self.gridLayout.addWidget(self.label_12, 11, 0, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.tab_3)
        self.label_14.setObjectName("label_14")
        self.gridLayout.addWidget(self.label_14, 2, 0, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.tab_3)
        self.label_13.setObjectName("label_13")
        self.gridLayout.addWidget(self.label_13, 12, 0, 1, 1)
        self.label_25 = QtWidgets.QLabel(self.tab_3)
        self.label_25.setObjectName("label_25")
        self.gridLayout.addWidget(self.label_25, 10, 0, 1, 1)
        self.antikanemele = QtWidgets.QLineEdit(self.tab_3)
        self.antikanemele.setObjectName("antikanemele")
        self.gridLayout.addWidget(self.antikanemele, 8, 1, 1, 1)
        self.label_37 = QtWidgets.QLabel(self.tab_3)
        self.label_37.setObjectName("label_37")
        self.gridLayout.addWidget(self.label_37, 8, 0, 1, 1)
        self.crystalfield = QtWidgets.QLineEdit(self.tab_3)
        self.crystalfield.setObjectName("crystalfield")
        self.gridLayout.addWidget(self.crystalfield, 1, 1, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.tab_3)
        self.label_17.setObjectName("label_17")
        self.gridLayout.addWidget(self.label_17, 1, 0, 1, 1)
        self.gridLayout_19.addLayout(self.gridLayout, 0, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_3, "")
        self.gridLayout_7.addWidget(self.tabWidget_2, 0, 0, 2, 1)
        self.tabWidget_3 = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget_3.setObjectName("tabWidget_3")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.tab_4)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.show_structure = QtWidgets.QPushButton(self.tab_4)
        self.show_structure.setObjectName("show_structure")
        self.gridLayout_9.addWidget(self.show_structure, 1, 0, 1, 1)
        self.show_structure_3d = QtWidgets.QPushButton(self.tab_4)
        self.show_structure_3d.setObjectName("show_structure_3d")
        self.gridLayout_9.addWidget(self.show_structure_3d, 1, 1, 1, 1)
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.nsuper_struct = QtWidgets.QLineEdit(self.tab_4)
        self.nsuper_struct.setObjectName("nsuper_struct")
        self.gridLayout_3.addWidget(self.nsuper_struct, 0, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.tab_4)
        self.label_7.setObjectName("label_7")
        self.gridLayout_3.addWidget(self.label_7, 0, 0, 1, 1)
        self.gridLayout_9.addLayout(self.gridLayout_3, 0, 0, 1, 2)
        self.tabWidget_3.addTab(self.tab_4, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.tab_5)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.show_bands = QtWidgets.QPushButton(self.tab_5)
        self.show_bands.setObjectName("show_bands")
        self.gridLayout_8.addWidget(self.show_bands, 0, 0, 1, 1)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_9 = QtWidgets.QLabel(self.tab_5)
        self.label_9.setObjectName("label_9")
        self.gridLayout_2.addWidget(self.label_9, 1, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.tab_5)
        self.label_15.setObjectName("label_15")
        self.gridLayout_2.addWidget(self.label_15, 0, 0, 1, 1)
        self.bands_color = QtWidgets.QComboBox(self.tab_5)
        self.bands_color.setObjectName("bands_color")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.bands_color.addItem("")
        self.gridLayout_2.addWidget(self.bands_color, 0, 1, 1, 1)
        self.nk_bands = QtWidgets.QLineEdit(self.tab_5)
        self.nk_bands.setObjectName("nk_bands")
        self.gridLayout_2.addWidget(self.nk_bands, 1, 1, 1, 1)
        self.gridLayout_8.addLayout(self.gridLayout_2, 1, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_5, "")
        self.tab_9 = QtWidgets.QWidget()
        self.tab_9.setObjectName("tab_9")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.tab_9)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gridLayout_11 = QtWidgets.QGridLayout()
        self.gridLayout_11.setObjectName("gridLayout_11")
        self.delta_kbands = QtWidgets.QLineEdit(self.tab_9)
        self.delta_kbands.setObjectName("delta_kbands")
        self.gridLayout_11.addWidget(self.delta_kbands, 0, 1, 1, 1)
        self.label_27 = QtWidgets.QLabel(self.tab_9)
        self.label_27.setObjectName("label_27")
        self.gridLayout_11.addWidget(self.label_27, 0, 0, 1, 1)
        self.ne_kbands = QtWidgets.QLineEdit(self.tab_9)
        self.ne_kbands.setObjectName("ne_kbands")
        self.gridLayout_11.addWidget(self.ne_kbands, 1, 1, 1, 1)
        self.label_28 = QtWidgets.QLabel(self.tab_9)
        self.label_28.setObjectName("label_28")
        self.gridLayout_11.addWidget(self.label_28, 1, 0, 1, 1)
        self.label_29 = QtWidgets.QLabel(self.tab_9)
        self.label_29.setObjectName("label_29")
        self.gridLayout_11.addWidget(self.label_29, 2, 0, 1, 1)
        self.window_kbands = QtWidgets.QLineEdit(self.tab_9)
        self.window_kbands.setObjectName("window_kbands")
        self.gridLayout_11.addWidget(self.window_kbands, 2, 1, 1, 1)
        self.label_30 = QtWidgets.QLabel(self.tab_9)
        self.label_30.setObjectName("label_30")
        self.gridLayout_11.addWidget(self.label_30, 3, 0, 1, 1)
        self.scale_kbands = QtWidgets.QLineEdit(self.tab_9)
        self.scale_kbands.setObjectName("scale_kbands")
        self.gridLayout_11.addWidget(self.scale_kbands, 3, 1, 1, 1)
        self.label_31 = QtWidgets.QLabel(self.tab_9)
        self.label_31.setObjectName("label_31")
        self.gridLayout_11.addWidget(self.label_31, 4, 0, 1, 1)
        self.nv_kbands = QtWidgets.QLineEdit(self.tab_9)
        self.nv_kbands.setObjectName("nv_kbands")
        self.gridLayout_11.addWidget(self.nv_kbands, 4, 1, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_11)
        self.show_dosbands = QtWidgets.QPushButton(self.tab_9)
        self.show_dosbands.setObjectName("show_dosbands")
        self.verticalLayout_2.addWidget(self.show_dosbands)
        self.tabWidget_3.addTab(self.tab_9, "")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.gridLayout_14 = QtWidgets.QGridLayout(self.tab_6)
        self.gridLayout_14.setObjectName("gridLayout_14")
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label_16 = QtWidgets.QLabel(self.tab_6)
        self.label_16.setObjectName("label_16")
        self.gridLayout_5.addWidget(self.label_16, 0, 0, 1, 1)
        self.dos_delta = QtWidgets.QLineEdit(self.tab_6)
        self.dos_delta.setObjectName("dos_delta")
        self.gridLayout_5.addWidget(self.dos_delta, 0, 1, 1, 1)
        self.label_40 = QtWidgets.QLabel(self.tab_6)
        self.label_40.setObjectName("label_40")
        self.gridLayout_5.addWidget(self.label_40, 1, 0, 1, 1)
        self.dos_nk = QtWidgets.QLineEdit(self.tab_6)
        self.dos_nk.setObjectName("dos_nk")
        self.gridLayout_5.addWidget(self.dos_nk, 1, 1, 1, 1)
        self.dos_ewindow = QtWidgets.QLineEdit(self.tab_6)
        self.dos_ewindow.setObjectName("dos_ewindow")
        self.gridLayout_5.addWidget(self.dos_ewindow, 2, 1, 1, 1)
        self.label_41 = QtWidgets.QLabel(self.tab_6)
        self.label_41.setObjectName("label_41")
        self.gridLayout_5.addWidget(self.label_41, 2, 0, 1, 1)
        self.gridLayout_14.addLayout(self.gridLayout_5, 1, 0, 1, 1)
        self.show_dos = QtWidgets.QPushButton(self.tab_6)
        self.show_dos.setObjectName("show_dos")
        self.gridLayout_14.addWidget(self.show_dos, 2, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_6, "")
        self.tab_7 = QtWidgets.QWidget()
        self.tab_7.setObjectName("tab_7")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.tab_7)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.gridLayout_22 = QtWidgets.QGridLayout()
        self.gridLayout_22.setObjectName("gridLayout_22")
        self.multildos_ewindow = QtWidgets.QLineEdit(self.tab_7)
        self.multildos_ewindow.setObjectName("multildos_ewindow")
        self.gridLayout_22.addWidget(self.multildos_ewindow, 0, 1, 1, 1)
        self.multildos_delta = QtWidgets.QLineEdit(self.tab_7)
        self.multildos_delta.setObjectName("multildos_delta")
        self.gridLayout_22.addWidget(self.multildos_delta, 3, 1, 1, 1)
        self.label_36 = QtWidgets.QLabel(self.tab_7)
        self.label_36.setObjectName("label_36")
        self.gridLayout_22.addWidget(self.label_36, 3, 0, 1, 1)
        self.label_42 = QtWidgets.QLabel(self.tab_7)
        self.label_42.setObjectName("label_42")
        self.gridLayout_22.addWidget(self.label_42, 0, 0, 1, 1)
        self.multildos_nk = QtWidgets.QLineEdit(self.tab_7)
        self.multildos_nk.setObjectName("multildos_nk")
        self.gridLayout_22.addWidget(self.multildos_nk, 1, 1, 1, 1)
        self.show_multildos = QtWidgets.QPushButton(self.tab_7)
        self.show_multildos.setObjectName("show_multildos")
        self.gridLayout_22.addWidget(self.show_multildos, 6, 0, 1, 2)
        self.label_43 = QtWidgets.QLabel(self.tab_7)
        self.label_43.setObjectName("label_43")
        self.gridLayout_22.addWidget(self.label_43, 1, 0, 1, 1)
        self.label_44 = QtWidgets.QLabel(self.tab_7)
        self.label_44.setObjectName("label_44")
        self.gridLayout_22.addWidget(self.label_44, 2, 0, 1, 1)
        self.multildos_nrep = QtWidgets.QLineEdit(self.tab_7)
        self.multildos_nrep.setObjectName("multildos_nrep")
        self.gridLayout_22.addWidget(self.multildos_nrep, 2, 1, 1, 1)
        self.basis_ldos = QtWidgets.QComboBox(self.tab_7)
        self.basis_ldos.setObjectName("basis_ldos")
        self.basis_ldos.addItem("")
        self.basis_ldos.addItem("")
        self.gridLayout_22.addWidget(self.basis_ldos, 4, 1, 1, 1)
        self.ratomic_ldos = QtWidgets.QLineEdit(self.tab_7)
        self.ratomic_ldos.setObjectName("ratomic_ldos")
        self.gridLayout_22.addWidget(self.ratomic_ldos, 5, 1, 1, 1)
        self.label_19 = QtWidgets.QLabel(self.tab_7)
        self.label_19.setObjectName("label_19")
        self.gridLayout_22.addWidget(self.label_19, 4, 0, 1, 1)
        self.label_20 = QtWidgets.QLabel(self.tab_7)
        self.label_20.setObjectName("label_20")
        self.gridLayout_22.addWidget(self.label_20, 5, 0, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout_22, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_7, "")
        self.tab_13 = QtWidgets.QWidget()
        self.tab_13.setObjectName("tab_13")
        self.gridLayout_13 = QtWidgets.QGridLayout(self.tab_13)
        self.gridLayout_13.setObjectName("gridLayout_13")
        self.gridLayout_23 = QtWidgets.QGridLayout()
        self.gridLayout_23.setObjectName("gridLayout_23")
        self.ldos_delta = QtWidgets.QLineEdit(self.tab_13)
        self.ldos_delta.setObjectName("ldos_delta")
        self.gridLayout_23.addWidget(self.ldos_delta, 2, 1, 1, 1)
        self.label_38 = QtWidgets.QLabel(self.tab_13)
        self.label_38.setObjectName("label_38")
        self.gridLayout_23.addWidget(self.label_38, 2, 0, 1, 1)
        self.label_45 = QtWidgets.QLabel(self.tab_13)
        self.label_45.setObjectName("label_45")
        self.gridLayout_23.addWidget(self.label_45, 0, 0, 1, 1)
        self.ldos_nk = QtWidgets.QLineEdit(self.tab_13)
        self.ldos_nk.setObjectName("ldos_nk")
        self.gridLayout_23.addWidget(self.ldos_nk, 1, 1, 1, 1)
        self.show_ldos = QtWidgets.QPushButton(self.tab_13)
        self.show_ldos.setObjectName("show_ldos")
        self.gridLayout_23.addWidget(self.show_ldos, 4, 0, 1, 2)
        self.ldos_ewindow = QtWidgets.QLineEdit(self.tab_13)
        self.ldos_ewindow.setObjectName("ldos_ewindow")
        self.gridLayout_23.addWidget(self.ldos_ewindow, 0, 1, 1, 1)
        self.label_46 = QtWidgets.QLabel(self.tab_13)
        self.label_46.setObjectName("label_46")
        self.gridLayout_23.addWidget(self.label_46, 1, 0, 1, 1)
        self.ldos_operator = QtWidgets.QComboBox(self.tab_13)
        self.ldos_operator.setObjectName("ldos_operator")
        self.ldos_operator.addItem("")
        self.ldos_operator.addItem("")
        self.gridLayout_23.addWidget(self.ldos_operator, 3, 1, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.tab_13)
        self.label_18.setObjectName("label_18")
        self.gridLayout_23.addWidget(self.label_18, 3, 0, 1, 1)
        self.gridLayout_13.addLayout(self.gridLayout_23, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_13, "")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_16 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_16.setObjectName("gridLayout_16")
        self.tabWidget_4 = QtWidgets.QTabWidget(self.tab)
        self.tabWidget_4.setObjectName("tabWidget_4")
        self.tab_12 = QtWidgets.QWidget()
        self.tab_12.setObjectName("tab_12")
        self.gridLayout_18 = QtWidgets.QGridLayout(self.tab_12)
        self.gridLayout_18.setObjectName("gridLayout_18")
        self.gridLayout_10 = QtWidgets.QGridLayout()
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.scf_initialization = QtWidgets.QComboBox(self.tab_12)
        self.scf_initialization.setObjectName("scf_initialization")
        self.scf_initialization.addItem("")
        self.scf_initialization.addItem("")
        self.scf_initialization.addItem("")
        self.gridLayout_10.addWidget(self.scf_initialization, 0, 1, 1, 1)
        self.label_22 = QtWidgets.QLabel(self.tab_12)
        self.label_22.setObjectName("label_22")
        self.gridLayout_10.addWidget(self.label_22, 0, 0, 1, 1)
        self.label_23 = QtWidgets.QLabel(self.tab_12)
        self.label_23.setObjectName("label_23")
        self.gridLayout_10.addWidget(self.label_23, 1, 0, 1, 1)
        self.hubbard = QtWidgets.QLineEdit(self.tab_12)
        self.hubbard.setObjectName("hubbard")
        self.gridLayout_10.addWidget(self.hubbard, 1, 1, 1, 1)
        self.label_34 = QtWidgets.QLabel(self.tab_12)
        self.label_34.setObjectName("label_34")
        self.gridLayout_10.addWidget(self.label_34, 2, 0, 1, 1)
        self.filling_scf = QtWidgets.QLineEdit(self.tab_12)
        self.filling_scf.setObjectName("filling_scf")
        self.gridLayout_10.addWidget(self.filling_scf, 2, 1, 1, 1)
        self.gridLayout_18.addLayout(self.gridLayout_10, 0, 0, 1, 2)
        self.do_scf = QtWidgets.QCheckBox(self.tab_12)
        self.do_scf.setObjectName("do_scf")
        self.gridLayout_18.addWidget(self.do_scf, 1, 0, 1, 1)
        self.solve_scf = QtWidgets.QPushButton(self.tab_12)
        self.solve_scf.setObjectName("solve_scf")
        self.gridLayout_18.addWidget(self.solve_scf, 1, 1, 1, 1)
        self.tabWidget_4.addTab(self.tab_12, "")
        self.tab_11 = QtWidgets.QWidget()
        self.tab_11.setObjectName("tab_11")
        self.gridLayoutWidget_12 = QtWidgets.QWidget(self.tab_11)
        self.gridLayoutWidget_12.setGeometry(QtCore.QRect(30, 30, 215, 161))
        self.gridLayoutWidget_12.setObjectName("gridLayoutWidget_12")
        self.gridLayout_12 = QtWidgets.QGridLayout(self.gridLayoutWidget_12)
        self.gridLayout_12.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_12.setObjectName("gridLayout_12")
        self.label_32 = QtWidgets.QLabel(self.gridLayoutWidget_12)
        self.label_32.setObjectName("label_32")
        self.gridLayout_12.addWidget(self.label_32, 0, 0, 1, 1)
        self.label_33 = QtWidgets.QLabel(self.gridLayoutWidget_12)
        self.label_33.setObjectName("label_33")
        self.gridLayout_12.addWidget(self.label_33, 1, 0, 1, 1)
        self.nk_scf = QtWidgets.QLineEdit(self.gridLayoutWidget_12)
        self.nk_scf.setObjectName("nk_scf")
        self.gridLayout_12.addWidget(self.nk_scf, 1, 1, 1, 1)
        self.mix_scf = QtWidgets.QLineEdit(self.gridLayoutWidget_12)
        self.mix_scf.setObjectName("mix_scf")
        self.gridLayout_12.addWidget(self.mix_scf, 0, 1, 1, 1)
        self.label_35 = QtWidgets.QLabel(self.gridLayoutWidget_12)
        self.label_35.setObjectName("label_35")
        self.gridLayout_12.addWidget(self.label_35, 2, 0, 1, 1)
        self.smearing_scf = QtWidgets.QLineEdit(self.gridLayoutWidget_12)
        self.smearing_scf.setObjectName("smearing_scf")
        self.gridLayout_12.addWidget(self.smearing_scf, 2, 1, 1, 1)
        self.tabWidget_4.addTab(self.tab_11, "")
        self.gridLayout_16.addWidget(self.tabWidget_4, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab, "")
        self.tab_8 = QtWidgets.QWidget()
        self.tab_8.setObjectName("tab_8")
        self.gridLayout_17 = QtWidgets.QGridLayout(self.tab_8)
        self.gridLayout_17.setObjectName("gridLayout_17")
        self.label_39 = QtWidgets.QLabel(self.tab_8)
        self.label_39.setObjectName("label_39")
        self.gridLayout_17.addWidget(self.label_39, 0, 0, 1, 1)
        self.magnetization_nrep = QtWidgets.QLineEdit(self.tab_8)
        self.magnetization_nrep.setObjectName("magnetization_nrep")
        self.gridLayout_17.addWidget(self.magnetization_nrep, 0, 1, 1, 1)
        self.show_magnetism = QtWidgets.QPushButton(self.tab_8)
        self.show_magnetism.setObjectName("show_magnetism")
        self.gridLayout_17.addWidget(self.show_magnetism, 1, 0, 1, 2)
        self.tabWidget_3.addTab(self.tab_8, "")
        self.gridLayout_7.addWidget(self.tabWidget_3, 0, 1, 1, 1)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setDocumentMode(False)
        self.tabWidget.setTabsClosable(False)
        self.tabWidget.setMovable(False)
        self.tabWidget.setObjectName("tabWidget")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_20 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_20.setObjectName("gridLayout_20")
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_8 = QtWidgets.QLabel(self.tab_2)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 0, 0, 1, 1)
        self.lattice = QtWidgets.QComboBox(self.tab_2)
        self.lattice.setObjectName("lattice")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.lattice.addItem("")
        self.gridLayout_4.addWidget(self.lattice, 0, 1, 1, 1)
        self.nsuper = QtWidgets.QLineEdit(self.tab_2)
        self.nsuper.setObjectName("nsuper")
        self.gridLayout_4.addWidget(self.nsuper, 1, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.tab_2)
        self.label_6.setObjectName("label_6")
        self.gridLayout_4.addWidget(self.label_6, 1, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.tab_2)
        self.label_2.setObjectName("label_2")
        self.gridLayout_4.addWidget(self.label_2, 2, 0, 1, 1)
        self.width = QtWidgets.QLineEdit(self.tab_2)
        self.width.setObjectName("width")
        self.gridLayout_4.addWidget(self.width, 2, 1, 1, 1)
        self.gridLayout_20.addLayout(self.gridLayout_4, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_10 = QtWidgets.QWidget()
        self.tab_10.setObjectName("tab_10")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.tab_10)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout_21 = QtWidgets.QGridLayout()
        self.gridLayout_21.setObjectName("gridLayout_21")
        self.remove_single_bonded = QtWidgets.QCheckBox(self.tab_10)
        self.remove_single_bonded.setChecked(True)
        self.remove_single_bonded.setObjectName("remove_single_bonded")
        self.gridLayout_21.addWidget(self.remove_single_bonded, 0, 0, 1, 1)
        self.remove_selected = QtWidgets.QCheckBox(self.tab_10)
        self.remove_selected.setObjectName("remove_selected")
        self.gridLayout_21.addWidget(self.remove_selected, 1, 0, 1, 1)
        self.select_atoms_removal = QtWidgets.QPushButton(self.tab_10)
        self.select_atoms_removal.setObjectName("select_atoms_removal")
        self.gridLayout_21.addWidget(self.select_atoms_removal, 2, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout_21)
        self.tabWidget.addTab(self.tab_10, "")
        self.gridLayout_7.addWidget(self.tabWidget, 1, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1017, 20))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget_2.setCurrentIndex(0)
        self.tabWidget_3.setCurrentIndex(0)
        self.bands_color.setCurrentIndex(0)
        self.tabWidget_4.setCurrentIndex(0)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "1D systems"))
        self.swave.setText(_translate("MainWindow", "0.0"))
        self.Bz.setText(_translate("MainWindow", "0.0"))
        self.Bx.setText(_translate("MainWindow", "0.0"))
        self.label_5.setText(_translate("MainWindow", "Zeeman Jz"))
        self.mAB.setText(_translate("MainWindow", "0.0"))
        self.By.setText(_translate("MainWindow", "0.0"))
        self.label.setText(_translate("MainWindow", "Fermi energy"))
        self.peierls.setText(_translate("MainWindow", "0.0"))
        self.kanemele.setText(_translate("MainWindow", "0.0"))
        self.label_24.setText(_translate("MainWindow", "Haldane"))
        self.antihaldane.setText(_translate("MainWindow", "0.0"))
        self.label_26.setText(_translate("MainWindow", "swave pairing"))
        self.fermi.setText(_translate("MainWindow", "0.0"))
        self.rashba.setText(_translate("MainWindow", "0.0"))
        self.haldane.setText(_translate("MainWindow", "0.0"))
        self.mAF.setText(_translate("MainWindow", "0.0"))
        self.label_10.setText(_translate("MainWindow", "Rashba"))
        self.label_11.setText(_translate("MainWindow", "Kane-Mele"))
        self.label_3.setText(_translate("MainWindow", "Zeeman Jx"))
        self.label_4.setText(_translate("MainWindow", "Zeeman Jy"))
        self.label_12.setText(_translate("MainWindow", "Sublattice imbalance"))
        self.label_14.setText(_translate("MainWindow", "Magnetic field"))
        self.label_13.setText(_translate("MainWindow", "Antiferromagnetism"))
        self.label_25.setText(_translate("MainWindow", "Anti-Haldane"))
        self.antikanemele.setText(_translate("MainWindow", "0.0"))
        self.label_37.setText(_translate("MainWindow", "Anti Kane-Mele"))
        self.crystalfield.setText(_translate("MainWindow", "0.0"))
        self.label_17.setText(_translate("MainWindow", "Crystal field"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_3), _translate("MainWindow", "Terms in the Hamiltonian"))
        self.show_structure.setText(_translate("MainWindow", "Show structure"))
        self.show_structure_3d.setText(_translate("MainWindow", "Show structure 3D"))
        self.nsuper_struct.setText(_translate("MainWindow", "5"))
        self.label_7.setText(_translate("MainWindow", "Supercell"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_4), _translate("MainWindow", "Structure"))
        self.show_bands.setText(_translate("MainWindow", "Band structure"))
        self.label_9.setText(_translate("MainWindow", "# kpoints"))
        self.label_15.setText(_translate("MainWindow", "Operator"))
        self.bands_color.setItemText(0, _translate("MainWindow", "None"))
        self.bands_color.setItemText(1, _translate("MainWindow", "y-position"))
        self.bands_color.setItemText(2, _translate("MainWindow", "Sx"))
        self.bands_color.setItemText(3, _translate("MainWindow", "Sy"))
        self.bands_color.setItemText(4, _translate("MainWindow", "Sz"))
        self.bands_color.setItemText(5, _translate("MainWindow", "Valley"))
        self.bands_color.setItemText(6, _translate("MainWindow", "IPR"))
        self.bands_color.setItemText(7, _translate("MainWindow", "Bulk"))
        self.bands_color.setItemText(8, _translate("MainWindow", "Edge"))
        self.nk_bands.setText(_translate("MainWindow", "100"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_5), _translate("MainWindow", "Bands"))
        self.delta_kbands.setText(_translate("MainWindow", "0.02"))
        self.label_27.setText(_translate("MainWindow", "Smearing"))
        self.ne_kbands.setText(_translate("MainWindow", "400"))
        self.label_28.setText(_translate("MainWindow", "# of energies"))
        self.label_29.setText(_translate("MainWindow", "Energy window"))
        self.window_kbands.setText(_translate("MainWindow", "3.0"))
        self.label_30.setText(_translate("MainWindow", "KPM scale"))
        self.scale_kbands.setText(_translate("MainWindow", "10.0"))
        self.label_31.setText(_translate("MainWindow", "# vectors"))
        self.nv_kbands.setText(_translate("MainWindow", "3"))
        self.show_dosbands.setToolTip(_translate("MainWindow", "This is equivalent to band structure calculation, but it can be applied for very large systems"))
        self.show_dosbands.setText(_translate("MainWindow", "Show DOS Bands"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_9), _translate("MainWindow", "DOS Bands"))
        self.label_16.setText(_translate("MainWindow", "Smearing"))
        self.dos_delta.setText(_translate("MainWindow", "0.01"))
        self.label_40.setText(_translate("MainWindow", "Number of k-points"))
        self.dos_nk.setText(_translate("MainWindow", "400"))
        self.dos_ewindow.setText(_translate("MainWindow", "4.0"))
        self.label_41.setText(_translate("MainWindow", "Energy window"))
        self.show_dos.setText(_translate("MainWindow", "Density of states"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_6), _translate("MainWindow", "DOS"))
        self.multildos_ewindow.setToolTip(_translate("MainWindow", "Energy window"))
        self.multildos_ewindow.setText(_translate("MainWindow", "1.5"))
        self.multildos_delta.setToolTip(_translate("MainWindow", "Energy smearing"))
        self.multildos_delta.setText(_translate("MainWindow", "0.03"))
        self.label_36.setText(_translate("MainWindow", "Smearing"))
        self.label_42.setText(_translate("MainWindow", "Energy window"))
        self.multildos_nk.setToolTip(_translate("MainWindow", "Number of kpoints used"))
        self.multildos_nk.setText(_translate("MainWindow", "50"))
        self.show_multildos.setText(_translate("MainWindow", "Show LDOS"))
        self.label_43.setText(_translate("MainWindow", "Number of kpoints"))
        self.label_44.setText(_translate("MainWindow", "Number of unit cells"))
        self.multildos_nrep.setToolTip(_translate("MainWindow", "Number of replicas of the unit cell to plot"))
        self.multildos_nrep.setText(_translate("MainWindow", "15"))
        self.basis_ldos.setToolTip(_translate("MainWindow", "Choose the basis of the LDOS, projection onto the tight binding basis, or directly in true real space assuming a certain atomic-like wavefunction"))
        self.basis_ldos.setItemText(0, _translate("MainWindow", "Tight binding"))
        self.basis_ldos.setItemText(1, _translate("MainWindow", "Real space atomic orbitals"))
        self.ratomic_ldos.setToolTip(_translate("MainWindow", "Radii of the atomic-like wavefunctions put in every site. Only affects the result for ht eBasis \"Real space atomic orbitals\""))
        self.ratomic_ldos.setText(_translate("MainWindow", "1.5"))
        self.label_19.setText(_translate("MainWindow", "Basis"))
        self.label_20.setText(_translate("MainWindow", "Local orbital radii"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_7), _translate("MainWindow", "LDOS"))
        self.ldos_delta.setToolTip(_translate("MainWindow", "Energy smearing"))
        self.ldos_delta.setText(_translate("MainWindow", "0.03"))
        self.label_38.setText(_translate("MainWindow", "Smearing"))
        self.label_45.setText(_translate("MainWindow", "Energy window"))
        self.ldos_nk.setToolTip(_translate("MainWindow", "Number of kpoints used"))
        self.ldos_nk.setText(_translate("MainWindow", "50"))
        self.show_ldos.setText(_translate("MainWindow", "Show LDOS"))
        self.ldos_ewindow.setToolTip(_translate("MainWindow", "Energy window"))
        self.ldos_ewindow.setText(_translate("MainWindow", "1.5"))
        self.label_46.setText(_translate("MainWindow", "Number of kpoints"))
        self.ldos_operator.setItemText(0, _translate("MainWindow", "None"))
        self.ldos_operator.setItemText(1, _translate("MainWindow", "current"))
        self.label_18.setText(_translate("MainWindow", "Operator"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_13), _translate("MainWindow", "E-y map"))
        self.scf_initialization.setItemText(0, _translate("MainWindow", "antiferro"))
        self.scf_initialization.setItemText(1, _translate("MainWindow", "ferro"))
        self.scf_initialization.setItemText(2, _translate("MainWindow", "random"))
        self.label_22.setText(_translate("MainWindow", "Initialization"))
        self.label_23.setText(_translate("MainWindow", "Hubbard"))
        self.hubbard.setText(_translate("MainWindow", "2.0"))
        self.label_34.setText(_translate("MainWindow", "Filling"))
        self.filling_scf.setText(_translate("MainWindow", "0.5"))
        self.do_scf.setText(_translate("MainWindow", "Include mean field"))
        self.solve_scf.setText(_translate("MainWindow", "Solve SCF"))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_12), _translate("MainWindow", "Basic"))
        self.label_32.setText(_translate("MainWindow", "Mixing"))
        self.label_33.setText(_translate("MainWindow", "# of kpoints"))
        self.nk_scf.setText(_translate("MainWindow", "10"))
        self.mix_scf.setText(_translate("MainWindow", "0.9"))
        self.label_35.setText(_translate("MainWindow", "Smearing"))
        self.smearing_scf.setText(_translate("MainWindow", "0.01"))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_11), _translate("MainWindow", "Convergence"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab), _translate("MainWindow", "SCF"))
        self.label_39.setText(_translate("MainWindow", "Number of unit cells"))
        self.magnetization_nrep.setText(_translate("MainWindow", "5"))
        self.show_magnetism.setText(_translate("MainWindow", "Show magnetism"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_8), _translate("MainWindow", "Magnetism"))
        self.label_8.setText(_translate("MainWindow", "Type of lattice"))
        self.lattice.setItemText(0, _translate("MainWindow", "Honeycomb zigzag"))
        self.lattice.setItemText(1, _translate("MainWindow", "Honeycomb armchair"))
        self.lattice.setItemText(2, _translate("MainWindow", "Square"))
        self.lattice.setItemText(3, _translate("MainWindow", "Chain"))
        self.lattice.setItemText(4, _translate("MainWindow", "Monochain"))
        self.lattice.setItemText(5, _translate("MainWindow", "Triangular"))
        self.lattice.setItemText(6, _translate("MainWindow", "Kagome"))
        self.lattice.setItemText(7, _translate("MainWindow", "Lieb"))
        self.nsuper.setText(_translate("MainWindow", "1"))
        self.label_6.setText(_translate("MainWindow", "Supercell"))
        self.label_2.setText(_translate("MainWindow", "Width"))
        self.width.setText(_translate("MainWindow", "10"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Geometry"))
        self.remove_single_bonded.setToolTip(_translate("MainWindow", "Remove atoms that have a single bond in the structure"))
        self.remove_single_bonded.setText(_translate("MainWindow", "Remove single bonds"))
        self.remove_selected.setText(_translate("MainWindow", "Remove selected atoms"))
        self.select_atoms_removal.setText(_translate("MainWindow", "Select atoms to remove"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_10), _translate("MainWindow", "Modify geometry"))

