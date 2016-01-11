import os
import unittest
import random
import math
import tempfile
import time
import numpy
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *

#
# TargetRegistrationError
#

class TargetRegistrationError(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "TargetRegistrationError" # TODO make this more human readable by adding spaces
    self.parent.categories = ["IGT"]
    self.parent.dependencies = []
    self.parent.contributors = ["Junichi Tokuda (Brigham and Women's Hospital)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This module calculates target registration error (TRE) based on two MarkUps nodes.
    """
    self.parent.acknowledgementText = """
    This module was developed using a template by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# TargetRegistrationErrorWidget
#

class TargetRegistrationErrorWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    # Instantiate and connect widgets ...

    #--------------------------------------------------
    # For debugging
    #
    # Reload and Test area
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    reloadCollapsibleButton.collapsed = True
    
    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "NeedleGuideTemlpate Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)
    #
    #--------------------------------------------------


    #
    # Fiducial Node
    #
    fiducialsCollapsibleButton = ctk.ctkCollapsibleButton()
    fiducialsCollapsibleButton.text = "Fiducials"
    self.layout.addWidget(fiducialsCollapsibleButton)

    # Layout within the dummy collapsible button
    fiducialsFormLayout = qt.QFormLayout(fiducialsCollapsibleButton)

    #
    # Fiducial node selector
    #
    self.fiducial1Selector = slicer.qMRMLNodeComboBox()
    self.fiducial1Selector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.fiducial1Selector.selectNodeUponCreation = True
    self.fiducial1Selector.addEnabled = True
    self.fiducial1Selector.removeEnabled = True
    self.fiducial1Selector.noneEnabled = True
    self.fiducial1Selector.renameEnabled = True
    self.fiducial1Selector.showHidden = False
    self.fiducial1Selector.showChildNodeTypes = False
    self.fiducial1Selector.setMRMLScene( slicer.mrmlScene )
    self.fiducial1Selector.setToolTip( "Pick the input to the algorithm." )
    fiducialsFormLayout.addRow("Fiducial 1: ", self.fiducial1Selector)

    self.fiducial2Selector = slicer.qMRMLNodeComboBox()
    self.fiducial2Selector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.fiducial2Selector.selectNodeUponCreation = True
    self.fiducial2Selector.addEnabled = True
    self.fiducial2Selector.removeEnabled = True
    self.fiducial2Selector.noneEnabled = True
    self.fiducial2Selector.renameEnabled = True
    self.fiducial2Selector.showHidden = False
    self.fiducial2Selector.showChildNodeTypes = False
    self.fiducial2Selector.setMRMLScene( slicer.mrmlScene )
    self.fiducial2Selector.setToolTip( "Pick the input to the algorithm." )
    fiducialsFormLayout.addRow("Fiducial 2: ", self.fiducial2Selector)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    fiducialsFormLayout.addRow(self.applyButton)

    #
    # Result
    #
    self.TREEdit = qt.QLineEdit()
    self.TREEdit.text = '--'
    self.TREEdit.readOnly = True
    self.TREEdit.frame = True
    self.TREEdit.styleSheet = "QLineEdit { background:transparent; }"
    self.TREEdit.cursor = qt.QCursor(qt.Qt.IBeamCursor)
    
    fiducialsFormLayout.addRow("TRE (mm):", self.TREEdit)

    self.errorRAEdit = qt.QLineEdit()
    self.errorRAEdit.text = '--'
    self.errorRAEdit.readOnly = True
    self.errorRAEdit.frame = True
    self.errorRAEdit.styleSheet = "QLineEdit { background:transparent; }"
    self.errorRAEdit.cursor = qt.QCursor(qt.Qt.IBeamCursor)
    
    fiducialsFormLayout.addRow("Error RA (mm):", self.errorRAEdit)

    self.errorAPEdit = qt.QLineEdit()
    self.errorAPEdit.text = '--'
    self.errorAPEdit.readOnly = True
    self.errorAPEdit.frame = True
    self.errorAPEdit.styleSheet = "QLineEdit { background:transparent; }"
    self.errorAPEdit.cursor = qt.QCursor(qt.Qt.IBeamCursor)
    
    fiducialsFormLayout.addRow("Error AP (mm):", self.errorAPEdit)

    self.errorSIEdit = qt.QLineEdit()
    self.errorSIEdit.text = '--'
    self.errorSIEdit.readOnly = True
    self.errorSIEdit.frame = True
    self.errorSIEdit.styleSheet = "QLineEdit { background:transparent; }"
    self.errorSIEdit.cursor = qt.QCursor(qt.Qt.IBeamCursor)
    
    fiducialsFormLayout.addRow("Error SI (mm):", self.errorSIEdit)
    
    
    # connections
    self.fiducial1Selector.connect("currentNodeChanged(vtkMRMLNode*)", self.onFiducialSelect)
    self.fiducial2Selector.connect("currentNodeChanged(vtkMRMLNode*)", self.onFiducialSelect)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Create logic
    self.logic = TargetRegistrationErrorLogic(None)

    # Enable buttons, if nodes are selected
    self.onFiducialSelect()

  def cleanup(self):
    pass

  def onFiducialSelect(self):
    #self.reconfigureButton.enabled = self.fiducial1Selector.currentNode()
    if self.fiducial1Selector.currentNode() and self.fiducial2Selector.currentNode():
      self.applyButton.enabled = True
    else:
      self.applyButton.enabled = False

  def onApplyButton(self):
    if self.logic.run(self.fiducial1Selector.currentNode(), self.fiducial2Selector.currentNode()):
      self.TREEdit.text = '%.3f' % self.logic.TRE
      self.errorRAEdit.text = '%.3f +/- %.3f' % (self.logic.meanError[0], self.logic.stdError[0])
      self.errorAPEdit.text = '%.3f +/- %.3f' % (self.logic.meanError[1], self.logic.stdError[1])
      self.errorSIEdit.text = '%.3f +/- %.3f' % (self.logic.meanError[2], self.logic.stdError[2])
    else:
      self.TREEdit.text = 'ERROR'
      self.errorRAEdit.text = 'ERROR'
      self.errorAPEdit.text = 'ERROR'
      self.errorSIEdit.text = 'ERROR'

  def onReload(self, moduleName="TargetRegistrationError"):
    # Generic reload method for any scripted module.
    # ModuleWizard will subsitute correct default moduleName.

    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)

#
# TargetRegistrationErrorLogic
#

class TargetRegistrationErrorLogic(ScriptedLoadableModuleLogic):

  def __init__(self, parent):
    ScriptedLoadableModuleLogic.__init__(self, parent)
    self.logFilePath = ''
    self.logFile = None

    self.TRE = 0.0
    self.meanError = [0.0, 0.0, 0.0]
    self.stdError = [0.0, 0.0, 0.0]
      
  def run(self, fiducial1Node, fiducial2Node):
    """
    Run the actual algorithm
    """
    
    if fiducial1Node == None or fiducial2Node == None:
      print "Fiducial node is not specified."
      return False
    
    nFid1 = fiducial1Node.GetNumberOfFiducials()
    nFid2 = fiducial2Node.GetNumberOfFiducials()

    if nFid1 != nFid2:
      print "The numbers of fiducial nodes do not match."
      return False

    sqrsum = 0.0

    arrayX = numpy.array([])
    arrayY = numpy.array([])
    arrayZ = numpy.array([])
    
    for m in range(0, nFid1):
      pos1 = [0.0, 0.0, 0.0]
      pos2 = [0.0, 0.0, 0.0]
      fiducial1Node.GetNthFiducialPosition(m, pos1)
      fiducial2Node.GetNthFiducialPosition(m, pos2)
      #lb = baseFiducial.GetNthFiducialLabel(m)

      p1 = numpy.array(pos1)
      p2 = numpy.array(pos2)
      d = p2-p1
      
      sqrsum = sqrsum + numpy.inner(d, d)
      arrayX = numpy.append(arrayX, d[0])
      arrayY = numpy.append(arrayY, d[1])
      arrayZ = numpy.append(arrayZ, d[2])
      
    self.TRE = numpy.sqrt(sqrsum / nFid1)
    self.meanError = [numpy.mean(arrayX), numpy.mean(arrayY), numpy.mean(arrayZ)]
    self.stdError = [numpy.std(arrayX), numpy.std(arrayY), numpy.std(arrayZ)]
    
    return True


