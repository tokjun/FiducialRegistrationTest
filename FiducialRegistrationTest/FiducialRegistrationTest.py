import os
import unittest
import random
import math
import tempfile
import time
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *

#
# FiducialRegistrationTest
#

class FiducialRegistrationTest(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "FiducialRegistrationTest" # TODO make this more human readable by adding spaces
    self.parent.categories = ["IGT"]
    self.parent.dependencies = []
    self.parent.contributors = ["Junichi Tokuda (Brigham and Women's Hospital)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This module tests image-based fiducial detection and registration using synthetic images. 
    """
    self.parent.acknowledgementText = """
    This module was developed using a template by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# FiducialRegistrationTestWidget
#

class FiducialRegistrationTestWidget(ScriptedLoadableModuleWidget):
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
    self.fiducialSelector = slicer.qMRMLNodeComboBox()
    self.fiducialSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.fiducialSelector.selectNodeUponCreation = True
    self.fiducialSelector.addEnabled = True
    self.fiducialSelector.removeEnabled = True
    self.fiducialSelector.noneEnabled = False
    self.fiducialSelector.renameEnabled = True
    self.fiducialSelector.showHidden = False
    self.fiducialSelector.showChildNodeTypes = False
    self.fiducialSelector.setMRMLScene( slicer.mrmlScene )
    self.fiducialSelector.setToolTip( "Pick the input to the algorithm." )
    fiducialsFormLayout.addRow("Input Volume: ", self.fiducialSelector)


    #
    # Reconfigure Button
    #
    self.radiusEdit = qt.QDoubleSpinBox()
    self.radiusEdit.setMinimum(0.0)
    self.radiusEdit.setMaximum(500.0)
    self.radiusEdit.setSingleStep(0.5)
    self.radiusEdit.setValue(50)

    self.numFiducialsEdit = qt.QSpinBox()
    self.numFiducialsEdit.setMinimum(0)
    self.numFiducialsEdit.setMaximum(100)
    self.numFiducialsEdit.setSingleStep(1)
    self.numFiducialsEdit.setValue(5)

    fiducialsFormLayout.addRow("Radius (mm):", self.radiusEdit)
    fiducialsFormLayout.addRow("# of fiducials:", self.numFiducialsEdit)

    self.reconfigureButton = qt.QPushButton("Reconfigure Fiducials")
    self.reconfigureButton.toolTip = "Reconfigure fiducial frame"
    self.reconfigureButton.enabled = False
    fiducialsFormLayout.addRow(self.reconfigureButton)


    #
    # Test Area
    #
    testCollapsibleButton = ctk.ctkCollapsibleButton()
    testCollapsibleButton.text = "Test"
    self.layout.addWidget(testCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(testCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Fiducial: ", self.inputSelector)

    #
    # reference volume selector
    #
    self.referenceSelector = slicer.qMRMLNodeComboBox()
    self.referenceSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.referenceSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 0 )
    self.referenceSelector.selectNodeUponCreation = False
    self.referenceSelector.addEnabled = True
    self.referenceSelector.removeEnabled = True
    self.referenceSelector.noneEnabled = False
    self.referenceSelector.showHidden = False
    self.referenceSelector.showChildNodeTypes = False
    self.referenceSelector.setMRMLScene( slicer.mrmlScene )
    self.referenceSelector.setToolTip( "Pick the reference to the algorithm." )
    parametersFormLayout.addRow("Reference Volume: ", self.referenceSelector)

    logFileLayout = qt.QHBoxLayout()
    self.logFileLineEdit = qt.QLineEdit()
    self.logFileLineEdit.text = ''
    self.logFileLineEdit.readOnly = True
    self.logFileLineEdit.frame = True
    self.logFileLineEdit.styleSheet = "QLineEdit { background:transparent; }"
    self.logFileLineEdit.cursor = qt.QCursor(qt.Qt.IBeamCursor)
    logFileLayout.addWidget(self.logFileLineEdit)

    self.logFileButton = qt.QPushButton("Choose File...")
    self.logFileButton.toolTip = "Choose log file from dialog box"
    logFileLayout.addWidget(self.logFileButton)

    parametersFormLayout.addRow("Log file:", logFileLayout)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.fiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onFiducialSelect)
    self.reconfigureButton.connect('clicked(bool)', self.onReconfigureButton)
    self.logFileButton.connect('clicked(bool)', self.onLogFileButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.referenceSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Create logic
    self.logic = FiducialRegistrationTestLogic(None)

    # Enable buttons, if nodes are selected
    self.onSelect()
    self.onFiducialSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector.currentNode() and self.referenceSelector.currentNode()

  def onFiducialSelect(self):
    self.reconfigureButton.enabled = self.fiducialSelector.currentNode()

  def onReconfigureButton(self):
    if self.fiducialSelector.currentNode():
      self.logic.configFiducialModel(self.fiducialSelector.currentNode(), self.radiusEdit.value, self.numFiducialsEdit.value, 20.0)

  def onLogFileButton(self):
 
    fileName = qt.QFileDialog.getSaveFileName(None, 'Open Log File', '', 'txt files (*.txt)')

    if fileName:
      self.logFileLineEdit.setText(fileName)

  def onApplyButton(self):
    self.logic.run(self.inputSelector.currentNode(), self.referenceSelector.currentNode(), self.logFileLineEdit.text)

  def onReload(self, moduleName="FiducialRegistrationTest"):
    # Generic reload method for any scripted module.
    # ModuleWizard will subsitute correct default moduleName.

    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)

#
# FiducialRegistrationTestLogic
#

class FiducialRegistrationTestLogic(ScriptedLoadableModuleLogic):

  def __init__(self, parent):
    ScriptedLoadableModuleLogic.__init__(self, parent)
    self.logFilePath = ''
    self.logFile = None


  def configFiducialModel(self, fiducialNode, radius, nFiducials, minDistance):

    # Create a circular fiducial with random spacing. (All spacings are unique.)
    name = fiducialNode.GetName()
    fiducialNode.RemoveAllMarkups()

    # The angle between the markers can be generated by shuffling
    #
    #     [ (m)*u, (m+1)*u, (m+2)*u, ..., (m+N-1)*u ]              (1)
    # 
    # where m is the minimum number of unites to keep the spacing greater than the minDistance (d),
    # u is the unit angle. The sum of the spacings must be 2*pi:
    #
    #     { m + (m+1) + (m+2) + ... + (m+N-1) } * u = 2*pi         (2)
    #     {N*m + (N-1) * N / 2} * u = 2*pi                         (3)
    #
    # To keep the spacing greater than the minDistance,
    #
    #     (m*u) * r > d                                            (4)
    #     u > d / (m * r)                                          (5)
    # 
    # From (5) and (3)
    #
    #     {N*m + (N-1) * N / 2} * d / (m*r) < 2*pi                 (6)
    #     m > (N-1)*N*d / (4*r*pi-2*N*d)                           (7)
    # 
    # m has to be an integer. Let the <m> = ceil(m) (smallest inegral value not less than x)
    # From (3),
    #
    #     u = 2*pi / {N*m + (N-1)*N/2}
    #
    # In the following code, we use 'unit' as u, 'minUnit' uas m,
    # 'nFiducials' as N, and 'radius' as r, 'minDistance' as d

    fnFiducials = float(nFiducials)
    minUnit = (fnFiducials-1.0)*fnFiducials*minDistance / (4.0*radius*math.pi - 2.0*fnFiducials*minDistance)
    minUnit = math.ceil(minUnit)
    unit = 2*math.pi / (fnFiducials*minUnit + (fnFiducials-1.0)*fnFiducials/2.0)

    spacing = range(int(minUnit), int(minUnit+nFiducials))
    random.shuffle(spacing)

    n = 0
    theta = 0.0
    pos = [0.0, 0.0, 0.0]
    for s in spacing:
      theta = theta + unit * s
      pos[0] = radius * math.sin(theta)
      pos[1] = radius * math.cos(theta)
      fiducialNode.AddFiducialFromArray(pos, "%s-%d" % (name, n))
      n = n + 1
      

  def generateRandomTransform(self, xRange, yRange, zRange, matrix):
    x = random.uniform(xRange[0], xRange[1]) 
    y = random.uniform(yRange[0], yRange[1]) 
    z = random.uniform(zRange[0], zRange[1]) 
    alpha = random.uniform(-math.pi,   math.pi) 
    beta  = random.uniform(-math.pi/2.0, math.pi/2.0)
    gamma = random.uniform(-math.pi,   math.pi) 

    s1 = math.sin(alpha)
    c1 = math.cos(alpha)
    s2 = math.sin(beta)
    c2 = math.cos(beta)
    s3 = math.sin(gamma)
    c3 = math.cos(gamma)

    # Rotation
    matrix.SetElement(0, 0, c2)
    matrix.SetElement(0, 1, -c3*s2)
    matrix.SetElement(0, 2, s2*s3)
    matrix.SetElement(1, 0, c1*s2)
    matrix.SetElement(1, 1, c1*c2*c3-s1*s3)
    matrix.SetElement(1, 2, -c3*s1-c1*c2*s3)
    matrix.SetElement(2, 0, s1*s2)
    matrix.SetElement(2, 1, c1*s3+c2*c3*s1)
    matrix.SetElement(2, 2, c1*c3-c2*s1*s3)

    # Translation
    matrix.SetElement(0, 3, x)
    matrix.SetElement(1, 3, y)
    matrix.SetElement(2, 3, z)

    matrix.SetElement(3, 0, 0.0)
    matrix.SetElement(3, 1, 0.0)
    matrix.SetElement(3, 2, 0.0)
    matrix.SetElement(3, 3, 1.0)


  def hasImageData(self,volumeNode):
    """This is a dummy logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      print('no volume node')
      return False
    if volumeNode.GetImageData() == None:
      print('no image data')
      return False
    return True

  
  def runRegistration(self, fiducialNode, volumeNode, testFiducialNode):

    logging.info('Processing started')

    # Get CLI modules
    fiducialDetectionCLI = slicer.modules.sphericalfiducialdetection
    circleFitCLI = slicer.modules.circlefit

    # Create temporary filename to store detected fiducials
    tmpImageFiducialFilename = tempfile.NamedTemporaryFile().name + ".fcsv"
    
    # Call fiducial detection
    detectionParameters = {}
    detectionParameters["inputVolume"] = volumeNode
    detectionParameters["outputFile"] = tmpImageFiducialFilename
    detectionParameters["threshold"] = 0.0
    detectionParameters["numberOfSpheres"] = fiducialNode.GetNumberOfFiducials()
    detectionParameters["sigmaGrad"] = 1.0
    detectionParameters["gradThreshold"] = 0.1
    detectionParameters["minRadius"] = 4.0
    detectionParameters["maxRadius"] = 6.0
    detectionParameters["variance"] = 1.0
    detectionParameters["outputThreshold"] = 0.5
    detectionParameters["sphereRadiusRatio"] = 1.0

    detectionParameters["alpha"] = 0.8
    detectionParameters["beta"] = 0.8
    detectionParameters["gamma"] = 0.8

    detectionParameters["minSigma"] = 3.0
    detectionParameters["maxSigma"] = 3.0
    detectionParameters["stepSigma"] = 1.0
    
    detectionParameters["debugSwitch"] = 0
    
    slicer.cli.run(fiducialDetectionCLI, None, detectionParameters, True)

    # Import fiducials in slicer scene
    (success, imageFiducialNode) = slicer.util.loadMarkupsFiducialList(tmpImageFiducialFilename, True)
    imageFiducialNode.SetName(slicer.mrmlScene.GenerateUniqueName('ImageFiducialsDetected'))
    
    
    #### Circle Fitting

    cfTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
    slicer.mrmlScene.AddNode(cfTransform)

    circleFitParameters = {}
    circleFitParameters["movingPoints"] = fiducialNode
    circleFitParameters["fixedPoints"] = imageFiducialNode
    circleFitParameters["radius"] = 50.0
    circleFitParameters["registration"] = cfTransform

    time0 = time.clock()
    slicer.cli.run(circleFitCLI, None, circleFitParameters, True)
    time1 = time.clock()

    matrix = vtk.vtkMatrix4x4()

    cfTransform.GetMatrixTransformToParent(matrix)
    self.printMatrixInLine("CircleFit", matrix)
    t = time1-time0
    self.printLog ("Time - CircleFit: %f\n" % t)

    self.evaluateRegistration(fiducialNode, testFiducialNode, matrix)
    
    slicer.mrmlScene.RemoveNode(cfTransform)

    #### ICP

    initialTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
    slicer.mrmlScene.AddNode(initialTransform)

    icpRegistrationCLI = slicer.modules.icpregistration

    icpTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
    slicer.mrmlScene.AddNode(icpTransform)

    icpRegistrationError = 0.0
    registrationParameters = {}
    registrationParameters["movingPoints"] = fiducialNode
    registrationParameters["fixedPoints"] = imageFiducialNode
    registrationParameters["initialTransform"] = initialTransform
    registrationParameters["registrationTransform"] = icpTransform
    
    registrationParameters["iterations"] = 2000
    registrationParameters["gradientTolerance"] = 0.0001
    registrationParameters["valueTolerance"] = 0.0001
    registrationParameters["epsilonFunction"] = 0.00001

    time0 = time.clock()
    cliNode = slicer.cli.run(icpRegistrationCLI, None, registrationParameters, True)
    time1 = time.clock()

    icpTransform.GetMatrixTransformToParent(matrix)
    self.printMatrixInLine("ICP", matrix)
    t = time1-time0
    self.printLog ("Time - ICP: %f\n" % t)

    self.evaluateRegistration(fiducialNode, testFiducialNode, matrix)

    # Cleanup
    slicer.mrmlScene.RemoveNode(icpTransform)
    

  def evaluateRegistration(self, fiducialNode, testFiducialNode, matrix):


    nFid = fiducialNode.GetNumberOfFiducials()
    for m in range(0, nFid):
      pos = [0.0, 0.0, 0.0]
      fiducialNode.GetNthFiducialPosition(m, pos)
      tpos
      matrix.MultiplyPoint(pos, tpos)



      testFiducialNode.AddFiducialFromArray(pos, lb)
    
      





  def generateFiducialImage(self, backgroundVolumeNode, outputVolumeNode, fiducialNode):

    cli = slicer.modules.fiducialimagemaker

    parameters = {}
    parameters["inputVolume"] = backgroundVolumeNode
    parameters["outputVolume"] = outputVolumeNode
    parameters["marker"] = fiducialNode
    parameters["radius"] = 5.0
    parameters["defaultVoxelValue"] = 200.0
    parameters["toleranceVolume"] = 0.01

    cliNode = slicer.cli.run(cli, None, parameters, True)

  def printMatrixInLine(self, name, matrix):
    if self.logFile:
      self.logFile.write("%s: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n" 
                         % (name, 
                            matrix.GetElement(0,0), matrix.GetElement(0,1), matrix.GetElement(0,2), matrix.GetElement(0,3),
                            matrix.GetElement(1,0), matrix.GetElement(1,1), matrix.GetElement(1,2), matrix.GetElement(1,3),
                            matrix.GetElement(2,0), matrix.GetElement(2,1), matrix.GetElement(2,2), matrix.GetElement(2,3),
                            matrix.GetElement(3,0), matrix.GetElement(3,1), matrix.GetElement(3,2), matrix.GetElement(3,3)))
  def printLog(self, text):
    if self.logFile:
      self.logFile.write(text)
    else:
      print text
      
  def run(self, baseFiducial, inputVolume, logPath=None):
    """
    Run the actual algorithm
    """

    if logPath:
      self.logFilePath = logPath
      self.logFile = open(logPath, 'a')

    srcMatrix = vtk.vtkMatrix4x4()
    
    xRange = [-50.0, 50.0]
    yRange = [-50.0, 50.0]
    zRange = [-50.0, 50.0]

    for n in range(1, 2):

      testFiducialVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
      slicer.mrmlScene.AddNode(testFiducialVolumeNode)
      testFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")
      slicer.mrmlScene.AddNode(testFiducialNode)

      #numMarkups = baseFiducial.GetNumberOfMarkups()
      #for m in range(0, numMarkups):
      #  markup = baseFiducial.GetNthMarkup(m)
      #  testFiducialNode.AddMarkup(markup)

      testFiducialNode.RemoveAllMarkups()
      nFid = baseFiducial.GetNumberOfFiducials()
      for m in range(0, nFid):
        pos = [0.0, 0.0, 0.0]
        baseFiducial.GetNthFiducialPosition(m, pos)
        lb = baseFiducial.GetNthFiducialLabel(m)
        testFiducialNode.AddFiducialFromArray(pos, lb)

      #testFiducialNode.Copy(baseFiducial)

      testFiducialNode.SetName("TestFiducial-%d" % n)
      testFiducialVolumeNode.SetName("TestImage-%d" % n)

      self.generateRandomTransform(xRange, yRange, zRange, srcMatrix)
      randomTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
      randomTransform.SetMatrixTransformToParent(srcMatrix)
      slicer.mrmlScene.AddNode(randomTransform)

      testFiducialNode.ApplyTransformMatrix(srcMatrix)

      self.generateFiducialImage(inputVolume, testFiducialVolumeNode, testFiducialNode)
      self.printMatrixInLine("src", srcMatrix)

      self.runRegistration(baseFiducial, testFiducialVolumeNode, testFiducialNode)


    if self.logFile:
      self.logFile.close()

    return True


