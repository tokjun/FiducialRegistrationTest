import os
import unittest
import random
import math
import tempfile
import time
import numpy
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *


### Parameters
workingDir = "/Users/junichi/Experiments/FiducialTest"

lt = time.localtime()
logFileName = "log-%04d-%02d-%02d-%02d-%02d-%02d.txt" % (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)

nTrialsPerCondition = 2

# Fiducial and volume parameters
radius = 50
nFiducials = 6
imageFOV = 300
pixelSpacing = 1.0
thicknessStep = 1.0
nThicknessSteps = 4

# Range for random transform
xRange = [-50.0, 50.0]
yRange = [-50.0, 50.0]
zRange = [-50.0, 50.0]

### Setup modules
testLogic = slicer.modules.FiducialRegistrationTestWidget.logic
imageMakerCLI = slicer.modules.imagemaker

testLogic.logFilePath = workingDir + '/' + logFileName
testLogic.logFile = open(testLogic.logFilePath, 'a')


### Volume node for template volume  (Volume that represents the size/resolution for marker images)
templateVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
slicer.mrmlScene.AddNode(templateVolumeNode)
templateVolumeNode.SetName("Template Volume")

### Generate a fiducial model
modelFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")
slicer.mrmlScene.AddNode(modelFiducialNode)
modelFiducialNode.SetName("Model-Fiducial-%d-%d" % (radius, nFiducials))
testLogic.configFiducialModel(modelFiducialNode, radius, nFiducials, 20.0)
slicer.util.saveNode(modelFiducialNode, workingDir+'/'+modelFiducialNode.GetName()+'.fcsv')


### Generate transform
referenceMatrix = vtk.vtkMatrix4x4()


# test rotation about z-axis

for zRot in numpy.arange(0.0,numpy.pi/2.0+0.001, numpy.pi/36.0): # every 5 degrees
    
    for trial in range(0, nTrialsPerCondition):
        x = random.uniform(-pixelSpacing, pixelSpacing) 
        y = random.uniform(-pixelSpacing, pixelSpacing) 
        z = random.uniform(-thicknessStep*nThicknessSteps, thicknessStep*nThicknessSteps) 
        
        testLogic.generateTransform([0.0,0.0,1.0], zRot, [x,y,z], referenceMatrix)
        
        testLogic.generateRandomTransform(xRange, yRange, zRange, referenceMatrix)
        randomTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
        randomTransform.SetMatrixTransformToParent(referenceMatrix)
        slicer.mrmlScene.AddNode(randomTransform)
        randomTransform.SetName("TestRandomTransform-zRot-%f-%03d" % (zRot, trial))

        testFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")  
        slicer.mrmlScene.AddNode(testFiducialNode)
        testFiducialNode.SetName("TestFiducial-zRot-%f-%03d" % (zRot, trial))
        testFiducialNode.RemoveAllMarkups()

        nFid = modelFiducialNode.GetNumberOfFiducials()
        for m in range(0, nFid):
            pos = [0.0, 0.0, 0.0]
            modelFiducialNode.GetNthFiducialPosition(m, pos)
            lb = modelFiducialNode.GetNthFiducialLabel(m)
            testFiducialNode.AddFiducialFromArray(pos, lb)

        testFiducialNode.ApplyTransformMatrix(referenceMatrix)
        testLogic.printLog("Condition: zRot=%f trial=%d\n" % (zRot, trial))
        testLogic.printMatrixInLine("Reference:", referenceMatrix)

        slicer.util.saveNode(testFiducialNode, workingDir+'/'+testFiducialNode.GetName()+'.fcsv')
        slicer.util.saveNode(randomTransform, workingDir+'/'+randomTransform.GetName()+'.h5')

        for thickness in numpy.arange(1.0, thicknessStep*nThicknessSteps+0.001, thicknessStep):

            testLogic.printLog("thickness: %f\n" % thickness)
            # Create template volume using the ImageMaker module
            imageMakerParameters = {}
            imageMakerParameters["OutputVolume"] = templateVolumeNode.GetID()
            imageMakerParameters["ScalarType"] = "unsigned_short"
            imageMakerParameters["NumberOfComponents"] = 1
            imageMakerParameters["Dimension"] = 3
            imageMakerParameters["Size"] = [imageFOV, imageFOV, imageFOV]
            imageMakerParameters["Origin"] = [-imageFOV/2.0, -imageFOV/2.0, -imageFOV/2.0]
            imageMakerParameters["Spacing"] = [pixelSpacing, pixelSpacing, thickness]
            imageMakerParameters["Direction"] = [1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00]
        
            slicer.cli.run(imageMakerCLI, None, imageMakerParameters, True)

            testFiducialVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
            slicer.mrmlScene.AddNode(testFiducialVolumeNode)
            testLogic.generateFiducialImage(templateVolumeNode, testFiducialVolumeNode, testFiducialNode)
            testFiducialVolumeNode.SetName("TestImage-zRot-%f-%03d-%f" % (zRot, trial, thickness))
            slicer.util.saveNode(testFiducialVolumeNode, workingDir+'/'+testFiducialVolumeNode.GetName()+'.nrrd')
        
            testLogic.runRegistration(modelFiducialNode, testFiducialVolumeNode, testFiducialNode)

            slicer.mrmlScene.RemoveNode(testFiducialVolumeNode)
        
        slicer.mrmlScene.RemoveNode(randomTransform)
        slicer.mrmlScene.RemoveNode(testFiducialNode)

if testLogic.logFile:
    testLogic.logFile.close()


slicer.mrmlScene.RemoveNode(templateVolumeNode)
