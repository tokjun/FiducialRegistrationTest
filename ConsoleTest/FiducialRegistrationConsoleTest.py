import os
import os.path
import unittest
import random
import math
import tempfile
import time
import numpy
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *


### Parameters
nFiducials = 6
noiseLevel = 0.0

workingDir = "/Users/junichi/Experiments/FiducialTest/Test-M%02d-N%1.2f" %(nFiducials, noiseLevel)
if not os.path.exists(workingDir): os.makedirs(workingDir)

lt = time.localtime()
logFileName = "log-%04d-%02d-%02d-%02d-%02d-%02d.txt" % (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)

nTrialsPerCondition = 5

# Fiducial and volume parameters
radius = 50
imageFOV = 300
pixelSpacing = 1.0
thicknessStep = 1.0
nThicknessSteps = 4

# Range for random transform
xRange = [-50.0, 50.0]
yRange = [-50.0, 50.0]
zRange = [-50.0, 50.0]

### Setup modules
slicer.util.selectModule('FiducialRegistrationTest')
testLogic = slicer.modules.FiducialRegistrationTestWidget.logic
imageMakerCLI = slicer.modules.imagemaker

testLogic.logFilePath = workingDir + '/' + logFileName
testLogic.logFile = open(testLogic.logFilePath, 'a')

testLogic.printLog("Console Test > zRot, trial, thickness \n")

### Volume node for template volume  (Volume that represents the size/resolution for marker images)
templateVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
slicer.mrmlScene.AddNode(templateVolumeNode)
templateVolumeNode.SetName("Template Volume")

### Generate a fiducial model
modelFiducialNodeName = "Model-Fiducial-%d-%d" % (radius, nFiducials)
modelFiducialNode = None

# If the file exists, load it. Otherwise, we generate one.
if os.path.isfile(workingDir+'/'+modelFiducialNodeName+'.fcsv'):
    (r, modelFiducialNode) = slicer.util.loadMarkupsFiducialList(workingDir+'/'+modelFiducialNodeName+'.fcsv', True)
else:
    modelFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")
    slicer.mrmlScene.AddNode(modelFiducialNode)
    modelFiducialNode.SetName(modelFiducialNodeName)
    testLogic.configFiducialModel(modelFiducialNode, radius, nFiducials, 20.0)
    slicer.util.saveNode(modelFiducialNode, workingDir+'/'+modelFiducialNodeName+'.fcsv')

### Generate transform
referenceMatrix = vtk.vtkMatrix4x4()


# test rotation about z-axis

for zRot in numpy.arange(0.0,numpy.pi/2.0+0.001, numpy.pi/36.0): # every 5 degrees
    
    for trial in range(0, nTrialsPerCondition):
        x = random.uniform(-pixelSpacing, pixelSpacing) 
        y = random.uniform(-pixelSpacing, pixelSpacing) 
        z = random.uniform(-thicknessStep*nThicknessSteps, thicknessStep*nThicknessSteps) 

        randomTransformName = "TestRandomTransform-zRot-%f-%03d" % (zRot, trial)
        randomTransform = None
        if os.path.isfile(workingDir+'/'+randomTransformName+'.h5'):
            (r, randomTransform) = slicer.util.loadTransform(workingDir+'/'+randomTransformName+'.h5', True)
            randomTransform.GetMatrixTransformToParent(referenceMatrix)
        else:
            testLogic.generateTransform([0.0,0.0,1.0], zRot, [x,y,z], referenceMatrix)
            testLogic.generateRandomTransform(xRange, yRange, zRange, referenceMatrix)
            randomTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
            randomTransform.SetMatrixTransformToParent(referenceMatrix)
            slicer.mrmlScene.AddNode(randomTransform)
            randomTransform.SetName(randomTransformName)

        testFiducialNodeName = "TestFiducial-zRot-%f-%03d" % (zRot, trial)
        testFiducialNode = None
        if os.path.isfile(workingDir+'/'+testFiducialNodeName+'.h5'):
            (r, testFiducialNode) = slicer.util.loadMarkupsFiducialList(workingDir+'/'+testFiducialNodeName+'.fcsv', True)
        else:
            testFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")  
            slicer.mrmlScene.AddNode(testFiducialNode)
            testFiducialNode.SetName(testFiducialNodeName)
            testFiducialNode.RemoveAllMarkups()
            nFid = modelFiducialNode.GetNumberOfFiducials()
            for m in range(0, nFid):
                pos = [0.0, 0.0, 0.0]
                modelFiducialNode.GetNthFiducialPosition(m, pos)
                lb = modelFiducialNode.GetNthFiducialLabel(m)
                testFiducialNode.AddFiducialFromArray(pos, lb)

            testFiducialNode.ApplyTransformMatrix(referenceMatrix)
            slicer.util.saveNode(testFiducialNode, workingDir+'/'+testFiducialNodeName+'.fcsv')
            slicer.util.saveNode(randomTransform, workingDir+'/'+randomTransformName+'.h5')

        #testLogic.printLog("Condition: zRot=%f trial=%d\n" % (zRot, trial))
        #testLogic.printMatrixInLine("Reference:", referenceMatrix)

        for thickness in numpy.arange(1.0, thicknessStep*nThicknessSteps+0.001, thicknessStep):

            testFiducialVolumeNodeName = "TestImage-zRot-%f-%03d-%f" % (zRot, trial, thickness)
            testFiducialVolumeNode = None
            
            #testLogic.printLog("thickness: %f\n" % thickness)

            if os.path.isfile(workingDir+'/'+testFiducialVolumeNodeName+'.h5'):
                (r, testFiducialVolumeNode) = slicer.util.loadVolume(workingDir+'/'+testFiducialVolumeNodeName+'.nrrd', {}, True )
            else:
                # Create template volume using the ImageMaker module
                imageMakerParameters = {}
                imageMakerParameters["OutputVolume"] = templateVolumeNode.GetID()
                imageMakerParameters["ScalarType"] = "unsigned_short"
                imageMakerParameters["NumberOfComponents"] = 1
                imageMakerParameters["Dimension"] = 3
                imageMakerParameters["Size"] = [int(imageFOV/pixelSpacing), int(imageFOV/pixelSpacing), int(imageFOV/thickness)]
                imageMakerParameters["Origin"] = [-imageFOV/2.0, -imageFOV/2.0, -imageFOV/2.0]
                imageMakerParameters["Spacing"] = [pixelSpacing, pixelSpacing, thickness]
                imageMakerParameters["Direction"] = [1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00]
                
                slicer.cli.run(imageMakerCLI, None, imageMakerParameters, True)

                testFiducialVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
                slicer.mrmlScene.AddNode(testFiducialVolumeNode)
                testLogic.generateFiducialImage(templateVolumeNode, testFiducialVolumeNode, testFiducialNode)
                testFiducialVolumeNode.SetName(testFiducialVolumeNodeName)
                slicer.util.saveNode(testFiducialVolumeNode, workingDir+'/'+testFiducialVolumeNode.GetName()+'.nrrd')

            testLogic.printLog("Console Test > %f, %d, %f \n" % (zRot, trial, thickness))
            testLogic.runRegistration(modelFiducialNode, testFiducialVolumeNode, testFiducialNode)
            slicer.mrmlScene.RemoveNode(testFiducialVolumeNode)
        
        slicer.mrmlScene.RemoveNode(randomTransform)
        slicer.mrmlScene.RemoveNode(testFiducialNode)
        slicer.mrmlScene.RemoveNode(testFiducialVolumeNode)

if testLogic.logFile:
    testLogic.logFile.close()


slicer.mrmlScene.RemoveNode(templateVolumeNode)
