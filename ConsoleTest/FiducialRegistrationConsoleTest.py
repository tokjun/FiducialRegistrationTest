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


def generateModel(modelFiducialName, radius, nFiducials, workingDir):
    modelFiducialNode = None

    # If the file exists, load it. Otherwise, we generate one.
    if os.path.isfile(workingDir+'/'+modelFiducialName+'.fcsv'):
        (r, modelFiducialNode) = slicer.util.loadMarkupsFiducialList(workingDir+'/'+modelFiducialName+'.fcsv', True)
    else:
        modelFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")
        slicer.mrmlScene.AddNode(modelFiducialNode)
        modelFiducialNode.SetName(modelFiducialName)
        testLogic.configFiducialModel(modelFiducialNode, radius, nFiducials, 20.0)
        slicer.util.saveNode(modelFiducialNode, workingDir+'/'+modelFiducialName+'.fcsv')
        
    return modelFiducialNode


# Generate a random transform. If the transform already exists in the working directory,
# the function will load it to the scene.
def generateRandomTransform(randomMatrix, randomTransformName, workingDir):
    
    randomTransform = None
    
    if os.path.isfile(workingDir+'/'+randomTransformName+'.h5'):
        (r, randomTransform) = slicer.util.loadTransform(workingDir+'/'+randomTransformName+'.h5', True)
        randomTransform.GetMatrixTransformToParent(randomMatrix)
        slicer.mrmlScene.RemoveNode(randomTransform)
    else:
        testLogic.generateRandomTransform(xRange, yRange, zRange, randomMatrix)
        randomTransform = slicer.mrmlScene.CreateNodeByClass("vtkMRMLLinearTransformNode")
        randomTransform.SetMatrixTransformToParent(randomMatrix)
        slicer.mrmlScene.AddNode(randomTransform)
        randomTransform.SetName(randomTransformName)
        slicer.util.saveNode(randomTransform, workingDir+'/'+randomTransformName+'.h5')
        slicer.mrmlScene.RemoveNode(randomTransform)

def generateTestFiducial(modelFiducialNode, randomMatrix):

    testFiducialNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLMarkupsFiducialNode")  
    slicer.mrmlScene.AddNode(testFiducialNode)
    testFiducialNode.RemoveAllMarkups()
    nFid = modelFiducialNode.GetNumberOfFiducials()
    for m in range(0, nFid):
        pos = [0.0, 0.0, 0.0]
        modelFiducialNode.GetNthFiducialPosition(m, pos)
        lb = modelFiducialNode.GetNthFiducialLabel(m)
        testFiducialNode.AddFiducialFromArray(pos, lb)
            
    testFiducialNode.ApplyTransformMatrix(randomMatrix)
    
    return testFiducialNode
        

def generateTestVolume(testFiducialNode, imageFOV, pixelSpacing, thickness, workingDir):

    testVolumeNode = None
    
    if os.path.isfile(workingDir+'/'+testVolumeName+'.h5'):
        (r, testVolumeNode) = slicer.util.loadVolume(workingDir+'/'+testVolumeName+'.nrrd', {}, True )
    else:
        ### Volume node for template volume  (Volume that represents the size/resolution for marker images)
        templateVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
        slicer.mrmlScene.AddNode(templateVolumeNode)
        templateVolumeNode.SetName("Template Volume")

        # Create template volume using the ImageMaker module
        imageMakerParameters = {}
        imageMakerParameters["OutputVolume"] = templateVolumeNode.GetID()
        imageMakerParameters["ScalarType"] = "unsigned_short"
        imageMakerParameters["NumberOfComponents"] = 1
        imageMakerParameters["Dimension"] = 3
        imageMakerParameters["Size"] = [int(imageFOV[0]/pixelSpacing), int(imageFOV[1]/pixelSpacing), int(imageFOV[2]/thickness)]
        imageMakerParameters["Origin"] = [-imageFOV[0]/2.0, -imageFOV[1]/2.0, -imageFOV[2]/2.0]
        imageMakerParameters["Spacing"] = [pixelSpacing, pixelSpacing, thickness]
        imageMakerParameters["Direction"] = [1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00]
        
        slicer.cli.run(imageMakerCLI, None, imageMakerParameters, True)
        
        testVolumeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLScalarVolumeNode")
        slicer.mrmlScene.AddNode(testVolumeNode)
        testLogic.generateFiducialImage(templateVolumeNode, testVolumeNode, testFiducialNode)
        testVolumeNode.SetName(testVolumeName)
        slicer.util.saveNode(testVolumeNode, workingDir+'/'+testVolumeNode.GetName()+'.nrrd')
        
        slicer.mrmlScene.RemoveNode(templateVolumeNode)

    return testVolumeNode


def computeEstimatedTRE(referenceMatrix, resultMatrix, needleLength):

    tipOffset = [0.0, 0.0, -needleLength, 1.0]
    referenceTipPos = [0.0, 0.0, 0.0, 1.0]
    registeredTipPos= [0.0, 0.0, 0.0, 1.0]

    referenceMatrix.MultiplyPoint(tipOffset, referenceTipPos)
    resultMatrix.MultiplyPoint(tipOffset, registeredTipPos)

    npReferenceTipPos = numpy.array(referenceTipPos[0:3])
    npRegisteredTipPos = numpy.array(registeredTipPos[0:3])

    errorVector = npReferenceTipPos - npRegisteredTipPos

    return numpy.linalg.norm(errorVector)


### Parameters
nFiducials = 6
noiseLevel = 0.0

#workingDir = "/Users/junichi/Experiments/FiducialTest/Test-M%02d-N%1.2f" %(nFiducials, noiseLevel)
workingDir = "/home/develop/Dropbox/Experiments/Canon/FiducialTest/Test-M%02d-N%1.2f" %(nFiducials, noiseLevel)
if not os.path.exists(workingDir): os.makedirs(workingDir)

lt = time.localtime()
logFileName = "log-%04d-%02d-%02d-%02d-%02d-%02d.txt" % (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)
csvFileName = "result-%04d-%02d-%02d-%02d-%02d-%02d.csv" % (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)

nTrialsPerCondition = 5

# Fiducial and volume parameters
radius = 50
imageFOV = [300, 255, 150]
pixelSpacing = 1.0
thicknessStep = 1.0
nThicknessSteps = 4

# Range for random transform
xRange = [-50.0, 50.0]
yRange = [-50.0, 50.0]
zRange = [-20.0, 20.0]

### Setup modules
slicer.util.selectModule('FiducialRegistrationTest')
testLogic = slicer.modules.FiducialRegistrationTestWidget.logic
imageMakerCLI = slicer.modules.imagemaker

testLogic.logFilePath = workingDir + '/' + logFileName
testLogic.logFile = open(testLogic.logFilePath, 'a')

testLogic.printLog("Console Test > trial, thickness \n")

csvFilePath = workingDir + '/' + csvFileName
csvFile = open(csvFilePath, 'a')

csvFile.write('thickness, trial, fre, fle, tre, nFid')

### Generate a fiducial model
modelFiducialName = "Model-Fiducial-%d-%d" % (radius, nFiducials)
modelFiducialNode = generateModel(modelFiducialName, radius, nFiducials, workingDir)

### Generate transform
randomMatrix = vtk.vtkMatrix4x4()

for trial in range(0, nTrialsPerCondition):

    randomTransformName = "TestRandomTransform-%03d" % (trial)
    generateRandomTransform(randomMatrix, randomTransformName, workingDir)

    testFiducialNode = generateTestFiducial(modelFiducialNode, randomMatrix)
    
    for thickness in numpy.arange(1.0, thicknessStep*nThicknessSteps+0.001, thicknessStep):

        testVolumeName = "TestImage-thickness-%d-%03d" % (thickness, trial)
        testVolumeNode = generateTestVolume(testFiducialNode, imageFOV, pixelSpacing, thickness, workingDir)

        testLogic.printLog("Console Test > %d, %f \n" % (trial, thickness))

        resultMatrix = vtk.vtkMatrix4x4()
        (fre, fle, nFidDetected) = testLogic.runRegistration(modelFiducialNode, testVolumeNode, testFiducialNode, resultMatrix)

        tre = computeEstimatedTRE(randomMatrix, resultMatrix, 150)
        csvFile.write('%f, %d, %f, %f, %f, %d\n' % (thickness, trial, fre, fle, tre, nFidDetected))

        slicer.mrmlScene.RemoveNode(testVolumeNode)
        
    slicer.mrmlScene.RemoveNode(testFiducialNode)
        
slicer.mrmlScene.RemoveNode(modelFiducialNode)

if testLogic.logFile:
    testLogic.logFile.close()

if csvFile:
    csvFile.close()
