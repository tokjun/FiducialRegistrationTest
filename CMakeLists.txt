cmake_minimum_required(VERSION 2.8.9)

project(FiducialRegistrationTest)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.example.com/Slicer/Extensions/FiducialRegistrationTest")
set(EXTENSION_CATEGORY "IGT")
set(EXTENSION_CONTRIBUTORS "Junichi Tokuda (Brigham and Women's Hospital)")
set(EXTENSION_DESCRIPTION "Module to test fiducial registration algorithms")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/FiducialRegistrationTest.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/FiducialRegistrationTest/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(FiducialRegistrationTest)
add_subdirectory(FiducialImageMaker)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
