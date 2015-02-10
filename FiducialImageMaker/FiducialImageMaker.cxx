/*=========================================================================

  Program:   Fiducial Image Maker CLI for 3D Slicer
  Module:    itkRenderSpatialObjectImageFilter.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  Copyright (c) Brigham and Women's Hosptial. All rights reserved.
  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "FiducialImageMakerCLP.h"


#include "itkPluginUtilities.h"

#include "itkRenderSpatialObjectImageFilter.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef    T InputPixelType;
  typedef    T OutputPixelType;

  typedef itk::Image<InputPixelType,  3> InputImageType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;

  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  /*
  typedef itk::SmoothingRecursiveGaussianImageFilter<
    InputImageType, OutputImageType>  FilterType;
  */
  typedef itk::RenderSpatialObjectImageFilter<
    InputImageType, OutputImageType>  FilterType;    

  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName( inputVolume.c_str() );

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  int n = marker.size();
  for (int i = 0; i < n; i ++)
    {
    double lps[3];
    lps[0] = -marker[i][0];
    lps[1] = -marker[i][1];
    lps[2] = marker[i][2];
    filter->AddFiducialCenter(lps);
    }
  
  filter->SetFiducialRadius(radius);
  filter->SetDefaultVoxelValue(defaultVoxelValue);
  filter->SetToleranceVolume(toleranceVolume);

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression(1);
  writer->Update();

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
