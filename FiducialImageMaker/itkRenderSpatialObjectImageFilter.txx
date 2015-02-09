/*=========================================================================

  Program:   LineMarkerRegistration CLI for 3D Slicer
  Module:    itkRenderSpatialObjectImageFilter.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  This code is based on vtkImageToImageFilter.h in ITK.

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRenderSpatialObjectImageFilter_txx
#define __itkRenderSpatialObjectImageFilter_txx

#include "itkRenderSpatialObjectImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkListSample.h"
#include "itkCovarianceSampleFilter.h"
//#include "itkCovarianceCalculator.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkCrossHelper.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < typename  TInput, typename TOutput  >
RenderSpatialObjectImageFilter< TInput, TOutput >
::RenderSpatialObjectImageFilter()
{
  this->m_Label = 1;
  this->m_Normal[0] = 0.0;
  this->m_Normal[1] = 0.0;
  this->m_Normal[2] = 1.0;
}


template < typename  TInput, typename TOutput  >
void 
RenderSpatialObjectImageFilter< TInput, TOutput >
::GenerateData()
{
  itkDebugMacro(<< "RenderSpatialObjectImageFilter generating data ");
  
  typename InputImageType::ConstPointer input = this->GetInput();
  typename OutputImageType::Pointer output = this->GetOutput();

  ImageRegionConstIterator<InputImageType> it;
  it = ImageRegionConstIterator<InputImageType>( input, input->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());

  // TODO: Should outout image be initialized?
  output->FillBuffer(static_cast<OutputPixelType>(0));
  
  typedef itk::Statistics::ListSample< VectorType > PointListType;
  typedef std::map< InputPixelType, PointListType::Pointer > PointListMapType;

  oit.GoToBegin();
  it.GoToBegin();

  // TODO: Get voxel size
  VectorType voxSize;
  //voxSize[0] =
  //voxSize[1] =
  //voxSize[2] =

  while (!it.IsAtEnd())
    {
    InputPixelType pix = it.Get();
    typename InputImageType::IndexType index = it.GetIndex();

    if (pix == this->m_Label)
      {
      found = 1;
      typename InputImageType::PointType point;
      input->TransformIndexToPhysicalPoint (index, point);

      double value = this->ComputeObjectVolumeInCube(point, voxsize);

      // TODO: How is the input image incorporated?
      oit.Set( static_cast<OutputPixelType>(value+it.Get()) );
      }
    ++it;
    ++oit;
    }
}

template < typename  TInput, typename TOutput  >
void
RenderSpatialObjectImageFilter< TInput, TOutput >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


template < typename  TInput, typename TOutput  >
void
RenderSpatialObjectImageFilter< TInput, TOutput >
::ComputeObjectVolumeInCube(InputImageType::IndexType index, double voxelSize) const
{
  
}


} // end namespace itk
  
#endif
