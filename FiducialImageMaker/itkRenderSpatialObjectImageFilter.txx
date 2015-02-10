/*=========================================================================

  Program:   Fiducial Image Maker CLI for 3D Slicer
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
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkListSample.h"
#include "itkCrossHelper.h"

#include "vnl/vnl_math.h"

#include <vector>

namespace itk
{

/**
 * Constructor
 */
template < typename  TInput, typename TOutput  >
RenderSpatialObjectImageFilter< TInput, TOutput >
::RenderSpatialObjectImageFilter()
{
  this->m_FiducialCenterList.clear();
  this->m_FiducialRadius = 10;
  this->m_DefaultVoxelValue = 100;
  this->m_ToleranceVolume = 0.01;
}

template < typename  TInput, typename TOutput  >
RenderSpatialObjectImageFilter< TInput, TOutput >
::~RenderSpatialObjectImageFilter()
{
}
  
template < typename  TInput, typename TOutput  >
void 
RenderSpatialObjectImageFilter< TInput, TOutput >
::GenerateData()
{
  itkDebugMacro(<< "RenderSpatialObjectImageFilter generating data ");
  
  typename InputImageType::ConstPointer input = this->GetInput();
  typename OutputImageType::Pointer output = this->GetOutput();

  this->AllocateOutputs();

  // TODO: Should outout image be initialized?
  output->FillBuffer(static_cast<OutputPixelType>(0));

  typename std::vector< typename InputImageType::PointType >::iterator fiter;


  // Define regions based on the locations of the markers
  for (fiter = this->m_FiducialCenterList.begin(); fiter != this->m_FiducialCenterList.end(); fiter ++)
    {
    ImageRegionConstIterator<InputImageType> it;
    it = ImageRegionConstIterator<InputImageType>(input, input->GetRequestedRegion());
    ImageRegionIterator<OutputImageType> oit;
    oit = ImageRegionIterator<OutputImageType>(output, output->GetRequestedRegion());

    oit.GoToBegin();
    it.GoToBegin();

  // TODO: Get voxel size
  typename InputImageType::SpacingType spacing;
  double voxelVolume;

  spacing = input->GetSpacing();
  voxelVolume = static_cast<double>(spacing[0]) * static_cast<double>(spacing[1]) * static_cast<double>(spacing[2]);

  // Indices for the 8 vertex of the given cubic region defined in the index space (i, j, k)
  // (not physical space (L, P, S))
  //     
  //            (v7)----+-----(v6)
  //            /|     /|     /|
  //           +------+------+ |
  //          /| |   /| |   /| |
  //        (v4)----+-----(v5)-+
  //         | |/|  | | |  | |/|
  //         | +----|-O----|-+ |
  //         |/| |  |/| |  |/| |
  //         +--(v3)+------+-|(v2)
  //         | |/   | |/   | |/
  //  2 1    | +----|-+----|-+
  //  |/     |/     |/     |/
  //  +--0  (v0)----+-----(v1)
  //
  //  The following vectors are defined:
  //     <iVec> = (<v1> - <v0>) / 2.0 (vector along i-direction (in the LPS coordinates)
  //     <jVec> = (<v3> - <v0>) / 2.0 (vector along j-direction (in the LPS coordinates)
  //     <kVec> = (<v4> - <v0>) / 2.0 (vector along k-direction (in the LPS coordinates)

  std::vector< VectorType > vertices;
  vertices.resize(8);

  // Get the vectors along the indices
  typename InputImageType::DirectionType direction = input->GetDirection();
  VectorType iVec;
  VectorType jVec;
  VectorType kVec;

  for (int i = 0; i < 3; i ++)
    {
    iVec[i] = direction[0][i] * spacing[0]/2.0;
    jVec[i] = direction[1][i] * spacing[1]/2.0;
    kVec[i] = direction[2][i] * spacing[2]/2.0;
    }

  VectorType ijVecPlus  = iVec + jVec;
  VectorType ijVecMinus = iVec - jVec;

  std::cerr << "Starting the loop" <<  std::endl;

  while (!it.IsAtEnd())
    {
    InputPixelType pix = it.Get();
    typename InputImageType::IndexType index = it.GetIndex();

    typename InputImageType::PointType point;
    input->TransformIndexToPhysicalPoint (index, point);
    VectorType pointVec = point.GetVectorFromOrigin();
    
    VectorType pointVecKMinus  = pointVec - kVec;
    VectorType pointVecKPlus   = pointVec + kVec;

    vertices[0] = pointVecKMinus - ijVecPlus;
    vertices[1] = pointVecKMinus + ijVecMinus;
    vertices[2] = pointVecKMinus + ijVecPlus;
    vertices[3] = pointVecKMinus - ijVecMinus;
    vertices[4] = pointVecKPlus - ijVecPlus;
    vertices[5] = pointVecKPlus + ijVecMinus;
    vertices[6] = pointVecKPlus + ijVecPlus;
    vertices[7] = pointVecKPlus - ijVecMinus;
    
    double partialVolume = this->ComputeObjectVolumeInCube(vertices, spacing);
    double percentage = partialVolume / voxelVolume;
   
    // TODO: How is the input image incorporated?
    //oit.Set( static_cast<OutputPixelType>(this->m_DefaultVoxelValue*percentage+it.Get()) );
    oit.Set( static_cast<OutputPixelType>(this->m_DefaultVoxelValue*percentage) );

    ++it;
    ++oit;
    
    if (percentage > 0.50)
      {
      std::cerr << "percentage = " << percentage << std::endl;
      }
    }
    
  vertices.clear();
  std::cerr << "End of loop." << std::endl;
  
}


template < typename  TInput, typename TOutput  >
void
RenderSpatialObjectImageFilter< TInput, TOutput >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


template < typename  TInput, typename TOutput  >
double
RenderSpatialObjectImageFilter< TInput, TOutput >
::ComputeObjectVolumeInCube(std::vector< VectorType >& vertices, typename InputImageType::SpacingType& spacing)
{

  double volume;
  if (vertices.size() != 8)
    {
    return 0.0;
    }

  // Assume the volume defined by vertices is a cube
  volume = static_cast<double>(spacing[0]) * static_cast<double>(spacing[1]) * static_cast<double>(spacing[2]);

  // Is the voxel is completely inside the object?
  int fInside = this->IsInsideObject(vertices);
  if (fInside == INSIDE_YES)
    {
    return volume;
    }

  else if (fInside == INSIDE_PARTIAL)
    {
    // Check if the divided voxel is still larger than the tolerance volume.
    if (volume < this->m_ToleranceVolume)
      {
      return volume;
      }
    else
      {
      // Divide the volume along the longest dimension and call this function recursively.

      // Find the largest dimension
      int index = 0;
      for (int i = 1; i < 3; i ++) index = (spacing[i] > spacing[index])? i : index;
      
      std::vector< VectorType > subvolumeVertices;
      subvolumeVertices.resize(8);

      int p0, p1, p2, p3, p4, p5, p6, p7;
      if (index == 0)
        {
        p0 = 0;
        p1 = 3;
        p2 = 7;
        p3 = 4;
        p4 = 1;
        p5 = 2;
        p6 = 6;
        p7 = 5;
        }
      else if (index == 1)
        {
        p0 = 0;
        p1 = 1;
        p2 = 5;
        p3 = 4;
        p4 = 3;
        p5 = 2;
        p6 = 6;
        p7 = 7;
        }
      else // if (index == 3)
        {
        p0 = 0;
        p1 = 1;
        p2 = 2;
        p3 = 3;
        p4 = 4;
        p5 = 5;
        p6 = 6;
        p7 = 7;
        }

      typename InputImageType::SpacingType subvolumeSpacing = spacing;
      subvolumeSpacing[index] /= 2.0;

      subvolumeVertices[p0] = vertices[p0];
      subvolumeVertices[p1] = vertices[p1];
      subvolumeVertices[p2] = vertices[p2];
      subvolumeVertices[p3] = vertices[p3];
      
      subvolumeVertices[p4] = (vertices[p0] + vertices[p4]) / 2.0;
      subvolumeVertices[p5] = (vertices[p1] + vertices[p5]) / 2.0;
      subvolumeVertices[p6] = (vertices[p2] + vertices[p6]) / 2.0;
      subvolumeVertices[p7] = (vertices[p3] + vertices[p7]) / 2.0;
      
      double volume1 = this->ComputeObjectVolumeInCube(subvolumeVertices, subvolumeSpacing);
      
      subvolumeVertices[p0] = subvolumeVertices[p4];
      subvolumeVertices[p1] = subvolumeVertices[p5];
      subvolumeVertices[p2] = subvolumeVertices[p6];
      subvolumeVertices[p3] = subvolumeVertices[p7];
      
      subvolumeVertices[p4] = vertices[p4];
      subvolumeVertices[p5] = vertices[p5];
      subvolumeVertices[p6] = vertices[p6];
      subvolumeVertices[p7] = vertices[p7];
      
      double volume2 = this->ComputeObjectVolumeInCube(subvolumeVertices, subvolumeSpacing);

      return (volume1 + volume2);

      }
    }
  else // fInside == INSIDE_NO
    {
    // Divide the volume in the longest dimension
    return 0.0;
    }

}


template < typename  TInput, typename TOutput  >
int
RenderSpatialObjectImageFilter< TInput, TOutput >
::IsInsideObject(std::vector< VectorType >& vertices)
{

  typename std::vector< typename InputImageType::PointType >::iterator fiter;

  // First, simply check if there is any vertex inside any spherical fiducial
  for (fiter = this->m_FiducialCenterList.begin(); fiter != this->m_FiducialCenterList.end(); fiter ++)
    {
    int count = 0;
    typename std::vector< VectorType >::iterator viter;
    for (viter = vertices.begin(); viter != vertices.end(); viter ++)
      {
      VectorType v = (*fiter).GetVectorFromOrigin() - *viter;
      if (v.GetNorm() <= this->m_FiducialRadius)
        {
        count ++;
        }
      }
    if (count == 6) // All vertices are in the sphere
      {
      return INSIDE_YES;
      }
    else if (count > 0) // At least one vertex is in the sphere
      {
      return INSIDE_PARTIAL;
      }
    }

  // If there is no vertex inside any spherical fiducial, check if there is any face of the cube that intersect the sphere.
  // Note that the vertices are ordered clock-wise about the normal vector.
  static const int faces[6][4] = {
    {0, 3, 2, 1},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {2, 3, 7, 6},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
  };

  for (fiter = this->m_FiducialCenterList.begin(); fiter != this->m_FiducialCenterList.end(); fiter ++)
    {
    for (int i = 0; i < 6; i ++)
      {
      VectorType& origin = vertices[faces[i][0]];
      typename InputImageType::PointType& fidCenter = *fiter;

      // Calculate the normal vector
      VectorType v1 = vertices[faces[i][1]] - origin;
      VectorType v2 = vertices[faces[i][3]] - origin;
      
      typedef itk::CrossHelper< VectorType > CrossType;
      CrossType cross;
      VectorType n = cross(v1, v2);
      n.Normalize();
      
      // Calculate the distance between the face plane (infinite plane) and the sphere center
      VectorType vc = fidCenter.GetVectorFromOrigin() - origin;
      VectorType::ValueType distance = vc * n;
      
      // If the distance is greater than the radius, the face does never intersect the sphere.
      // If the distance is less than the radius, the projection of the sphere center should be on the face.
      // (We have already examined the case, where one or more vertex is inside the sphere.)
      if (distance <= this->m_FiducialRadius)
        {
        // Calculate the projection of the sphere center onto the face, and the vector from the origin to the projected point
        VectorType p = fidCenter.GetVectorFromOrigin() - distance*n;
        VectorType vp = p - origin;
        
        // If <p> is on the face, the inner product of <p> and <v1> and that of <p> and <v2> should be positive,
        // but less than the length of <v1> and <v2> respectively. (Assuming that the face is a rectangle)
        VectorType::ValueType ipv1 = p*v1;
        VectorType::ValueType ipv2 = p*v2;
        if (ipv1 > 0 && ipv1 < v1.GetNorm() &&
            ipv2 > 0 && ipv2 < v2.GetNorm())
          {
          return INSIDE_PARTIAL;
          }
        }
      }
    }

  return INSIDE_NO;
}

  
} // end namespace itk
  
#endif
