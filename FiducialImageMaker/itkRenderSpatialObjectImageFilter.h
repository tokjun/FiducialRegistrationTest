/*=========================================================================

  Program:   Fiducial Image Maker CLI for 3D Slicer
  Module:    itkRenderSpatialObjectImageFilter.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  Copyright (c) Brigham and Women's Hosptial All. rights reserved.
  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRenderSpatialObjectImageFilter_h
#define __itkRenderSpatialObjectImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkCenteredAffineTransform.h"

namespace itk
{
/** \class RenderSpatialObjectImageFilter
 * \brief 
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
  
template < typename  TInput, typename TOutput  >
class ITK_EXPORT RenderSpatialObjectImageFilter : public
ImageToImageFilter< TInput, TOutput >
{
public:
  /** Standard class typedefs. */
  typedef RenderSpatialObjectImageFilter Self;
  typedef ImageToImageFilter<
          TInput,
          TOutput >                               Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  
  typedef typename Superclass::InputImageType            InputImageType;
  typedef typename Superclass::OutputImageType           OutputImageType;
  typedef typename InputImageType::PixelType             InputPixelType;
  typedef typename OutputImageType::PixelType            OutputPixelType;
  

  typedef typename itk::CenteredAffineTransform< float, 3 >      LineTransformType;

  typedef itk::Vector< double, 3 > VectorType;

  // Returned values for IsInsideObject(position, spacing)
  enum {
    INSIDE_NO,
    INSIDE_YES,
    INSIDE_PARTIAL,
  };

public:

  /** Run-time type information (and related methods).   */
  itkTypeMacro( RenderSpatialObjectImageFilter, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (Concept::Convertible<double, OutputPixelType>));
  /** End concept checking */
#endif

  void ClearFiducialCenterList() { this->m_FiducialCenterList.clear(); }
  
  void AddFiducialCenter(double p[3])
  {
    typename InputImageType::PointType point;
    for (int i = 0; i < 3; i ++) point[i] = p[i];
    this->m_FiducialCenterList.push_back(point);
  }
    
  itkSetMacro(FiducialRadius, double);
  itkGetConstMacro(FiducialRadius, double);

  itkSetMacro(DefaultVoxelValue, double);
  itkGetConstMacro(DefaultVoxelValue, double);

  itkSetMacro(ToleranceVolume, double);
  itkGetConstMacro(ToleranceVolume, double);

protected:
  RenderSpatialObjectImageFilter();
  ~RenderSpatialObjectImageFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  void SetupForVertexComputation(typename InputImageType::ConstPointer& input, typename InputImageType::SpacingType& spacing);
  void ComputeVerticesOfCubicRegion(typename InputImageType::PointType& center,
                                    std::vector< VectorType >& vertices);

  // Compute the volume of the object in the specified cube.
  // 'objectID' may be specified, if only one object is used.
  // (All object will be used when objectID < 0)
  double ComputeObjectVolumeInCube(std::vector< VectorType >& vertices,
                                   typename InputImageType::SpacingType& spacing,
                                   int objectID=-1);

  // Function to determine if the given region is completely inside (INSIDE_YES),
  // partially inside (INSIDE_PARTIAL), or completely outside (INSIDE_NO).
  // 'objectID' may be specified, if only one object is used.
  // (All object will be used when objectID < 0)
  virtual int IsInsideObject(std::vector< VectorType >& vertice,
                             int objectID=-1);

private:
  RenderSpatialObjectImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector< typename InputImageType::PointType > m_FiducialCenterList;
  double m_FiducialRadius;

  double m_DefaultVoxelValue;
  double m_ToleranceVolume;

  // Parameters to calculate vertex
  VectorType m_iVec;
  VectorType m_jVec;
  VectorType m_kVec;
  VectorType m_ijVecPlus;
  VectorType m_ijVecMinus;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRenderSpatialObjectImageFilter.txx"
#endif
  
#endif
