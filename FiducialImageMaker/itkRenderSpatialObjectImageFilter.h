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
  ~RenderSpatialObjectImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  double ComputeObjectVolumeInCube(std::vector< VectorType >& vertices);

  // Function to determine if the given region is completely inside (INSIDE_YES),
  // partially inside (INSIDE_PARTIAL), or completely outside (INSIDE_NO).
  virtual int IsInsideObject(std::vector< VectorType >& vertices);


private:
  RenderSpatialObjectImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector< typename InputImageType::PointType > m_FiducialCenterList;
  double m_FiducialRadius;

  double m_DefaultVoxelValue;
  double m_ToleranceVolume;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRenderSpatialObjectImageFilter.txx"
#endif
  
#endif
