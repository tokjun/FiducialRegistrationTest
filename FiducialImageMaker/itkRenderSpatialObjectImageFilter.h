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

  /** Run-time type information (and related methods).   */
  itkTypeMacro( RenderSpatialObjectImageFilter, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(Label, int);
  itkGetConstMacro(Label, int);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (Concept::Convertible<double, OutputPixelType>));
  /** End concept checking */
#endif


protected:
  RenderSpatialObjectImageFilter();
  ~RenderSpatialObjectImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

private:
  RenderSpatialObjectImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRenderSpatialObjectImageFilter.txx"
#endif
  
#endif
