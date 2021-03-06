/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkGrayscaleDilateImageFilter_h
#define __itkGrayscaleDilateImageFilter_h

#include "itkKernelImageFilter.h"
#include "itkMovingHistogramDilateImageFilter.h"
#include "itkBasicDilateImageFilter.h"
#include "itkAnchorDilateImageFilter.h"
#include "itkVanHerkGilWermanDilateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkNeighborhood.h"

namespace itk
{
/**
 * \class GrayscaleDilateImageFilter
 * \brief gray scale dilation of an image
 *
 * Dilate an image using grayscale morphology. Dilation takes the
 * maximum of all the pixels identified by the structuring element.
 *
 * The structuring element is assumed to be composed of binary
 * values (zero or one). Only elements of the structuring element
 * having values > 0 are candidates for affecting the center pixel.
 *
 * \sa MorphologyImageFilter, GrayscaleFunctionDilateImageFilter, BinaryDilateImageFilter
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 * \ingroup ITK-MathematicalMorphology
 *
 * \wiki
 * \wikiexample{Morphology/GrayscaleDilateImageFilter,Dilate a grayscale image}
 * \endwiki
 */

template< class TInputImage, class TOutputImage, class TKernel >
class ITK_EXPORT GrayscaleDilateImageFilter:
  public KernelImageFilter< TInputImage, TOutputImage, TKernel >
{
public:
  /** Standard class typedefs. */
  typedef GrayscaleDilateImageFilter                              Self;
  typedef KernelImageFilter< TInputImage, TOutputImage, TKernel > Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(GrayscaleDilateImageFilter,
               KernelImageFilter);

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Image related typedefs. */
  typedef TInputImage                                InputImageType;
  typedef TOutputImage                               OutputImageType;
  typedef typename TInputImage::RegionType           RegionType;
  typedef typename TInputImage::SizeType             SizeType;
  typedef typename TInputImage::IndexType            IndexType;
  typedef typename TInputImage::PixelType            PixelType;
  typedef typename TInputImage::OffsetType           OffsetType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  typedef MovingHistogramDilateImageFilter< TInputImage, TOutputImage, TKernel >
  HistogramFilterType;
  typedef BasicDilateImageFilter< TInputImage, TOutputImage, TKernel >
  BasicFilterType;

  typedef FlatStructuringElement< itkGetStaticConstMacro(ImageDimension) > FlatKernelType;

  typedef AnchorDilateImageFilter< TInputImage, FlatKernelType >           AnchorFilterType;
  typedef VanHerkGilWermanDilateImageFilter< TInputImage, FlatKernelType > VHGWFilterType;
  typedef CastImageFilter< TInputImage, TOutputImage >                     CastFilterType;

  /** Typedef for boundary conditions. */
  typedef ImageBoundaryCondition< InputImageType > *      ImageBoundaryConditionPointerType;
  typedef ImageBoundaryCondition< InputImageType > const *ImageBoundaryConditionConstPointerType;
  typedef ConstantBoundaryCondition< InputImageType >     DefaultBoundaryConditionType;

  /** Kernel typedef. */
  typedef TKernel KernelType;
//   typedef typename KernelType::Superclass KernelSuperClass;
//   typedef Neighborhood< typename KernelType::PixelType, ImageDimension >
// KernelSuperClass;

  /** Set kernel (structuring element). */
  void SetKernel(const KernelType & kernel);

  /** Set/Get the boundary value. */
  void SetBoundary(const PixelType value);

  itkGetConstMacro(Boundary, PixelType);

  /** Set/Get the backend filter class. */
  void SetAlgorithm(int algo);

  itkGetConstMacro(Algorithm, int);

  /** GrayscaleDilateImageFilter need to set its internal filters as modified */
  virtual void Modified() const;

  /** define values used to determine which algorithm to use */
  enum {
    BASIC = 0,
    HISTO = 1,
    ANCHOR = 2,
    VHGW = 3
    } AlgorithmChoice;

  void SetNumberOfThreads(int nb);

protected:
  GrayscaleDilateImageFilter();
  ~GrayscaleDilateImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  void GenerateData();

private:
  GrayscaleDilateImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented

  PixelType m_Boundary;

  // the filters used internally
  typename HistogramFilterType::Pointer m_HistogramFilter;

  typename BasicFilterType::Pointer m_BasicFilter;

  typename AnchorFilterType::Pointer m_AnchorFilter;

  typename VHGWFilterType::Pointer m_VHGWFilter;

  // and the name of the filter
  int m_Algorithm;

  // the boundary condition need to be stored here
  DefaultBoundaryConditionType m_BoundaryCondition;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGrayscaleDilateImageFilter.txx"
#endif

#endif
