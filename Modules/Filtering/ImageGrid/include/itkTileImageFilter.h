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
#ifndef __itkTileImageFilter_h
#define __itkTileImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkFixedArray.h"

namespace itk
{
/** \class TileImageFilter
 * \brief Tile multiple input images into a single output image.
 *
 * This filter will tile multiple images using a user-specified
 * layout. The tile sizes will be large enough to accommodate the
 * largest image for each tile. The layout is specified with the
 * SetLayout method. The layout has the same dimension as the output
 * image. If all entries of the layout are positive, the tiled output
 * will contain the exact number of tiles. If the layout contains a 0
 * in the last dimension, the filter will compute a size that will
 * accommodate all of the images. Empty tiles are filled with the
 * value specified with the SetDefault value method. The input images
 * must have a dimension less than or equal to the output image. The
 * output image have a larger dimension than the input images. This
 * filter can be used to create a volume from a series of inputs by
 * specifying a layout of 1,1,0.
 * \ingroup ITK-ImageGrid
 *
 * \wiki
 * \wikiexample{ImageProcessing/TileImageFilter,Tile multiple images into another image}
 * \endwiki
 */

template< class TInputImage, class TOutputImage >
class ITK_EXPORT TileImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard Self typedef */
  typedef TileImageFilter                                 Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(TileImageFilter, ImageToImageFilter);

  /** Image pixel value typedef. */
  typedef typename TInputImage::PixelType  InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;

  /** Image related typedefs. */
  typedef typename TInputImage::Pointer  InputImagePointer;
  typedef typename TOutputImage::Pointer OutputImagePointer;

  typedef typename TInputImage::SizeType    InputSizeType;
  typedef typename TInputImage::IndexType   InputIndexType;
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TOutputImage::SizeType   OutputSizeType;
  typedef typename TOutputImage::IndexType  OutputIndexType;
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image related typedefs. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  /** \class TileInfo
   * Define a tile structure
   * \ingroup ITK-ImageGrid
   */
  class TileInfo
  {
public:
    int                   m_ImageNumber;
    OutputImageRegionType m_Region;
    TileInfo():m_ImageNumber(-1) {}
  };

  typedef Image< TileInfo, itkGetStaticConstMacro(OutputImageDimension) > TileImageType;

  /** LayoutArray type. */
  typedef FixedArray< unsigned int, itkGetStaticConstMacro(OutputImageDimension) > LayoutArrayType;

  /** Set/Get the layout of the tiles. If the last Layout value is 0,
   * the filter will compute a value that will acoomodate all of the
   * images. */
  itkSetMacro(Layout, LayoutArrayType);
  itkGetConstMacro(Layout, LayoutArrayType);

  /** Set the pixel value for locations that are not covered by an
   * input image. The default default pixel value is Zero. */
  itkSetMacro(DefaultPixelValue, OutputPixelType);

  /** Get the pixel value for locations that are not covered by an
   * input image. */
  itkGetConstMacro(DefaultPixelValue, OutputPixelType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputEqualityComparableCheck,
                   ( Concept::EqualityComparable< OutputPixelType > ) );
  itkConceptMacro( SameTypeCheck,
                   ( Concept::SameType< InputPixelType, OutputPixelType > ) );
  itkConceptMacro( OutputOStreamWritableCheck,
                   ( Concept::OStreamWritable< OutputPixelType > ) );
  /** End concept checking */
#endif
protected:
  TileImageFilter();
  ~TileImageFilter(){}
  void PrintSelf(std::ostream & os, Indent indent) const;

  void GenerateInputRequestedRegion();

  void GenerateOutputInformation();

  void  GenerateData();

private:
  TileImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  typename TileImageType::Pointer m_TileImage;

  OutputPixelType m_DefaultPixelValue;

  LayoutArrayType m_Layout;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTileImageFilter.txx"
#endif

#endif
