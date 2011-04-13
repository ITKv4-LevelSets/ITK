/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1550 $  // Revision of last commit
  Date: $Date: 2010-06-06 23:50:34 -0400 (Sun, 06 Jun 2010) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef __itkLevelSetDomainMapImageFilter_h
#define __itkLevelSetDomainMapImageFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <list>
#include <vector>

namespace itk
{
template < class TInputImage, class TOutputImage >
class ITK_EXPORT LevelSetDomainMapImageFilter : public ImageToImageFilter<
  TInputImage, TOutputImage >
{
  public:
    typedef LevelSetDomainMapImageFilter                      Self;
    typedef ImageToImageFilter< TInputImage,TOutputImage >    Superclass;
    typedef SmartPointer< Self >                              Pointer;
    typedef SmartPointer< const Self >                        ConstPointer;

    itkStaticConstMacro ( ImageDimension, unsigned int,
                          TInputImage::ImageDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( LevelSetDomainMapImageFilter, ImageToImageFilter );

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

    typedef TInputImage                           InputImageType;
    typedef typename InputImageType::Pointer      InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::PixelType    InputImagePixelType;
    typedef typename InputImageType::RegionType   InputImageRegionType;
    typedef typename InputImageType::SizeType     InputImageSizeType;
    typedef typename InputImageSizeType::SizeValueType
                                                  InputImageSizeValueType;
    typedef typename InputImageType::SpacingType  InputImageSpacingType;
    typedef typename InputImageType::IndexType    InputImageIndexType;
    typedef typename InputImageType::PointType    InputImagePointType;

    typedef TOutputImage                           OutputImageType;
    typedef typename OutputImageType::Pointer      OutputImagePointer;
    typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
    typedef typename OutputImageType::IndexType    OutputImageIndexType;
    typedef typename OutputImageType::PixelType    OutputImagePixelType;

    typedef ImageRegionConstIteratorWithIndex< InputImageType >
                                                           InputConstIteratorType;
    typedef ImageRegionIteratorWithIndex< InputImageType > InputIndexIteratorType;
    typedef ImageRegionIterator< InputImageType >          InputIteratorType;

    typedef ImageRegionConstIteratorWithIndex< OutputImageType >
                                                            OutputConstIteratorType;
    typedef ImageRegionIteratorWithIndex< OutputImageType > OutputIndexIteratorType;
    typedef ImageRegionIterator< OutputImageType >          OutputIteratorType;

//    struct Something // ~ kind of cache to speed up computations
//      {
//      IdentifierType m_Id;
//      InputImageRegionType m_Region;
//      InputImagePixelType m_List;
//      };

    // IdentifierType = 0 means it is background, there is no level sets
    // nothing to be done at this pixel...
    // Note that the identifier for a given index is given by the output image
    // of this filter
    std::map< IdentifierType, InputImageRegionType > m_SetOfRegions;
    std::map< IdentifierType, InputImagePixelType >  m_LevelSetList;

  protected:
    LevelSetDomainMapImageFilter();
    ~LevelSetDomainMapImageFilter() {}

    void GenerateData();

  private:
    LevelSetDomainMapImageFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkLevelSetDomainMapImageFilter.txx"
#endif
