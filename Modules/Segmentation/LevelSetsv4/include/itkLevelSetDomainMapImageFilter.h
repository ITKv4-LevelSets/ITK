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

#ifndef __itkLevelSetDomainMapImageFilter_h
#define __itkLevelSetDomainMapImageFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkLabelImageToLabelMapFilter.h"

namespace itk
{
/**
  \class LevelSetDomainMapImageFilter
  \tparam TInputImage Image where the pixel type is a container of ids
  \tparam TOutputImage Image where the pixel type is an integer to split the region
*/
template < class TInputImage >
class ITK_EXPORT LevelSetDomainMapImageFilter :
    public LabelImageToLabelMapFilter< TInputImage >
{
  public:
    typedef LevelSetDomainMapImageFilter              Self;
    typedef LabelImageToLabelMapFilter< TInputImage > Superclass;
    typedef SmartPointer< Self >                      Pointer;
    typedef SmartPointer< const Self >                ConstPointer;

    itkStaticConstMacro ( ImageDimension, unsigned int,
                          TInputImage::ImageDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( LevelSetDomainMapImageFilter, LabelImageToLabelMapFilter );

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

    typedef typename Superclass::InputImageType         InputImageType;
    typedef typename Superclass::OutputImageType        OutputImageType;
    typedef typename Superclass::InputImagePointer      InputImagePointer;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    typedef typename Superclass::InputImageRegionType   InputImageRegionType;
    typedef typename Superclass::InputImagePixelType    InputImagePixelType;
    typedef typename Superclass::IndexType              IndexType;

    typedef typename Superclass::OutputImagePointer       OutputImagePointer;
    typedef typename Superclass::OutputImageConstPointer  OutputImageConstPointer;
    typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
    typedef typename Superclass::OutputImagePixelType     OutputImagePixelType;
    typedef typename Superclass::LabelObjectType          LabelObjectType;
    typedef typename Superclass::LengthType               LengthType;

    /** ImageDimension constants */
    itkStaticConstMacro( InputImageDimension, unsigned int,
                         TInputImage::ImageDimension);

 protected:
    LevelSetDomainMapImageFilter();
    ~LevelSetDomainMapImageFilter() {}


  private:
    LevelSetDomainMapImageFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#endif
