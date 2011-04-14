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

#ifndef ITKLEVELSETIMAGEBASE_H
#define ITKLEVELSETIMAGEBASE_H

#include "itkObject.h"

namespace itk
{
template<class TImage>
class LevelSetImageBase : public Object
  {
public:
  typedef LevelSetImageBase Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object Superclass;

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( LevelSetImageBase, Object );

  typedef TImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::PixelType OutputType;
  typedef typename ImageType::IndexType InputType;

  itkSetObjectMacro( Image, ImageType );
  itkGetObjectMacro( Image, ImageType );

  OutputType Evaluate( const InputType& iP ) const
    {
    return m_Image->GetPixel( iP );
    }

protected:
  LevelSetImageBase()
    {
    m_Image = ImageType::New();
    }

  ~LevelSetImageBase() {}

  ImagePointer m_Image;

private:
  LevelSetImageBase( const Self& );
  void operator = ( const Self& );
  };
}
#endif // ITKLEVELSETIMAGEBASE_H
