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

#ifndef __itkLevelSetImageBase_h
#define __itkLevelSetImageBase_h

#include "itkLevelSetBase.h"

namespace itk
{
template< class TImage >
class LevelSetImageBase :
  public LevelSetBase<
      typename TImage::IndexType,
      TImage::ImageDimension,
      typename TImage::PixelType >
  {
public:
  typedef TImage                        ImageType;
  typedef typename ImageType::Pointer   ImagePointer;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::PixelType PixelType;

  typedef LevelSetImageBase          Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef LevelSetBase< IndexType,
    ImageType::ImageDimension, PixelType >
                                     Superclass;

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( LevelSetImageBase, LevelSetBase );

  typedef typename Superclass::InputType      InputType;
  typedef typename Superclass::OutputType     OutputType;
  typedef typename Superclass::OutputRealType OutputRealType;
  typedef typename Superclass::GradientType   GradientType;
  typedef typename Superclass::HessianType    HessianType;

  itkSetObjectMacro( Image, ImageType );
  itkGetObjectMacro( Image, ImageType );

  virtual OutputType Evaluate( const InputType& iP ) const
    {
    return m_Image->GetPixel( iP );
    }

  virtual GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  virtual HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }

  struct LevelSetDataType
    {
    OutputType Value;
    GradientType Gradient;
    HessianType Hessian;
    };

  virtual LevelSetDataType* GetAllData( const InputType& iP ) const
    {
    LevelSetDataType* oData = new LevelSetDataType;
    oData->Value = m_Image->GetPixel( iP );

    return oData;
    }

  virtual void Initialize()
    {
    Superclass::Initialize();

    m_Image = 0;
    }

  virtual void CopyInformation(const DataObject *data)
    {
    Superclass::CopyInformation( data );

    const Self *LevelSet = NULL;

    try
      {
      LevelSet = dynamic_cast< const Self * >( data );
      }
    catch ( ... )
      {
      // LevelSet could not be cast back down
      itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::Mesh::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }
    }

  virtual void Graft( const DataObject* data )
    {
    Superclass::Graft( data );
    const Self *LevelSet = 0;

    try
      {
      LevelSet = dynamic_cast< const Self* >( data );
      }
    catch( ... )
      {
      // mesh could not be cast back down
      itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    this->m_Image = LevelSet->m_Image;
    }

protected:
  LevelSetImageBase() : Superclass() {}

  virtual ~LevelSetImageBase() {}

  ImagePointer m_Image;

private:
  LevelSetImageBase( const Self& );
  void operator = ( const Self& );
  };
}
#endif // __itkLevelSetImageBase_h
