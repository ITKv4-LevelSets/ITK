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

#ifndef __itkShiSparseLevelSetBase_h
#define __itkShiSparseLevelSetBase_h

#include "itkLevelSetBase.h"
#include "itkIndex.h"

#include "itkLabelObject.h"
#include "itkLabelMap.h"

namespace itk
{
template< unsigned int VDimension >
class ShiSparseLevelSetBase :
    public LevelSetBase< Index< VDimension >,
                         VDimension,
                         char,
                         ImageBase< VDimension > >
//    public LevelSetImageBase< Image< char, VDimension > >
{
public:
  typedef Index< VDimension >                     IndexType;
  typedef char                                    OutputType;
  typedef ImageBase< VDimension >                 ImageBaseType;

  typedef ShiSparseLevelSetBase                   Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef LevelSetBase< IndexType, VDimension,
                        OutputType, ImageBaseType >
                                                  Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShiSparseLevelSetBase, LevelSetBase);

  itkStaticConstMacro ( Dimension, unsigned int,
                        VDimension );

  typedef typename Superclass::InputType      InputType;
  typedef typename Superclass::OutputRealType OutputRealType;
  typedef typename Superclass::GradientType   GradientType;
  typedef typename Superclass::HessianType    HessianType;

  typedef LabelObject< char, VDimension >       LabelObjectType;
  typedef typename LabelObjectType::Pointer     LabelObjectPointer;
  typedef typename LabelObjectType::LengthType  LabelObjectLengthType;
  typedef typename LabelObjectType::LineType    LabelObjectLineType;

  typedef LabelMap< LabelObjectType >         LabelMapType;
  typedef typename LabelMapType::Pointer      LabelMapPointer;

  typedef std::map< IndexType, OutputType,
                    Functor::IndexLexicographicCompare< VDimension > >
                                                  LayerType;
  typedef typename LayerType::iterator            LayerIterator;
  typedef typename LayerType::const_iterator      LayerConstIterator;

  typedef std::map< char, LayerType >             LayerMapType;
  typedef typename LayerMapType::iterator         LayerMapIterator;
  typedef typename LayerMapType::const_iterator   LayerMapConstIterator;

  virtual char Status( const InputType& iP ) const
  {
    return m_LabelMap->GetPixel( iP );
  }

  virtual OutputType Evaluate( const InputType& iP ) const
  {
    LayerMapConstIterator layerIt = m_Layers.begin();

    while( layerIt != m_Layers.end() )
      {
      LayerConstIterator it = ( layerIt->second ).find( iP );
      if( it != ( layerIt->second ).end() )
        {
        return it->second;
        }

      ++layerIt;
      }

    if( m_LabelMap->GetLabelObject( -3 )->HasIndex( iP ) )
      {
      return -3;
      }
    else
      {
      char status = m_LabelMap->GetPixel( iP );
      if( status == 3 )
        {
        return 3;
        }
      else
        {
        itkGenericExceptionMacro( <<"status "
                                  << static_cast< int >( status )
                                  << " should be 3 or -3" );
        return 3;
        }
      }
  }

  virtual GradientType EvaluateGradient( const InputType& iP ) const
    {
    itkWarningMacro( <<"to be implemented" );
    return GradientType();
    }

  virtual HessianType EvaluateHessian( const InputType& iP ) const
    {
    itkWarningMacro( <<"to be implemented" );
    return HessianType();
    }

  const LayerType& GetLayer( char iVal ) const
    {
    LayerMapConstIterator it = m_Layers.find( iVal );
    if( it == m_Layers.end() )
    {
      itkGenericExceptionMacro( <<"This layer does not exist" );
    }
    return it->second;
    }

  LayerType& GetLayer( char iVal )
    {
    LayerMapIterator it = m_Layers.find( iVal );
    if( it == m_Layers.end() )
    {
      itkGenericExceptionMacro( <<"This layer does not exist" );
    }
    return it->second;
    }

  void SetLayer( char iVal, const LayerType& iLayer )
  {
    LayerMapIterator it = m_Layers.find( iVal );
    if( it != m_Layers.end() )
    {
      it->second = iLayer;
    }
    else
    {
      itkGenericExceptionMacro( <<iVal << "is out of bounds" );
    }
  }

  virtual void Initialize()
    {
    Superclass::Initialize();

    m_LabelMap = 0;
    m_Layers.clear();
    this->InitializeLayers();
    }

  virtual void CopyInformation( const DataObject* data )
    {
    Superclass::CopyInformation( data );

    const Self *LevelSet = NULL;

    try
      {
      LevelSet = dynamic_cast< const Self* >( data );
      }
    catch( ... )
      {
      // LevelSet could not be cast back down
      itkExceptionMacro( << "itk::ShiSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::ShiSparseLevelSetBase::CopyInformation() cannot cast "
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
      itkExceptionMacro( << "itk::ShiSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::ShiSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    this->m_LabelMap->Graft( LevelSet->m_LabelMap );
    this->m_Layers = LevelSet->m_Layers;
    }

  template< class TLabel >
  LabelObject< TLabel, Dimension >* GetAsLabelObject()
    {
    typedef LabelObject< TLabel, Dimension > OutputLabelObjectType;
    typename OutputLabelObjectType::Pointer object =
        OutputLabelObjectType::New();

    LabelObjectPointer labelObject = m_LabelMap->GetLabelObject( -3 );

    for( SizeValueType i = 0; i < labelObject->GetNumberOfLines(); i++ )
      {
      object->AddLine( labelObject->GetLine( i ) );
      }

    labelObject = m_LabelMap->GetLabelObject( -1 );

    for( SizeValueType i = 0; i < labelObject->GetNumberOfLines(); i++ )
      {
      object->AddLine( labelObject->GetLine( i ) );
      }
    object->Optimize();

    return object.GetPointer();
    }

  itkSetObjectMacro( LabelMap, LabelMapType );
  itkGetObjectMacro( LabelMap, LabelMapType );

protected:

  ShiSparseLevelSetBase() : Superclass(), m_LabelMap( 0 )
    {
    InitializeLayers();
    }
  virtual ~ShiSparseLevelSetBase() {}

  LayerMapType     m_Layers;
  LabelMapPointer  m_LabelMap;


  void InitializeLayers()
    {
    this->m_Layers.clear();
    this->m_Layers[ -1 ] = LayerType();
    this->m_Layers[  1 ] = LayerType();
    }

private:
  ShiSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkShiSparseLevelSetBase_h
