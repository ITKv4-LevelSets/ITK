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

#ifndef __itkWhitakerSparseLevelSetBase_h
#define __itkWhitakerSparseLevelSetBase_h

#include "itkLevelSetBase.h"
#include "itkIndex.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkLabelObject.h"
#include "itkLabelMap.h"

namespace itk
{
/**
 *  \class WhitakerSparseLevelSetBase
 *  \brief Derived class for the sparse-field representation of level-set function
 *
 *  \tparam TOutput Output type (float or double) of the level set function
 *  \tparam VDimension Dimension of the input space
 */
template< typename TOutput,
          unsigned int VDimension >
class WhitakerSparseLevelSetBase :
    public LevelSetBase< Index< VDimension >,
                         VDimension,
                         TOutput,
                         ImageBase< VDimension > >
{
public:
  typedef Index< VDimension >                       IndexType;
  typedef TOutput                                   OutputType;
  typedef ImageBase< VDimension >                   ImageBaseType;

  typedef WhitakerSparseLevelSetBase                Self;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;
  typedef LevelSetBase< IndexType, VDimension, OutputType, ImageBaseType >
                                                    Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(WhitakerSparseLevelSetBase, LevelSetBase);

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

  /** Returns the layer affiliation of a given location iP */
  virtual char Status( const InputType& iP ) const
  {
    if( m_LabelMap.IsNotNull() )
      {
      return m_LabelMap->GetPixel( iP );
      }
    else
      {
      LayerMapConstIterator layerIt = m_Layers.begin();

      while( layerIt != m_Layers.end() )
        {
        LayerConstIterator it = ( layerIt->second ).find( iP );
        if( it != ( layerIt->second ).end() )
          {
          return layerIt->first;
          }

        ++layerIt;
        }

      return 3;
      }
  }

  /** Returns the value of the level set function at a given location iP */
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

    if( m_LabelMap.IsNotNull() )
      {
      if( m_LabelMap->GetLabelObject( -3 )->HasIndex( iP ) )
        {
        return -3;
        }
      else
        {
        char status = m_LabelMap->GetPixel( iP );
        if( status == 3 )
          {
          return 3.;
          }
        else
          {
          itkGenericExceptionMacro( <<"status "
                                    << static_cast< int >( status )
                                    << " should be 3 or -3" );
          return 3.;
          }
        }
      }
    else
      {
      std::cout << "*** Note: m_LabelMap is NULL ***" << std::endl;
      return 3;
      }
  }

  /** Returns the gradient of the level set function at a given location iP
   * \todo to be implemented */
  virtual GradientType EvaluateGradient( const InputType& iP ) const
    {
    itkWarningMacro( <<"to be implemented" );
    return GradientType();
    }

  /** Returns the Hessian of the level set function at a given location iP
   * \todo to be implemented */
  virtual HessianType EvaluateHessian( const InputType& iP ) const
    {
    itkWarningMacro( <<"to be implemented" );
    return HessianType();
    }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */

  itkConceptMacro( DoubleConvertible,
                    ( Concept::Convertible< OutputRealType, OutputType > ) );

  /** End concept checking */
#endif // ITK_USE_CONCEPT_CHECKING

  /** Initialize the label map point and the sparse-field layers */
  virtual void Initialize()
    {
    Superclass::Initialize();

    m_LabelMap = 0;
    m_Layers.clear();
    this->InitializeLayers();
    }

  /** Copy level set information from data object */
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
      itkExceptionMacro( << "itk::WhitakerSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::WhitakerSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }
    }

  /** Graft data object as level set object */
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
      itkExceptionMacro( << "itk::WhitakerSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::WhitakerSparseLevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    this->m_LabelMap->Graft( LevelSet->m_LabelMap );
    this->m_Layers = LevelSet->m_Layers;
    }

  /** Return the const pointer to a layer map with given id  */
  const LayerType& GetLayer( char iVal ) const
    {
    LayerMapConstIterator it = m_Layers.find( iVal );
    if( it == m_Layers.end() )
    {
      itkGenericExceptionMacro( <<"This layer does not exist" );
    }
    return it->second;
    }

  /** Return the pointer to a layer map with given id  */
  LayerType& GetLayer( char iVal )
    {
    LayerMapIterator it = m_Layers.find( iVal );
    if( it == m_Layers.end() )
    {
      itkGenericExceptionMacro( <<"This layer does not exist" );
    }
    return it->second;
    }

  /** Set a layer map with id to the given layer pointer */
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

  /** Return the label object pointer with a given id */
  template< class TLabel >
  LabelObject< TLabel, Dimension >* GetAsLabelObject()
    {
    typedef LabelObject< TLabel, Dimension > OutputLabelObjectType;
    typename OutputLabelObjectType::Pointer object =
        OutputLabelObjectType::New();

    for( char status = -3; status < 1; status++ )
      {
      LabelObjectPointer labelObject = m_LabelMap->GetLabelObject( status );

      for( SizeValueType i = 0; i < labelObject->GetNumberOfLines(); i++ )
        {
        object->AddLine( labelObject->GetLine( i ) );
        }
      }
    object->Optimize();

    return object.GetPointer();
    }

  /** Set/Get the label map for computing the sparse representation */
  itkSetObjectMacro( LabelMap, LabelMapType );
  itkGetObjectMacro( LabelMap, LabelMapType );

protected:

  WhitakerSparseLevelSetBase() : Superclass(), m_LabelMap( 0 )
  {
    InitializeLayers();
  }
  virtual ~WhitakerSparseLevelSetBase() {}

  LayerMapType     m_Layers;
  LabelMapPointer  m_LabelMap;

  /** Initialize the sparse field layers */
  void InitializeLayers()
    {
    this->m_Layers.clear();
    this->m_Layers[ -2 ] = LayerType();
    this->m_Layers[ -1 ] = LayerType();
    this->m_Layers[  0 ] = LayerType();
    this->m_Layers[  1 ] = LayerType();
    this->m_Layers[  2 ] = LayerType();
    }

private:
  WhitakerSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkWhitakerSparseLevelSetBase_h
