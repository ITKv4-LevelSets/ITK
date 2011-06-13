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

#include "itkImage.h"
#include "itkLevelSetImageBase.h"

namespace itk
{
template< unsigned int VDimension >
class ShiSparseLevelSetBase :
    public LevelSetImageBase< Image< char, VDimension > >
{
public:
  typedef char                                    OutputType;
  typedef Image< OutputType, VDimension >         ImageType;
  typedef typename ImageType::Pointer             ImagePointer;

  typedef ShiSparseLevelSetBase                   Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef LevelSetImageBase< ImageType >          Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShiSparseLevelSetBase, LevelSetImageBase);

  typedef typename Superclass::InputType      InputType;
  typedef typename Superclass::OutputRealType OutputRealType;
  typedef typename Superclass::GradientType   GradientType;
  typedef typename Superclass::HessianType    HessianType;

  typedef std::pair< InputType, OutputType >        NodePairType;
  typedef std::list< NodePairType >                 NodeListType;
  typedef typename NodeListType::iterator           NodeListIterator;
  typedef typename NodeListType::const_iterator     NodeListConstIterator;

  typedef std::map< OutputType, NodeListType >        SparseLayerMapType;
  typedef typename SparseLayerMapType::iterator       SparseLayerMapIterator;
  typedef typename SparseLayerMapType::const_iterator SparseLayerMapConstIterator;

  /*
  GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }
  */

  NodeListType* GetListNode( const OutputType& iId )
    {
    typename SparseLayerMapType::iterator it = m_LayerList.find( iId );
    if( it != m_LayerList.end() )
      {
      return & (it->second);
      }
    else
      {
      itkGenericExceptionMacro( << "this layer " << iId << " does not exist" );
      return NULL;
      }
    }

  virtual void Initialize()
    {
    Superclass::Initialize();

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

    this->m_LayerList = LevelSet->m_LayerList;
    }

protected:

  ShiSparseLevelSetBase() : Superclass()
    {
    InitializeLayers();
    }
  virtual ~ShiSparseLevelSetBase() {}

  SparseLayerMapType m_LayerList;

  void InitializeLayers()
    {
    this->m_LayerList.clear();
    this->m_LayerList[ -1 ] = NodeListType();
    this->m_LayerList[  1 ] = NodeListType();
    }

private:
  ShiSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkShiSparseLevelSetBase_h
