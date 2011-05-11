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

#ifndef __itkSparseLevelSetBase_h
#define __itkSparseLevelSetBase_h

#include "itkLevelSetBase.h"

#include "itkImage.h"
#include "itkIndex.h"

#include <list>
#include <map>

namespace itk
{
/**
\class SparseLevelSetBase
\brief Abstract class for sparse level set representation
*/
template< typename TOutput,
          unsigned int VDimension >
class SparseLevelSetBase :
    public LevelSetBase< Index< VDimension >,
                         VDimension,
                         TOutput >
{
public:
  typedef Index< VDimension >           IndexType;
  typedef TOutput                       OutputType;

  typedef SparseLevelSetBase            Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef LevelSetBase< IndexType, VDimension, OutputType >
                                        Superclass;

  typedef typename Superclass::InputType    InputType;
  typedef typename Superclass::GradientType GradientType;
  typedef typename Superclass::HessianType  HessianType;

  typedef std::pair< IndexType, NodeAttributeType > NodePairType;
  typedef std::list< NodePairType >                 NodeListType;
  typedef typename NodeListType::iterator           NodeListIterator;
  typedef typename NodeListType::const_iterator     NodeListConstIterator;

  typedef std::map< NodeStatusType, NodeListType >    SparseLayerMapType;
  typedef typename SparseLayerMapType::iterator       SparseLayerMapIterator;
  typedef typename SparseLayerMapType::const_iterator SparseLayerMapConstIterator;

  typedef char NodeStatusType;

  struct NodeAttributeType
    {
    /** status of a given node (its value also define in which layer it is)*/
    NodeStatusType  m_Status;

    /** level set value for a given node */
    OutputType      m_Value;
    };

  typedef Image< NodeAttributeType, VDimension >  SparseImageType;
  typedef typename SparseImageType::Pointer       SparseImagePointer;

  char GetStatus( const InputType& iP ) const
    {
    NodeAttributeType temp = m_Image->GetPixel( iP );
    return temp.m_Status;
    }

  virtual OutputType Evaluate( const InputType& iP ) const
    {
    NodeAttributeType temp = m_Image->GetPixel( iP );
    return temp.m_Value;
    }

  virtual GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  virtual HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }

  NodeListType* GetListNode( const int& iId )
    {
    typename SparseLayerMapType::iterator it = m_LayerList.find( iId );
    if( it != m_LayerList.end() )
      {
      return & (it->second);
      }
    else
      {
      itkGenericExceptionMacro( << "this layer does not exist" );
      return NULL;
      }
    }

  itkSetObjectMacro( Image, SparseImageType );
  itkGetObjectMacro( Image, SparseImageType );

protected:
  SparseLevelSetBase()
  {
    this->InitializeLayers();
  }
  virtual ~SparseLevelSetBase() {}

  SparseImagePointer m_Image;
  SparseLayerMapType m_LayerList;

  virtual void InitializeLayers() = 0;

private:
  SparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkSparseLevelSetBase_h
