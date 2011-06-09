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

#include "itkImage.h"
#include "itkIndex.h"

namespace itk
{
template< typename TOutput,
          unsigned int VDimension >
class WhitakerSparseLevelSetBase :
    public LevelSetBase< Index< VDimension >,
                         VDimension,
                         TOutput >
{
public:
  typedef Index< VDimension >                       IndexType;
  typedef TOutput                                   OutputType;

  typedef WhitakerSparseLevelSetBase                Self;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;
  typedef LevelSetBase< IndexType, VDimension, OutputType >
                                                    Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(WhitakerSparseLevelSetBase, LevelSetBase);

  typedef typename Superclass::InputType      InputType;
  typedef typename Superclass::OutputRealType OutputRealType;
  typedef typename Superclass::GradientType   GradientType;
  typedef typename Superclass::HessianType    HessianType;

  typedef char NodeStatusType;

  struct NodeAttributeType
    {
    /** status of a given node (its value also define in which layer it is)*/
    NodeStatusType  m_Status;

    /** level set value for a given node */
    OutputType      m_Value;
    };

  typedef std::pair< IndexType, NodeAttributeType >   NodePairType;
  typedef std::list< NodePairType >                   NodeListType;
  typedef typename NodeListType::iterator             NodeListIterator;
  typedef typename NodeListType::const_iterator       NodeListConstIterator;

  typedef std::map< NodeStatusType, NodeListType >    SparseLayerMapType;
  typedef typename SparseLayerMapType::iterator       SparseLayerMapIterator;
  typedef typename SparseLayerMapType::const_iterator SparseLayerMapConstIterator;

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
      itkGenericExceptionMacro( << "this layer " << iId << " does not exist" );
      return NULL;
      }
    }

  itkSetObjectMacro( Image, SparseImageType );
  itkGetObjectMacro( Image, SparseImageType );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */

  itkConceptMacro( DoubleConvertible,
                    ( Concept::Convertible< OutputRealType, OutputType > ) );

  /** End concept checking */
#endif // ITK_USE_CONCEPT_CHECKING

protected:

  WhitakerSparseLevelSetBase() : Superclass()
  {
    InitializeLayers();
  }
  ~WhitakerSparseLevelSetBase() {}

  SparseImagePointer m_Image;
  SparseLayerMapType m_LayerList;

  void InitializeLayers()
    {
    this->m_LayerList[ -2 ] = NodeListType();
    this->m_LayerList[ -1 ] = NodeListType();
    this->m_LayerList[  0 ] = NodeListType();
    this->m_LayerList[  1 ] = NodeListType();
    this->m_LayerList[  2 ] = NodeListType();
    }

private:
  WhitakerSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkWhitakerSparseLevelSetBase_h
