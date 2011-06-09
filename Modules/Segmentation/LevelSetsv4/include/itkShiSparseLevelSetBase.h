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
#include "itkIndex.h"
#include "itkLevelSetBase.h"

namespace itk
{
template< unsigned int VDimension >
class ShiSparseLevelSetBase :
    public LevelSetBase< Index< VDimension >,
                         VDimension,
                         char >
{
public:
  typedef Index< VDimension >                     InputType;
  typedef char                                    OutputType;

  typedef ShiSparseLevelSetBase                   Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef LevelSetBase< InputType,
                        VDimension,
                        OutputType >              Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShiSparseLevelSetBase, LevelSetBase);

  typedef typename Superclass::GradientType GradientType;
  typedef typename Superclass::HessianType  HessianType;

  typedef std::pair< InputType, OutputType >        NodePairType;
  typedef std::list< NodePairType >                 NodeListType;
  typedef typename NodeListType::iterator           NodeListIterator;
  typedef typename NodeListType::const_iterator     NodeListConstIterator;

  typedef std::map< OutputType, NodeListType >        SparseLayerMapType;
  typedef typename SparseLayerMapType::iterator       SparseLayerMapIterator;
  typedef typename SparseLayerMapType::const_iterator SparseLayerMapConstIterator;

  typedef Image< OutputType, VDimension >         SparseImageType;
  typedef typename SparseImageType::Pointer       SparseImagePointer;

  OutputType Evaluate( const InputType& iP ) const
    {
    return m_Image->GetPixel( iP );
    }

  GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }

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

  itkSetObjectMacro( Image, SparseImageType );
  itkGetObjectMacro( Image, SparseImageType );

protected:

  ShiSparseLevelSetBase() : Superclass()
    {
    InitializeLayers();
    }
  ~ShiSparseLevelSetBase() {}

  SparseImagePointer m_Image;
  SparseLayerMapType m_LayerList;

  void InitializeLayers()
    {
    this->m_LayerList[ -1 ] = NodeListType();
    this->m_LayerList[  1 ] = NodeListType();
    }

private:
  ShiSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkShiSparseLevelSetBase_h
