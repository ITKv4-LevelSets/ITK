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

#include "itkLevelSetImageBase.h"

#include <list>
#include <map>

namespace itk
{
template< class TImage >
class SparseLevelSetBase : public LevelSetImageBase< TImage >
{
public:
  typedef SparseLevelSetBase            Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef LevelSetImageBase< TImage >   Superclass;

  typedef TImage                        ImageType;
  typedef typename ImageType::Pointer   ImagePointer;

  typedef typename ImageType::IndexType NodeType;
  typedef std::list< NodeType >         NodeListType;

  typedef std::map< int, NodeListType > SparseLayerMapType;

protected:
  SparseLevelSetBase() {}
  virtual ~SparseLevelSetBase() {}

  virtual void InitializeLayers() = 0;

  SparseLayerMapType m_LayerList;

private:
  SparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkSparseLevelSetBase_h
