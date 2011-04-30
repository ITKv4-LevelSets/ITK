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

#include "itkSparseLevelSetBase.h"

namespace itk
{
template< class TImage >
class WhitakerSparseLevelSetBase :
    public SparseLevelSetBase< TImage >
{
public:
  typedef WhitakerSparseLevelSetBase    Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef SparseLevelSetBase< TImage >  Superclass;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */

  itkConceptMacro( DoubleConvertible,
                    ( Concept::Convertible< double, OutputType > ) );

  /** End concept checking */
#endif // ITK_USE_CONCEPT_CHECKING

protected:

  WhitakerSparseLevelSetBase() : Superclass() {}
  ~WhitakerSparseLevelSetBase() {}

  void InitializeLayers()
    {
    m_LayerList[ -2 ] = NodeListType();
    m_LayerList[ -1 ] = NodeListType();
    m_LayerList[  0 ] = NodeListType();
    m_LayerList[  1 ] = NodeListType();
    m_LayerList[  2 ] = NodeListType();
    }

private:
  WhitakerSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkWhitakerSparseLevelSetBase_h
