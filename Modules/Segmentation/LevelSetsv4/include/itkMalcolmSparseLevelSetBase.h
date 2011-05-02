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

#ifndef __itkMalcolmSparseLevelSetBase_h
#define __itkMalcolmSparseLevelSetBase_h

#include "itkSparseLevelSetBase.h"

namespace itk
{
template< unsigned int VDimension >
class MalcolmSparseLevelSetBase :
    public SparseLevelSetBase< char, VDimension >
{
public:
  typedef MalcolmSparseLevelSetBase                 Self;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;
  typedef SparseLevelSetBase< TOutput, VDimension > Superclass;

protected:

  MalcolmSparseLevelSetBase() : Superclass() {}
  ~MalcolmSparseLevelSetBase() {}

  void InitializeLayers()
    {
    this->m_LayerList[ -1 ] = NodeListType();
    this->m_LayerList[  0 ] = NodeListType();
    this->m_LayerList[  1 ] = NodeListType();
    }

private:
  MalcolmSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkMalcolmSparseLevelSetBase_h
