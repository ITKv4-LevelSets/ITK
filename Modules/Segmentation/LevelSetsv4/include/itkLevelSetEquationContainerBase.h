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

#ifndef __itkLevelSetEquationContainerBase_h
#define __itkLevelSetEquationContainerBase_h

#include "itkObject.h"

namespace itk
{
template< class TTermContainer >
class LevelSetEquationContainerBase : public Object
{
public:
  typedef LevelSetEquationContainerBase Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationContainerBase, Object );

  typedef TTermContainer                            TermContainerType;
  typedef typename TermContainerType::Pointer       TermContainerPointer;

  typedef typename TermContainerType::InputImageType    InputImageType;
  typedef typename TermContainerType::InputImagePointer InputImagePointer;

  typedef typename TermContainerType::LevelSetOutputPixelType LevelSetOutputPixelType;

  void AddEquation( const unsigned int& iId, TermContainerPointer iEquation )
    {
    if ( iEquation.IsNotNull() )
      {
      m_Container[iId] = iEquation;
      if( iEquation->GetInput() )
        {
        m_Input = iEquation->GetInput();
        }
      this->Modified();
      }
    else
      {
      itkGenericExceptionMacro( <<"Term supplied is null" );
      }
    }

  TermContainerPointer GetEquation( const unsigned int& iId )
    {
    typename std::map< unsigned int, TermContainerPointer >::iterator
        it = m_Container.find( iId );
    if( it != m_Container.end() )
      {
      return it->second;
      }
    else
      {
      itkGenericExceptionMacro( <<"this equation does not exist" );
      return TermContainerPointer();
      }
    }

  void Update()
    {
    typedef typename std::map< unsigned int, TermContainerPointer >::iterator
        ContainerIterator;

    for( ContainerIterator it = m_Container.begin();
         it != m_Container.end();
         ++it )
      {
      (it->second )->Update();
      }
    }

  LevelSetOutputPixelType GetCFLContribution()
    {
    typedef typename std::map< unsigned int, TermContainerPointer >::iterator
        ContainerIterator;

    LevelSetOutputPixelType oValue = NumericTraits< LevelSetOutputPixelType >::max();

    for( ContainerIterator it = m_Container.begin();
         it != m_Container.end();
         ++it )
      {
      oValue = vnl_math_min( oValue, ( it->second )->GetCFLContribution() );
      }

    return oValue;
    }

  itkSetObjectMacro( Input, InputImageType );
  itkGetObjectMacro( Input, InputImageType );

protected:

  LevelSetEquationContainerBase() : m_Input( NULL ) {}
  ~LevelSetEquationContainerBase() {}

  std::map< unsigned int, TermContainerPointer >  m_Container;
  InputImagePointer                               m_Input;


private:
  LevelSetEquationContainerBase( const Self& );
  void operator = ( const Self& );

};
}

#endif // __itkLevelSetEquationContainerBase_h
