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

#ifndef __itkLevelSetEquationTermContainerBase_h
#define __itkLevelSetEquationTermContainerBase_h

#include "itkLevelSetEquationTermBase.h"
#include "itkObject.h"

namespace itk
{
template< class TInput,
          class TLevelSetContainer >
class LevelSetEquationTermContainerBase : public Object
{
public:
  typedef LevelSetEquationTermContainerBase Self;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;
  typedef Object                            Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationTermContainerBase,
                Object );

  typedef TInput                      InputType;
  typedef typename InputType::Pointer InputPointer;

  typedef TLevelSetContainer                           LevelSetContainerType;
  typedef typename LevelSetContainerType::Pointer      LevelSetContainerPointer;
  typedef typename LevelSetContainerType::OutputType   LevelSetOutputType;
  typedef typename LevelSetContainerType::InputType    LevelSetInputType;
  typedef typename LevelSetContainerType::GradientType GradientType;
  typedef typename LevelSetContainerType::HessianType  HessianType;

  typedef LevelSetEquationTermBase< InputType, LevelSetContainerType > TermType;
  typedef typename TermType::Pointer                                   TermPointer;

  itkSetObjectMacro( Input, InputType );
  itkGetObjectMacro( Input, InputType );

  void AddTerm( const unsigned int& iId, TermPointer iTerm )
    {
    if ( iTerm.IsNotNull() )
      {
      if( iTerm->GetInput() == NULL )
        {
        if( m_Input.IsNotNull() )
          {
          iTerm->SetInput( m_Input );
          }
        else
          {
          itkGenericExceptionMacro( <<"m_Input and iTerm->GetInput are NULL" );
          }
        }
      m_Container[iId] = iTerm;
      this->Modified();
      }
    else
      {
      itkGenericExceptionMacro( <<"Term supplied is null" );
      }
    }

  TermPointer GetTerm( const unsigned int& iId )
    {
    typename std::map< unsigned int, TermPointer >::iterator
        it = m_Container.find( iId );

    if( it != m_Container.end() )
      {
      return it->second;
      }
    else
      {
      itkGenericExceptionMacro( <<"this term does not exist" );
      return TermPointer();
      }
    }

  LevelSetOutputType Evaluate( const LevelSetInputType& iP )
    {
    typename std::map< unsigned int, TermPointer >::iterator
        term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
        term_end = m_Container.end();

    LevelSetOutputType oValue = NumericTraits< LevelSetOutputType >::Zero;

    while( term_it != term_end )
      {
      oValue += ( term_it->second )->Evaluate( iP );
      ++term_it;
      }

    return oValue;
    }

  void Update()
    {
    typename std::map< unsigned int, TermPointer >::iterator
        term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
        term_end = m_Container.end();

    while( term_it != term_end )
      {
      ( term_it->second )->Update();
      ++term_it;
      }
    }

protected:
  LevelSetEquationTermContainerBase() : Superclass(), m_Input( NULL ) {}

  virtual ~LevelSetEquationTermContainerBase() {}

  std::map< unsigned int, TermPointer > m_Container;
  InputPointer                          m_Input;

private:
  LevelSetEquationTermContainerBase( const Self& );
  void operator = ( const Self& );
};

}
#endif // __itkLevelSetEquationTermContainerBase_h
