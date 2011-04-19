/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-11

 Copyright (c) 2009-11, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

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
    return m_Container[iId];
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
      oValue += term_it->Evaluate( iP );
      ++term_it;
      }

    return iP;
    }

protected:
  LevelSetEquationTermContainerBase() : Superclass() {}

  virtual ~LevelSetEquationTermContainerBase() {}

  std::map< unsigned int, TermPointer > m_Container;
  InputPointer                          m_Input;

private:
  LevelSetEquationTermContainerBase( const Self& );
  void operator = ( const Self& );
};

}
#endif // __itkLevelSetEquationTermContainerBase_h
