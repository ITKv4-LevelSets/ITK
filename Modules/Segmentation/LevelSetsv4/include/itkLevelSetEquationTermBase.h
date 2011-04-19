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


#ifndef __itkLevelSetEquationTermBase_h
#define __itkLevelSetEquationTermBase_h

#include "itkObject.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationTermBase : public Object
{
public:
  typedef LevelSetEquationTermBase   Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  typedef TInput                      InputType;
  typedef typename InputType::Pointer InputPointer;

  typedef TLevelSetContainer                             LevelSetContainerType;
  typedef typename LevelSetContainerType::IdentifierType LevelSetIdentifierType;
  typedef typename LevelSetContainerType::Pointer        LevelSetContainerPointer;
  typedef typename LevelSetContainerType::OutputType     LevelSetOutputType;
  typedef typename LevelSetContainerType::InputType      LevelSetInputType;
  typedef typename LevelSetContainerType::GradientType   GradientType;
  typedef typename LevelSetContainerType::HessianType    HessianType;

  itkSetMacro( Coefficient, LevelSetOutputType );
  itkGetMacro( Coefficient, LevelSetOutputType );

  itkSetMacro( CurrentLevelSet, LevelSetIdentifierType );
  itkGetMacro( CurrentLevelSet, LevelSetIdentifierType );

  itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  itkSetObjectMacro( Input, InputType );
  itkGetObjectMacro( Input, InputType );

  virtual LevelSetOutputType Evaluate( const LevelSetInputType& iP )
    {
    return m_Coefficient * this->Value( iP );
    }

  itkSetStringMacro( TermName );
  itkGetStringMacro( TermName );

  virtual void Update() = 0;

protected:
  LevelSetEquationTermBase() : Superclass(),
    m_Coefficient( NumericTraits< LevelSetOutputType >::One )
  {}

  virtual ~LevelSetEquationTermBase() {}

  virtual void SetDefaultTermName() = 0;
  virtual LevelSetOutputType Value( const LevelSetInputType& iP ) = 0;

  InputPointer             m_Input;
  LevelSetContainerPointer m_LevelSetContainer;
  LevelSetIdentifierType   m_CurrentLevelSet;
  LevelSetOutputType       m_Coefficient;
  std::string              m_TermName;

private:
  LevelSetEquationTermBase( const Self& );
  void operator = ( const Self& );
};
}
#endif
