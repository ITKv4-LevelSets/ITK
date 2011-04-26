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

  itkGetMacro( CFLContribution, LevelSetOutputType );

  itkSetStringMacro( TermName );
  itkGetStringMacro( TermName );

  virtual void Update() = 0;

protected:
  LevelSetEquationTermBase() : Superclass(),
    m_Coefficient( NumericTraits< LevelSetOutputType >::One ),
    m_CFLContribution( NumericTraits< LevelSetOutputType >::Zero )
  {}

  virtual ~LevelSetEquationTermBase() {}

  virtual void SetDefaultTermName() = 0;
  virtual LevelSetOutputType Value( const LevelSetInputType& iP ) = 0;

  InputPointer             m_Input;
  LevelSetContainerPointer m_LevelSetContainer;
  LevelSetIdentifierType   m_CurrentLevelSet;
  LevelSetOutputType       m_Coefficient;
  LevelSetOutputType       m_CFLContribution;
  std::string              m_TermName;

private:
  LevelSetEquationTermBase( const Self& );
  void operator = ( const Self& );
};
}
#endif
