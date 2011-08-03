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
#include "itkHeavisideStepFunctionBase.h"
#include "itkNumericTraits.h"

namespace itk
{
template< class TInputImage, // Input image
          class TLevelSetContainer >
class LevelSetEquationTermBase : public Object
{
public:
  typedef LevelSetEquationTermBase   Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationTermBase, Object );

  typedef TInputImage                                     InputImageType;
  typedef typename InputImageType::Pointer                InputImagePointer;
  typedef typename InputImageType::PixelType              InputPixelType;
  typedef typename NumericTraits< InputPixelType >::RealType
                                                          InputPixelRealType;

  typedef TLevelSetContainer                              LevelSetContainerType;
  typedef typename LevelSetContainerType::IdentifierType  LevelSetIdentifierType;
  typedef typename LevelSetContainerType::Pointer         LevelSetContainerPointer;
  typedef typename LevelSetContainerType::LevelSetType    LevelSetType;
  typedef typename LevelSetContainerType::LevelSetPointer LevelSetPointer;
  //typedef typename LevelSetContainerType::OutputPixelType LevelSetOutputPixelType;
  typedef typename LevelSetContainerType::OutputType      LevelSetOutputPixelType;
  typedef typename LevelSetContainerType::OutputRealType  LevelSetOutputRealType;
  typedef typename LevelSetContainerType::InputIndexType  LevelSetInputIndexType;
  typedef typename LevelSetContainerType::GradientType    LevelSetGradientType;
  typedef typename LevelSetContainerType::HessianType     LevelSetHessianType;

  typedef HeavisideStepFunctionBase< LevelSetOutputRealType,
                                     LevelSetOutputRealType >
                                          HeavisideType;
  typedef typename HeavisideType::Pointer HeavisidePointer;

  /** Set/Get the image to be segmented */
  itkSetObjectMacro( Input, InputImageType );
  itkGetObjectMacro( Input, InputImageType );

  itkSetMacro( Coefficient, LevelSetOutputRealType );
  itkGetMacro( Coefficient, LevelSetOutputRealType );

  itkSetMacro( CurrentLevelSet, LevelSetIdentifierType );
  itkGetMacro( CurrentLevelSet, LevelSetIdentifierType );

//   itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  void SetLevelSetContainer( LevelSetContainerPointer ptr )
  {
    m_LevelSetContainer = ptr;
    m_Heaviside = ptr->GetHeaviside();
    this->Modified();
  }

  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  virtual LevelSetOutputRealType Evaluate( const LevelSetInputIndexType& iP )
    {
    return m_Coefficient * this->Value( iP );
    }

  virtual void Initialize( const LevelSetInputIndexType& iP ) = 0;

  virtual void InitializeParameters() = 0;

  virtual void UpdatePixel( const LevelSetInputIndexType& iP,
                           const LevelSetOutputRealType & oldValue,
                           const LevelSetOutputRealType & newValue ) = 0;
//    {
//    itkWarningMacro( << "This method has to be reimplemented in the inherited classes.");
//    }

  itkGetMacro( CFLContribution, LevelSetOutputRealType );

  itkSetStringMacro( TermName );
  itkGetStringMacro( TermName );

  virtual void Update() = 0;

protected:
  LevelSetEquationTermBase() : Superclass(),
    m_Coefficient( NumericTraits< LevelSetOutputRealType >::One ),
    m_CFLContribution( NumericTraits< LevelSetOutputRealType >::Zero )
  {}

  virtual ~LevelSetEquationTermBase() {}

  virtual void SetDefaultTermName() = 0;
  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP ) = 0;

  InputImagePointer        m_Input;
  LevelSetContainerPointer m_LevelSetContainer;
  LevelSetIdentifierType   m_CurrentLevelSet;
  LevelSetOutputRealType   m_Coefficient;
  LevelSetOutputRealType   m_CFLContribution;
  HeavisidePointer         m_Heaviside;
  std::string              m_TermName;

private:
  LevelSetEquationTermBase( const Self& );
  void operator = ( const Self& );
};
}
#endif