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

#ifndef __itkLevelSetEquationChanAndVeseInternalTerm_h
#define __itkLevelSetEquationChanAndVeseInternalTerm_h

#include "itkLevelSetEquationTermBase.h"
#include "itkNumericTraits.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationChanAndVeseInternalTerm :
    public LevelSetEquationTermBase< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationChanAndVeseInternalTerm         Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationTermBase< TInput,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationChanAndVeseInternalTerm,
                LevelSetEquationTermBase );

  typedef TInput                                  InputType;
  typedef typename InputType::Pointer             InputPointer;
  typedef typename InputType::PixelType           InputPixelType;
  typedef typename NumericTraits< InputPixelType >::RealType
                                                  InputPixelRealType;

  typedef TLevelSetContainer                              LevelSetContainerType;
  typedef typename LevelSetContainerType::Pointer         LevelSetContainerPointer;
  typedef typename LevelSetContainerType::LevelSetType    LevelSetType;
  typedef typename LevelSetContainerType::LevelSetPointer LevelSetPointer;
  typedef typename LevelSetContainerType::OutputType      LevelSetOutputType;
  typedef typename LevelSetContainerType::InputType       LevelSetInputType;
  typedef typename LevelSetContainerType::GradientType    GradientType;
  typedef typename LevelSetContainerType::HessianType     HessianType;

  typedef typename Superclass::HeavisideType HeavisideType;
  typedef typename HeavisideType::Pointer    HeavisidePointer;

  itkSetMacro( Mean, InputPixelRealType );
  itkGetMacro( Mean, InputPixelRealType );

  virtual void Update()
  {
    if( m_TotalH > NumericTraits< LevelSetOutputType >::epsilon() )
      {
      m_Mean = m_TotalValue / m_TotalH;
      }
    else
      {
      m_Mean = NumericTraits< InputPixelRealType >::Zero;
      }

    std::cout << m_TotalValue << '/' << m_TotalH << '=' << m_Mean << std::endl;
    m_TotalValue = NumericTraits< InputPixelRealType >::Zero;
    m_TotalH = NumericTraits< LevelSetOutputType >::Zero;
    this->m_CFLContribution = NumericTraits< LevelSetOutputType >::Zero;
  }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual void Initialize( const LevelSetInputType& iP )
  {
    if( m_CurrentLevelSetPointer.IsNull() )
    {
      m_CurrentLevelSetPointer =
      this->m_LevelSetContainer->GetLevelSet( this->m_CurrentLevelSet );

      if( m_CurrentLevelSetPointer.IsNull() )
      {
        itkWarningMacro(
        << "m_CurrentLevelSet does not exist in the level set container" );
      }
    }

    if( this->m_Heaviside.IsNotNull() )
    {
      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputType prod;
      this->ComputeProduct( iP, prod );
      this->Accumulate( pixel, prod );
    }
    else
    {
      itkWarningMacro( << "m_Heaviside is NULL" );
    }
  }

  virtual void ComputeProduct( const LevelSetInputType& iP, LevelSetOutputType& prod )
  {
    LevelSetOutputType value = m_CurrentLevelSetPointer->Evaluate( iP );
    prod = this->m_Heaviside->Evaluate( -value );
  }

  virtual void ComputeProductTerm( const LevelSetInputType& iP, LevelSetOutputType& prod )
  {}

protected:
  LevelSetEquationChanAndVeseInternalTerm() : Superclass(),
    m_CurrentLevelSetPointer( NULL ),
    m_Mean( NumericTraits< InputPixelRealType >::Zero ),
    m_TotalH( NumericTraits< LevelSetOutputType >::Zero ),
    m_TotalValue( NumericTraits< InputPixelRealType >::Zero )
  {}

  virtual ~LevelSetEquationChanAndVeseInternalTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "Internal Chan And Vese term";
    }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputType Value( const LevelSetInputType& iP )
    {
    if( this->m_Heaviside.IsNotNull() )
      {
      LevelSetOutputType value = m_CurrentLevelSetPointer->Evaluate( iP );
      LevelSetOutputType d_val = this->m_Heaviside->EvaluateDerivative( -value );

      InputPixelType pixel = this->m_Input->GetPixel( iP );
      LevelSetOutputType prod = 1;

      ComputeProductTerm( iP, prod );
      LevelSetOutputType oValue = d_val * prod *
        static_cast< LevelSetOutputType >( ( pixel - m_Mean ) * ( pixel - m_Mean ) );

      return oValue;
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }
    return NumericTraits< LevelSetOutputType >::Zero;
    }


  void Accumulate( const InputPixelType& iPix,
                   const LevelSetOutputType& iH )
    {
    m_TotalValue += static_cast< InputPixelRealType >( iPix ) * static_cast< InputPixelRealType >( iH );
    m_TotalH += static_cast< InputPixelRealType >( iH );
    }

  LevelSetPointer     m_CurrentLevelSetPointer;
  InputPixelRealType  m_Mean;
  LevelSetOutputType  m_TotalH;
  InputPixelRealType  m_TotalValue;

//   InputPointer m_Input;
//   LevelSetContainerPointer m_LevelSetContainer;
//   LevelSetOutputType m_Cofficient;
//   std::string m_TermName;

private:
  LevelSetEquationChanAndVeseInternalTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
