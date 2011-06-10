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

#ifndef __itkLevelSetEquationChanAndVeseExternalTerm_h
#define __itkLevelSetEquationChanAndVeseExternalTerm_h

#include "itkLevelSetEquationTermBase.h"
#include "itkNumericTraits.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationChanAndVeseExternalTerm :
    public LevelSetEquationTermBase< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationChanAndVeseExternalTerm         Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationTermBase< TInput,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationChanAndVeseExternalTerm,
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
  typedef typename LevelSetContainerType::OutputRealType  LevelSetOutputRealType;
  typedef typename LevelSetContainerType::InputType       LevelSetInputType;
  typedef typename LevelSetContainerType::GradientType    GradientType;
  typedef typename LevelSetContainerType::HessianType     HessianType;
  typedef typename LevelSetContainerType::IdentifierType  IdentifierType;

  typedef std::list< IdentifierType >                    IdListType;
  typedef typename IdListType::iterator                  IdListIterator;

  typedef typename Superclass::HeavisideType    HeavisideType;
  typedef typename Superclass::HeavisidePointer HeavisidePointer;

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
  }

  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual void Initialize( const LevelSetInputType& iP )
  {
    IdentifierType id = this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
    IdListType lout = this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

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
      LevelSetOutputType value = m_CurrentLevelSetPointer->Evaluate( iP );

      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputType prod = 1.;
      for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
      {
        if ( *lIt-1 != this->m_CurrentLevelSet )
        {
          LevelSetPointer levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1 );
          value = levelSet->Evaluate( iP );
          prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
        }
      }
      this->Accumulate( pixel, prod );
    }
    else
    {
      itkWarningMacro( << "m_Heaviside is NULL" );
    }
  }

//   virtual LevelSetOutputType Evaluate( const LevelSetInputType& iP )
//     {
//     return m_Coefficient * this->Value( iP );
//     }

protected:
  LevelSetEquationChanAndVeseExternalTerm() : Superclass(),
    m_CurrentLevelSetPointer( NULL ),
    m_Mean( NumericTraits< InputPixelRealType >::Zero ),
    m_TotalH( NumericTraits< LevelSetOutputType >::Zero ),
    m_TotalValue( NumericTraits< InputPixelRealType >::Zero )
  {}

  virtual ~LevelSetEquationChanAndVeseExternalTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "External Chan And Vese term";
    }

  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputRealType Value( const LevelSetInputType& iP )
    {
    if( this->m_Heaviside.IsNotNull() )
      {
      IdentifierType id = this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
      IdListType lout = this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

      LevelSetOutputType value = m_CurrentLevelSetPointer->Evaluate( iP );
      LevelSetOutputType d_val = this->m_Heaviside->EvaluateDerivative( -value );

      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputRealType prod = NumericTraits< LevelSetOutputRealType >::One;

      for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
        {
        LevelSetPointer levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1);
        value = levelSet->Evaluate( iP );
        prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
        }
      LevelSetOutputRealType oValue = -d_val * prod *
        static_cast< LevelSetOutputRealType >( ( pixel - m_Mean ) * ( pixel - m_Mean ) );

//       this->Accumulate( pixel, 1 - h_val );

      return oValue;
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }

    return NumericTraits< LevelSetOutputRealType >::Zero;
    }

  void Accumulate( const InputPixelType& iPix,
                   const LevelSetOutputRealType& iH )
    {
    m_TotalValue += static_cast< InputPixelRealType >( iPix ) *
        static_cast< InputPixelRealType >( iH );
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
  LevelSetEquationChanAndVeseExternalTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
