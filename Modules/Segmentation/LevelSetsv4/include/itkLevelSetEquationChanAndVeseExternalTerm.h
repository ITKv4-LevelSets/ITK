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

  typedef typename Superclass::InputImageType     InputImageType;
  typedef typename Superclass::InputImagePointer  InputImagePointer;
  typedef typename Superclass::InputPixelType     InputPixelType;
  typedef typename Superclass::InputPixelRealType InputPixelRealType;

  typedef typename Superclass::LevelSetContainerType      LevelSetContainerType;
  typedef typename Superclass::LevelSetContainerPointer   LevelSetContainerPointer;
  typedef typename Superclass::LevelSetType               LevelSetType;
  typedef typename Superclass::LevelSetPointer            LevelSetPointer;
  typedef typename Superclass::LevelSetOutputPixelType    LevelSetOutputPixelType;
  typedef typename Superclass::LevelSetOutputRealType     LevelSetOutputRealType;
  typedef typename Superclass::LevelSetInputIndexType     LevelSetInputIndexType;
  typedef typename Superclass::LevelSetGradientType       LevelSetGradientType;
  typedef typename Superclass::LevelSetHessianType        LevelSetHessianType;
  typedef typename Superclass::LevelSetIdentifierType     LevelSetIdentifierType;

  typedef std::list< LevelSetIdentifierType >            IdListType;
  typedef typename IdListType::iterator                  IdListIterator;

  typedef typename Superclass::HeavisideType    HeavisideType;
  typedef typename Superclass::HeavisidePointer HeavisidePointer;

  itkSetMacro( Mean, InputPixelRealType );
  itkGetMacro( Mean, InputPixelRealType );

  virtual void Update()
  {
    if( m_TotalH > NumericTraits< LevelSetOutputRealType >::epsilon() )
      {
      LevelSetOutputRealType inv_total_h = 1. / m_TotalH;

      // depending on the pixel type, it may be more efficient to do
      // a multiplication than to do a division
      m_Mean = m_TotalValue * inv_total_h;
      }
    else
      {
      m_Mean = NumericTraits< InputPixelRealType >::Zero;
      }

    std::cout << m_TotalValue << '/' << m_TotalH << '=' << m_Mean << std::endl;
    m_TotalValue = NumericTraits< InputPixelRealType >::Zero;
    m_TotalH = NumericTraits< LevelSetOutputRealType >::Zero;
    this->m_CFLContribution = NumericTraits< LevelSetOutputRealType >::Zero;
  }

  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual void Initialize( const LevelSetInputIndexType& iP )
  {
    LevelSetIdentifierType id =
        this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );

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
      LevelSetOutputPixelType value = m_CurrentLevelSetPointer->Evaluate( iP );

      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputPixelType prod = 1.;
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
    m_TotalH( NumericTraits< LevelSetOutputRealType >::Zero ),
    m_TotalValue( NumericTraits< InputPixelRealType >::Zero )
  {}

  virtual ~LevelSetEquationChanAndVeseExternalTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "External Chan And Vese term";
    }

  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP )
    {
    if( this->m_Heaviside.IsNotNull() )
      {
      LevelSetIdentifierType id =
          this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
      IdListType lout =
          this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

      LevelSetOutputRealType value =
          static_cast< LevelSetOutputRealType  >( m_CurrentLevelSetPointer->Evaluate( iP ) );
      LevelSetOutputRealType d_val = this->m_Heaviside->EvaluateDerivative( -value );

      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputRealType prod = NumericTraits< LevelSetOutputRealType >::One;

      for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
        {
        LevelSetPointer levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1);
        value = static_cast< LevelSetOutputRealType  >( levelSet->Evaluate( iP ) );
        prod *= ( 1 - this->m_Heaviside->Evaluate( -value ) );
        }
      LevelSetOutputRealType oValue = -d_val * prod *
        static_cast< LevelSetOutputRealType >( ( pixel - m_Mean ) * ( pixel - m_Mean ) );

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

  LevelSetPointer         m_CurrentLevelSetPointer;
  InputPixelRealType      m_Mean;
  LevelSetOutputRealType  m_TotalH;
  InputPixelRealType      m_TotalValue;

private:
  LevelSetEquationChanAndVeseExternalTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
