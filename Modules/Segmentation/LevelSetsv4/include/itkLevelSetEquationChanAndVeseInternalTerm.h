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
  }

  virtual void InitializeParameters()
  {
    m_TotalValue = NumericTraits< InputPixelRealType >::Zero;
    m_TotalH = NumericTraits< LevelSetOutputRealType >::Zero;
    this->SetUp();
  }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual void Initialize( const LevelSetInputIndexType& iP )
  {
    if( this->m_Heaviside.IsNotNull() )
      {
      InputPixelType pixel = this->m_Input->GetPixel( iP );

      LevelSetOutputRealType prod;
      this->ComputeProduct( iP, prod );
      this->Accumulate( pixel, prod );
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }
  }

  virtual void ComputeProduct( const LevelSetInputIndexType& iP,
                              LevelSetOutputRealType& prod )
  {
    LevelSetOutputRealType value = this->m_CurrentLevelSetPointer->Evaluate( iP );
    prod = this->m_Heaviside->Evaluate( -value );
  }

  virtual void ComputeProductTerm( const LevelSetInputIndexType& iP,
                                  LevelSetOutputRealType& prod )
  {}


  /* Performs the narrow-band update of the Heaviside function for each voxel. The
  characteristic function of each region is recomputed. Using the
  new H values, the previous c_i are updated. Used by only the sparse image
  filter */
  virtual void UpdatePixel( const LevelSetInputIndexType& iP,
                           const LevelSetOutputRealType & oldValue,
                           const LevelSetOutputRealType & newValue )
  {
    // For each affected h val: h val = new hval (this will dirty some cvals)
    InputPixelType input = this->m_Input->GetPixel( iP );

    LevelSetOutputRealType oldH = this->m_Heaviside->Evaluate( -oldValue );
    LevelSetOutputRealType newH = this->m_Heaviside->Evaluate( -newValue );
    LevelSetOutputRealType change = newH - oldH;

    // update the foreground constant for current level-set function
    this->m_TotalH += change;
    this->m_TotalValue += input * change;
  }


protected:
  LevelSetEquationChanAndVeseInternalTerm() : Superclass(),
    m_Mean( NumericTraits< InputPixelRealType >::Zero ),
    m_TotalValue( NumericTraits< InputPixelRealType >::Zero ),
    m_TotalH( NumericTraits< LevelSetOutputRealType >::Zero )
  {}

  virtual ~LevelSetEquationChanAndVeseInternalTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "Internal Chan And Vese term";
    }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP )
    {
    if( this->m_Heaviside.IsNotNull() )
      {
      LevelSetOutputRealType value =
          static_cast< LevelSetOutputRealType >( this->m_CurrentLevelSetPointer->Evaluate( iP ) );

      LevelSetOutputRealType d_val = this->m_Heaviside->EvaluateDerivative( -value );

      InputPixelType pixel = this->m_Input->GetPixel( iP );
      LevelSetOutputRealType prod = 1;

      ComputeProductTerm( iP, prod );
      LevelSetOutputRealType oValue = d_val * prod *
        static_cast< LevelSetOutputRealType >( ( pixel - m_Mean ) * ( pixel - m_Mean ) );

//       std::cout << value << ' ' << int(pixel) << ' ' << d_val << ' ' << prod << ' ' << oValue << std::endl;
      return oValue;
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }
    return NumericTraits< LevelSetOutputPixelType >::Zero;
    }


  void Accumulate( const InputPixelType& iPix,
                   const LevelSetOutputRealType& iH )
    {
    m_TotalValue += static_cast< InputPixelRealType >( iPix ) *
        static_cast< LevelSetOutputRealType >( iH );
    m_TotalH += static_cast< LevelSetOutputRealType >( iH );
    }

  InputPixelRealType      m_Mean;
  InputPixelRealType      m_TotalValue;
  LevelSetOutputRealType  m_TotalH;

private:
  LevelSetEquationChanAndVeseInternalTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
