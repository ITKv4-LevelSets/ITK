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

#ifndef __itkLevelSetEquationChanAndVeseInternalTermBase_h
#define __itkLevelSetEquationChanAndVeseInternalTermBase_h

#include "itkLevelSetEquationTermBase.h"
#include "itkHeavisideStepFunctionBase.h"
#include "itkNumericTraits.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationChanAndVeseInternalTermBase :
    public LevelSetEquationTermBase< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationChanAndVeseInternalTermBase         Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationTermBase< TInput,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationChanAndVeseInternalTermBase,
                LevelSetEquationTermBase );

  typedef TInput                                  InputType;
  typedef typename Superclass::InputPointer       InputPointer;
  typedef typename Superclass::InputPixelType     InputPixelType;
  typedef typename Superclass::InputPixelRealType InputPixelRealType;

  typedef TLevelSetContainer                              LevelSetContainerType;
  typedef typename Superclass::LevelSetContainerPointer   LevelSetContainerPointer;
  typedef typename Superclass::LevelSetType               LevelSetType;
  typedef typename Superclass::LevelSetPointer            LevelSetPointer;
  typedef typename Superclass::OutputType                 LevelSetOutputType;
  typedef typename Superclass::InputType                  LevelSetInputType;
  typedef typename Superclass::GradientType               GradientType;
  typedef typename Superclass::HessianType                HessianType;

  typedef HeavisideStepFunctionBase< LevelSetOutputType, LevelSetOutputType >
                                          HeavisideType;
  typedef typename HeavisideType::Pointer HeavisidePointer;

  itkSetObjectMacro( Heaviside, HeavisideType );
  itkGetObjectMacro( Heaviside, HeavisideType );

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
    m_TotalValue = NumericTraits< InputPixelRealType >::Zero;
    m_TotalH = NumericTraits< LevelSetOutputType >::Zero;
  }

protected:
  LevelSetEquationChanAndVeseInternalTermBase() : Superclass(),
    m_Heaviside( NULL ),
    m_CurrentLevelSetPointer( NULL ),
    m_Mean( NumericTraits< InputPixelRealType >::Zero ),
    m_TotalH( NumericTraits< LevelSetOutputType >::Zero ),
    m_TotalValue( NumericTraits< InputPixelRealType >::Zero )
  {}

  virtual ~LevelSetEquationChanAndVeseInternalTermBase() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "Internal Chan And Vese term";
    }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputType Value( const LevelSetInputType& iP )
    {
    if( m_CurrentLevelSetPointer.IsNull() )
      {
      m_CurrentLevelSetPointer =
          this->m_LevelSetContainer->GetLevelSet( this->m_CurrentLevelSet );

      if( m_CurrentLevelSetPointer.IsNull() )
        {
        itkWarningMacro(
              <<"m_CurrentLevelSet does not exist in the level set container" );

        return NumericTraits< LevelSetOutputType >::Zero;
        }
      }

    if( m_Heaviside.IsNotNull() )
      {
      LevelSetOutputType value = m_CurrentLevelSetPointer->Evaluate( iP );
      LevelSetOutputType h_val = this->m_Heaviside->Evaluate( value );

      if( h_val > 0.5 * NumericTraits< LevelSetOutputType >::One )
        {
        InputPixelType pixel = this->m_Input->GetPixel( iP );

        this->Accumulate( pixel, h_val );

        return  h_val *
            static_cast< LevelSetOutputType >( ( pixel - m_Mean ) * ( pixel - m_Mean ) );
        }
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }

    return NumericTraits< LevelSetOutputType >::Zero;
    }

      // This should be in Iteration class
      // Pointer to DomainMap with the two maps giving region and Id
      // Obtain the region of RegionMap[this->m_CurrentLevelSet]
      // Obtain the list of neighbor Levelsets ListOfLevelSets[this->m_CurrentLevelSet]
      // Iterate through the region
      // Iterate through all terms in the TermContainer
      // Compute internal and external and overlap terms
//     }

  void Accumulate( const InputPixelType& iPix,
                   const LevelSetOutputType& iH )
    {
    m_TotalValue += static_cast< InputPixelRealType >( iPix ) * static_cast< InputPixelRealType >( iH );
    m_TotalH += static_cast< InputPixelRealType >( iH );
    }

  HeavisidePointer    m_Heaviside;
  LevelSetPointer     m_CurrentLevelSetPointer;
  InputPixelRealType  m_Mean;
  LevelSetOutputType  m_TotalH;
  InputPixelRealType  m_TotalValue;

//   InputPointer m_Input;
//   LevelSetContainerPointer m_LevelSetContainer;
//   LevelSetOutputType m_Cofficient;
//   std::string m_TermName;

private:
  LevelSetEquationChanAndVeseInternalTermBase( const Self& );
  void operator = ( const Self& );
};
}
#endif
