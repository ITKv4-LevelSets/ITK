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

#ifndef __itkLevelSetEquationAdvectionTermBase_h
#define __itkLevelSetEquationAdvectionTermBase_h

#include "itkImage.h"
#include "itkLevelSetEquationTermBase.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorCastImageFilter.h"

namespace itk
{
template< class TInputImage,
          class TAdvectionImage,
          class TLevelSetContainer >
class LevelSetEquationAdvectionTermBase :
    public LevelSetEquationTermBase< TInputImage, TLevelSetContainer >
{
public:
  typedef LevelSetEquationAdvectionTermBase               Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationTermBase< TInputImage,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationAdvectionTermBase,
                LevelSetEquationTermBase );

  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

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

  /** The vector type that will be used in the calculations. */
  typedef FixedArray< LevelSetOutputRealType, ImageDimension > VectorType;

  typedef Image< VectorType, ImageDimension >     AdvectionImageType;
  typedef typename AdvectionImageType::Pointer    AdvectionImagePointer;

  typedef VectorLinearInterpolateImageFunction< AdvectionImageType >
    AdvectionInterpolatorType;
  typedef typename AdvectionInterpolatorType::Pointer
    AdvectionInterpolatorPointer;

  virtual void Update() {}

  virtual void Initialize( const LevelSetInputIndexType& iP )
    {
    if( m_AdvectionImage.IsNotNull() )
      {
      if( !m_Initialized )
        {
        m_AdvectionImage->SetRequestedRegion(
              this->m_Input->GetRequestedRegion() );
        m_AdvectionImage->SetBufferedRegion(
              this->m_Input->GetBufferedRegion() );
        m_AdvectionImage->SetLargestPossibleRegion(
              this->m_Input->GetLargestPossibleRegion() );
        m_AdvectionImage->Allocate();

        m_VectorInterpolator->SetInputImage( m_AdvectionImage );

        m_Initialized = true;
        }
      }
    else
      {
      itkGenericExceptionMacro( <<"m_AdvectionImage is NULL" );
      }
    }

protected:
  LevelSetEquationAdvectionTermBase() : Superclass()
  {
    m_VectorInterpolator = AdvectionInterpolatorType::New();
    m_AdvectionImage = 0;
    m_Initialized = false;
  }
  virtual ~LevelSetEquationAdvectionTermBase() {}

  AdvectionInterpolatorPointer m_VectorInterpolator;

  /** The image holding the advection field for front propation */
  AdvectionImagePointer        m_AdvectionImage;

  /** A casting functor to convert between vector types.  */
  Functor::VectorCast< typename AdvectionInterpolatorType::OutputType,
                       VectorType > m_VectorCast;


  bool m_Initialized;

  virtual void SetDefaultTerm()
    {
    this->m_TermName = "Advection Term";
    }

  virtual VectorType AdvectionField( const LevelSetInputIndexType& iP )
    {
    if( m_VectorInterpolator->IsInsideBuffer( iP ) )
      {
      return m_VectorCast( m_VectorInterpolator->EvaluateAtContinuousIndex( iP ) );
      }
    return m_VectorCast( m_AdvectionImage->GetPixel( iP ) );
    }

  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP )
    {
    VectorType advection = this->AdvectionField( iP );
    LevelSetGradientType grad = this->m_CurrentLevelSet->EvaluateGradient1( iP );

    return advection * grad;
    }

private:
};
}
#endif // __itkLevelSetEquationAdvectionTermBase_h
