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

#ifndef __itkLevelSetImageBase_hxx
#define __itkLevelSetImageBase_hxx

#include "itkLevelSetImageBase.h"

namespace itk
{
// ----------------------------------------------------------------------------
template< class TImage >
LevelSetImageBase< TImage >
::LevelSetImageBase() : Superclass(), m_Image( NULL )
{
  m_NeighborhoodScales.Fill( NumericTraits< OutputRealType >::One );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::SetImage( ImageType* iImage )
{
  m_Image = iImage;
  typename ImageType::SpacingType spacing = m_Image->GetSpacing();

  for( unsigned int dim = 0; dim < Dimension; dim++ )
    {
    m_NeighborhoodScales[dim] =
        NumericTraits< OutputRealType >::One / static_cast< OutputRealType >( spacing[dim ] );
    }
  this->Modified();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
typename LevelSetImageBase< TImage >::OutputType
LevelSetImageBase< TImage >::Evaluate( const InputType& iP ) const
{
  return m_Image->GetPixel( iP );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
typename LevelSetImageBase< TImage >::GradientType
LevelSetImageBase< TImage >::EvaluateGradient( const InputType& iP ) const
{
  OutputRealType center_value =
      static_cast< OutputRealType >( this->Evaluate( iP ) );

  InputType pA, pB;
  pA = pB = iP;

  GradientType dx_forward;
  GradientType dx_backward;

  for( unsigned int dim = 0; dim < Dimension; dim++ )
    {
    pA[dim] += 1;
    pB[dim] -= 1;

    OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
    OutputRealType valueB = static_cast< OutputRealType >( this->Evaluate( pB ) );
    OutputRealType scale = m_NeighborhoodScales[dim];

    dx_forward[dim] = ( valueA - center_value ) * scale;
    dx_backward[dim] = ( center_value - valueB ) * scale;

    pA[dim] = pB[dim] = iP[dim];
    }

  return dx_forward;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
typename LevelSetImageBase< TImage >::HessianType
LevelSetImageBase< TImage >
::EvaluateHessian( const InputType& iP ) const
{
  HessianType oHessian;

  OutputRealType center_value =
      static_cast< OutputRealType >( this->Evaluate( iP ) );

  InputType pA, pB, pAa, pBa, pCa, pDa;
  pA = pB = iP;

  for( unsigned int dim1 = 0; dim1 < Dimension; dim1++ )
    {
    pA[dim1] += 1;
    pB[dim1] -= 1;

    OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
    OutputRealType valueB = static_cast< OutputRealType >( this->Evaluate( pB ) );

    oHessian[dim1][dim1] = ( valueA + valueB - 2.0 * center_value )
        * vnl_math_sqr( m_NeighborhoodScales[dim1] );

    pAa = pB;
    pBa = pB;

    pCa = pA;
    pDa = pA;

    for( unsigned int dim2 = dim1 + 1; dim2 < Dimension; dim2++ )
      {
      pAa[dim2] -= 1;
      pBa[dim2] += 1;

      pCa[dim2] -= 1;
      pDa[dim2] += 1;

      OutputRealType valueAa = static_cast< OutputRealType >( this->Evaluate( pAa ) );
      OutputRealType valueBa = static_cast< OutputRealType >( this->Evaluate( pBa ) );
      OutputRealType valueCa = static_cast< OutputRealType >( this->Evaluate( pCa ) );
      OutputRealType valueDa = static_cast< OutputRealType >( this->Evaluate( pDa ) );

      oHessian[dim1][dim2] = oHessian[dim2][dim1] =
          0.25 * ( valueAa - valueBa - valueCa + valueDa )
          * m_NeighborhoodScales[dim1] * m_NeighborhoodScales[dim2];

      pAa[dim2] = pB[dim2];
      pBa[dim2] = pB[dim2];

      pCa[dim2] = pA[dim2];
      pDa[dim2] = pA[dim2];
      }
    pA[dim1] = iP[dim1];
    pB[dim1] = iP[dim1];
    }

  return oHessian;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::Evaluate( const InputType& iP, LevelSetDataType& ioData ) const
{
  // if it has not already been computed before
  if( !ioData.Value.m_Computed )
    {
    ioData.Value.m_Computed = true;
    ioData.Value.m_Value = m_Image->GetPixel( iP );
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::EvaluateGradient( const InputType& iP, LevelSetDataType& ioData ) const
{
  // if it has not already been computed before
  if( !ioData.Gradient.m_Computed )
    {
    ioData.Gradient.m_Computed = true;
    // compute the gradient

    if( !ioData.Value.m_Computed )
      {
      ioData.Value.m_Computed = true;
      ioData.Value.m_Value = this->Evaluate( iP );
      }

    OutputRealType center_value =
        static_cast< OutputRealType >( ioData.Value.m_Value );

    InputType pA, pB;
    pA = pB = iP;

    GradientType dx_forward;
    GradientType dx_backward;

    for( unsigned int dim = 0; dim < Dimension; dim++ )
      {
      pA[dim] += 1;
      pB[dim] -= 1;

      OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
      OutputRealType valueB = static_cast< OutputRealType >( this->Evaluate( pB ) );
      OutputRealType scale = m_NeighborhoodScales[dim];

      ioData.Gradient.m_Value[dim] = ( valueA - valueB ) * scale;

      pA[dim] = pB[dim] = iP[dim];
      }
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::EvaluateHessian( const InputType& iP, LevelSetDataType& ioData ) const
{
  if( !ioData.Hessian.m_Computed )
    {
    ioData.Hessian.m_Computed = true;

    if( !ioData.Value.m_Computed )
      {
      ioData.Value.m_Computed = true;
      ioData.Value.m_Value = this->Evaluate( iP );
      }

    // compute the hessian
    OutputRealType center_value =
        static_cast< OutputRealType >( ioData.Value.m_Value );

    InputType pA, pB, pAa, pBa, pCa, pDa;
    pA = pB = iP;

    bool backward = ioData.BackwardGradient.m_Computed;
    bool forward = ioData.ForwardGradient.m_Computed;

    for( unsigned int dim1 = 0; dim1 < Dimension; dim1++ )
      {
      pA[dim1] += 1;
      pB[dim1] -= 1;

      OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
      OutputRealType valueB = static_cast< OutputRealType >( this->Evaluate( pB ) );

      ioData.Hessian.m_Value[dim1][dim1] =
          ( valueA + valueB - 2.0 * center_value ) * vnl_math_sqr( m_NeighborhoodScales[dim1] );

      if( !backward )
        {
        ioData.BackwardGradient.m_Computed = true;
        ioData.BackwardGradient.m_Value[i] =
            ( center_value - valueB ) * m_NeighborhoodScales[i];
        }
      if( !forward )
        {
        ioData.ForwardGradient.m_Computed = true;
        ioData.ForwardGradient.m_Value[i]  =
            ( valueA - center_value ) * m_NeighborhoodScales[i];
        }

      pAa = pB;
      pBa = pB;

      pCa = pA;
      pDa = pA;

      for( unsigned int dim2 = dim1 + 1; dim2 < Dimension; dim2++ )
        {
        pAa[dim2] -= 1;
        pBa[dim2] += 1;

        pCa[dim2] -= 1;
        pDa[dim2] += 1;

        OutputRealType valueAa = static_cast< OutputRealType >( this->Evaluate( pAa ) );
        OutputRealType valueBa = static_cast< OutputRealType >( this->Evaluate( pBa ) );
        OutputRealType valueCa = static_cast< OutputRealType >( this->Evaluate( pCa ) );
        OutputRealType valueDa = static_cast< OutputRealType >( this->Evaluate( pDa ) );

        ioData.Hessian.m_Value[dim1][dim2] =
            ioData.Hessian.m_Value[dim2][dim1] =
            0.25 * ( valueAa - valueBa - valueCa + valueDa )
            * m_NeighborhoodScales[dim1] * m_NeighborhoodScales[dim2];

        pAa[dim2] = pB[dim2];
        pBa[dim2] = pB[dim2];

        pCa[dim2] = pA[dim2];
        pDa[dim2] = pA[dim2];
        }
      pA[dim1] = iP[dim1];
      pB[dim1] = iP[dim1];
      }
    }
}
// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::EvaluateForwardGradient( const InputType& iP, LevelSetDataType& ioData ) const
{
  if( !ioData.ForwardGradient.m_Computed )
    {
    ioData.ForwardGradient.m_Computed = true;

    // compute the gradient
    if( !ioData.Value.m_Computed )
      {
      ioData.Value.m_Computed = true;
      ioData.Value.m_Value = this->Evaluate( iP );
      }

    OutputRealType center_value =
        static_cast< OutputRealType >( ioData.Value.m_Value );

    InputType pA;
    pA = iP;

    GradientType dx;

    for( unsigned int dim = 0; dim < Dimension; dim++ )
      {
      pA[dim] += 1;

      OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
      OutputRealType scale = m_NeighborhoodScales[dim];

      dx[dim] = ( valueA - center_value ) * scale;

      pA[dim] = iP[dim];
      }
    ioData.ForwardGradient.m_Value = dx;
    }
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::EvaluateBackwardGradient( const InputType& iP, LevelSetDataType& ioData ) const
{
  if( !ioData.BackwardGradient.m_Computed )
    {
    ioData.BackwardGradient.m_Computed = true;

    // compute the gradient
    if( !ioData.Value.m_Computed )
      {
      ioData.Value.m_Computed = true;
      ioData.Value.m_Value = this->Evaluate( iP );
      }

    OutputRealType center_value =
        static_cast< OutputRealType >( ioData.Value.m_Value );

    InputType pA;
    pA = iP;

    GradientType dx;

    for( unsigned int dim = 0; dim < Dimension; dim++ )
      {
      pA[dim] -= 1;

      OutputRealType valueA = static_cast< OutputRealType >( this->Evaluate( pA ) );
      OutputRealType scale = m_NeighborhoodScales[dim];

      dx[dim] = ( center_value - valueA ) * scale;

      pA[dim] = iP[dim];
      }
    ioData.BackwardGradient.m_Value = dx;
    }
}

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::Initialize()
{
  Superclass::Initialize();

  m_Image = NULL;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::CopyInformation(const DataObject *data)
{
  Superclass::CopyInformation( data );

  const Self *LevelSet = NULL;

  try
    {
    LevelSet = dynamic_cast< const Self * >( data );
    }
  catch ( ... )
    {
    // LevelSet could not be cast back down
    itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }

  if ( !LevelSet )
    {
    // pointer could not be cast back down
    itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template< class TImage >
void
LevelSetImageBase< TImage >
::Graft( const DataObject* data )
{
  Superclass::Graft( data );
  const Self *LevelSet = NULL;

  try
    {
    LevelSet = dynamic_cast< const Self* >( data );
    }
  catch( ... )
    {
    // image could not be cast back down
    itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
    }

  if ( !LevelSet )
    {
    // pointer could not be cast back down
    itkExceptionMacro( << "itk::LevelSetImageBase::CopyInformation() cannot cast "
                       << typeid( data ).name() << " to "
                       << typeid( Self * ).name() );
    }

  this->m_Image = LevelSet->m_Image;
  this->m_NeighborhoodScales = LevelSet->m_NeighborhoodScales;
}
// ----------------------------------------------------------------------------

}
#endif // __itkLevelSetImageBase_hxx
