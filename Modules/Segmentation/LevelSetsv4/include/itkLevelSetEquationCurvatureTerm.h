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

#ifndef __itkLevelSetEquationCurvatureTerm_h
#define __itkLevelSetEquationCurvatureTerm_h

#include "itkLevelSetEquationTermBase.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationCurvatureTerm :
    public LevelSetEquationTermBase< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationCurvatureTerm         Self;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;
  typedef LevelSetEquationTermBase< TInput, TLevelSetContainer >
                                                Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationCurvatureTerm,
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
  {}

  virtual void InitializeParameters()
  {
    this->m_CFLContribution = NumericTraits< LevelSetOutputRealType >::Zero;
  }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual void Initialize( const LevelSetInputIndexType& iP )
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

    if( !this->m_Heaviside.IsNotNull() )
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }
  }

  virtual void UpdatePixel( LevelSetInputIndexType& iP, LevelSetOutputRealType & oldValue, LevelSetOutputRealType & newValue )
  {}

protected:
  LevelSetEquationCurvatureTerm() : Superclass(),
    m_CurrentLevelSetPointer( NULL )
  {}

  virtual ~LevelSetEquationCurvatureTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "Curvature term";
    }


  // this will work for scalars and vectors. For matrices, one may have to reimplement
  // his specialized term
  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP )
    {
    if( this->m_Heaviside.IsNotNull() )
      {
      LevelSetOutputRealType value =
          static_cast< LevelSetOutputRealType >( m_CurrentLevelSetPointer->Evaluate( iP ) );
      LevelSetOutputRealType d_val = this->m_Heaviside->EvaluateDerivative( -value );
      InputPixelType pixel = this->m_Input->GetPixel( iP );
      LevelSetOutputRealType oValue = NumericTraits< LevelSetOutputRealType >::Zero;

      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          if ( j != i )
            {
            oValue -= gd->m_dx[i] * gd->m_dx[j] * gd->m_dxy[i][j];
            oValue += gd->m_dxy[j][j] * gd->m_dx[i] * gd->m_dx[i];
            }
          }
        }

      if ( gd->m_GradMag > vnl_math::eps )
        {
        oValue /= gd->m_GradMag * gd->m_GradMag * gd->m_GradMag;
        }
      else
        {
        oValue /= 1 + gd->m_GradMagSqr;
        }


//       std::cout << value << ' ' << int(pixel) << ' ' << d_val << ' ' << prod << ' ' << oValue << std::endl;
      return oValue;
      }
    else
      {
      itkWarningMacro( << "m_Heaviside is NULL" );
      }
    return NumericTraits< LevelSetOutputPixelType >::Zero;
    }


  LevelSetPointer         m_CurrentLevelSetPointer;

private:
  LevelSetEquationCurvatureTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
