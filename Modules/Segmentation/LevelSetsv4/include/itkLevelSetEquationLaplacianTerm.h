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

#ifndef __itkLevelSetEquationLaplacianTerm_h
#define __itkLevelSetEquationLaplacianTerm_h

#include "itkLevelSetEquationTermBase.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkVector.h"
#include "vnl/vnl_matrix_fixed.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationLaplacianTerm :
    public LevelSetEquationTermBase< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationLaplacianTerm       Self;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;
  typedef LevelSetEquationTermBase< TInput, TLevelSetContainer >
                                                Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationLaplacianTerm,
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

  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  /** Neighborhood radius type */
  typedef ZeroFluxNeumannBoundaryCondition< InputImageType > DefaultBoundaryConditionType;
  typedef typename ConstNeighborhoodIterator< InputImageType >::RadiusType RadiusType;
  typedef ConstNeighborhoodIterator< InputImageType, DefaultBoundaryConditionType > NeighborhoodType;

  typedef Vector< LevelSetOutputRealType, itkGetStaticConstMacro(ImageDimension) > NeighborhoodScalesType;

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
  LevelSetEquationLaplacianTerm() : Superclass(),
    m_CurrentLevelSetPointer( NULL )
  {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      m_NeighborhoodScales[i] = 1.0;
      }
  }

  virtual ~LevelSetEquationLaplacianTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "Laplacian term";
    }

  LevelSetOutputRealType LaplacianSpeed( const LevelSetInputIndexType& iP ) const
  {
    return ( static_cast< LevelSetOutputRealType >( this->m_Input->GetPixel(iP) ) );
  }

  virtual LevelSetOutputRealType Value( const LevelSetInputIndexType& iP )
  {
    LevelSetInputIndexType pA, pB;
    LevelSetInputIndexType pAa, pBa, pCa, pDa;
    LevelSetOutputRealType valueAa, valueBa, valueCa, valueDa;
    LevelSetOutputRealType ZERO = NumericTraits< LevelSetOutputRealType >::Zero;

    vnl_matrix_fixed< LevelSetOutputRealType,
                    itkGetStaticConstMacro(ImageDimension),
                    itkGetStaticConstMacro(ImageDimension) > m_dxy;

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      pA = pB = iP;
      pA[i] += 1;
      pB[i] -= 1;

      for ( unsigned int j = i + 1; j < ImageDimension; j++ )
        {
        pAa = pB;
        pAa[j] -= 1;

        pBa = pB;
        pBa[j] += 1;

        pCa = pA;
        pCa[j] -= 1;

        pDa = pA;
        pDa[j] += 1;

        valueAa = static_cast< LevelSetOutputRealType >( m_CurrentLevelSetPointer->Evaluate( pAa ) );
        valueBa = static_cast< LevelSetOutputRealType >( m_CurrentLevelSetPointer->Evaluate( pBa ) );
        valueCa = static_cast< LevelSetOutputRealType >( m_CurrentLevelSetPointer->Evaluate( pCa ) );
        valueDa = static_cast< LevelSetOutputRealType >( m_CurrentLevelSetPointer->Evaluate( pDa ) );

        m_dxy[i][j] = m_dxy[j][i] = 0.25 * ( valueAa - valueBa - valueCa + valueDa )
                                         * m_NeighborhoodScales[i] * m_NeighborhoodScales[j];
        }
      }

    LevelSetOutputRealType laplacianSpeed = this->LaplacianSpeed( iP );
    LevelSetOutputRealType laplacian = ZERO;

    // Compute the laplacian using the existing second derivative values
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      laplacian += m_dxy[i][i];
      }

    laplacian = laplacian * laplacianSpeed;
//     std::cout << iP << ' ' << laplacian << std::endl;

    return laplacian;
  }

  LevelSetPointer         m_CurrentLevelSetPointer;
  LevelSetOutputRealType  m_NeighborhoodScales[ImageDimension];

private:
  LevelSetEquationLaplacianTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
