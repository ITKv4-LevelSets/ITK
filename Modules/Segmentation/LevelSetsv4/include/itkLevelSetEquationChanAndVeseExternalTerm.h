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

#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkNumericTraits.h"

namespace itk
{
template< class TInput, // Input image
          class TLevelSetContainer >
class LevelSetEquationChanAndVeseExternalTerm :
    public LevelSetEquationChanAndVeseInternalTerm< TInput, TLevelSetContainer >
{
public:
  typedef LevelSetEquationChanAndVeseExternalTerm         Self;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  typedef LevelSetEquationChanAndVeseInternalTerm< TInput,
                                    TLevelSetContainer >  Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationChanAndVeseExternalTerm,
                LevelSetEquationChanAndVeseInternalTerm );

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

  virtual void ComputeProduct( const LevelSetInputIndexType& iP,
                               LevelSetOutputRealType& prod )
  {
    LevelSetIdentifierType id =
        this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );

    IdListType lout = this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

    LevelSetPointer levelSet;
    LevelSetOutputRealType value;
    prod = 1;
    for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
      {
      levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1 );
      value = levelSet->Evaluate( iP );
      prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
      }
  }

  virtual void ComputeProductTerm( const LevelSetInputIndexType& iP,
                                   LevelSetOutputRealType& prod )
  {
    prod = -1.;
    LevelSetIdentifierType id =
        this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
    IdListType lout =
        this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

    LevelSetPointer levelSet;
    LevelSetOutputRealType value;
    for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
      {
      if( *lIt-1 != this->m_CurrentLevelSet )
        {
        levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1);
        value = levelSet->Evaluate( iP );
        prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
        }
      }
  }

protected:
  LevelSetEquationChanAndVeseExternalTerm() : Superclass()
  {}

  virtual ~LevelSetEquationChanAndVeseExternalTerm() {}

  virtual void SetDefaultTermName()
    {
    this->m_TermName = "External Chan And Vese term";
    }

private:
  LevelSetEquationChanAndVeseExternalTerm( const Self& );
  void operator = ( const Self& );
};

}
#endif
