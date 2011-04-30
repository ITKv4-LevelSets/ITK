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
  typedef typename LevelSetContainerType::IdentifierType  IdentifierType;

  typedef std::list< IdentifierType >                    IdListType;
  typedef typename IdListType::iterator                  IdListIterator;

  typedef typename Superclass::HeavisideType HeavisideType;
  typedef typename HeavisideType::Pointer    HeavisidePointer;

  virtual void ComputeProduct( const LevelSetInputType& iP, LevelSetOutputType& prod )
  {
    IdentifierType id = this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
    IdListType lout = this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

    LevelSetPointer levelSet;
    LevelSetOutputType value;
    prod = 1;
    for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
    {
      if ( *lIt-1 != this->m_CurrentLevelSet )
      {
        levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1 );
        value = levelSet->Evaluate( iP );
        prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
      }
    }
  }

  virtual void ComputeProductTerm( const LevelSetInputType& iP, LevelSetOutputType& prod )
  {
    prod = -1.;
    IdentifierType id = this->m_LevelSetContainer->GetDomainMapFilter()->GetOutput()->GetPixel( iP );
    IdListType lout = this->m_LevelSetContainer->GetDomainMapFilter()->m_LevelSetMap[id].m_List;

    LevelSetPointer levelSet;
    LevelSetOutputType value;
    for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
    {
      levelSet = this->m_LevelSetContainer->GetLevelSet( *lIt - 1);
      value = levelSet->Evaluate( iP );
      prod *= (1 - this->m_Heaviside->Evaluate( -value ) );
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
