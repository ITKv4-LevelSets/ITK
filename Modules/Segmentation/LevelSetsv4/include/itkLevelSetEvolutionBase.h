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

#ifndef __itkLevelSetEvolutionBase_h
#define __itkLevelSetEvolutionBase_h

#include "itkObject.h"

namespace itk
{
template< class TTermContainer >
class LevelSetEvolutionBase : public Object
{
public:
  typedef LevelSetEvolutionBase   Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  // set the input image (image to be segmented)
  void SetInput();

  // create another class which contains all the equations
  // i.e. it is a container of term container :-):
  // set the i^th term container
  void SetLevelSetEquations( EquationContainer );

  void Update()
    {
    this->GenerateData();
    }

protected:

  void GenerateData()
    {
    for( unsigned int iter = 0; iter < m_NumberOfIterations; iter++ )
      {
      // one iteration over all container
      // update each level set based on the different equations provided
      OneIteration();

      //ComputeCFL();

      //ComputeDtForNextIteration();

      // at first we would not reinitialize, then reinitialize at each iteration
      // then improve this later
      //ReInitialize();
      }
    }

  void OneIteration()
    {
    DomainMapImageFilterType::Pointer domainMapFilter = DomainMapImageFilterType::New();
    domainMapFilter->SetInput( id_image );
    domainMapFilter->Update();
    CacheImageType::Pointer output = domainMapFilter->GetOutput();
    std::cout << "Domain partition computed" << std::endl;

    typedef std::map< itk::IdentifierType, DomainMapImageFilterType::NounToBeDefined >::iterator DomainMapIterator;
    DomainMapIterator map_it = domainMapFilter->m_LevelSetMap.begin();
    DomainMapIterator map_end = domainMapFilter->m_LevelSetMap.end();

    LevelSetType::Pointer levelSet;
    ChanAndVeseTermType::Pointer eqTerm;
    while( map_it != map_end )
      {
      IdListImageType::RegionType temp_region = map_it->second.m_Region;

      itk::ImageRegionConstIteratorWithIndex<IdListImageType >
          temp_it( id_image, temp_region );
      temp_it.GoToBegin();

      while( !temp_it.IsAtEnd() )
        {
        std::cout << temp_it.GetIndex() << std::endl;
        IdListType lout = map_it->second.m_List;

        if( lout.empty() )
          {
          return EXIT_FAILURE;
          }

        for( IdListType::iterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
          {
            std::cout << *lIt << " ";
            levelSet = lscontainer->GetLevelSet( *lIt - 1);
            std::cout << levelSet->Evaluate( temp_it.GetIndex() ) << std::endl;

            if ( *lIt - 1 == 0 )
            {
              eqTerm =  dynamic_cast< ChanAndVeseTermType* >( termContainer->GetTerm( *lIt - 1 ).GetPointer() );
              std::cout << eqTerm->Evaluate( temp_it.GetIndex() ) << std::endl;
            }
          }
        std::cout << std::endl;
        ++temp_it;
        }
      ++map_it;
      }
    }

private:
};
}
#endif // __itkLevelSetEvolutionBase_h
