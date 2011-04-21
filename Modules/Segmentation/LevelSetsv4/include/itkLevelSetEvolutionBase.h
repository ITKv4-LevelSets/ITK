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

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include <list>
#include "itkObject.h"

namespace itk
{
template< class TTermContainer >
class LevelSetEvolutionBase : public Object
{
public:
  typedef LevelSetEvolutionBase      Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEvolutionBase, Object );

  typedef TTermContainer                      TermContainerType;
  typedef typename TermContainerType::Pointer TermContainerPointer;

  typedef typename TermContainerType::TermType TermType;
  typedef typename TermType::Pointer           TermPointer;

  typedef typename TermContainerType::InputType InputImageType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int, InputImageType::ImageDimension );

  typedef typename TermContainerType::LevelSetContainerType LevelSetContainerType;
  typedef typename LevelSetContainerType::IdentifierType    IdentifierType;
  typedef typename LevelSetContainerType::Pointer           LevelSetContainerPointer;

  typedef typename LevelSetContainerType::LevelSetType LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;

  typedef ImageRegionIteratorWithIndex< InputImageType > InputImageIteratorType;

  typedef std::list< IdentifierType >                    IdListType;
  typedef typename IdListType::iterator                  IdListIterator;
  typedef Image< IdListType, ImageDimension >            IdListImageType;
  typedef Image< short, ImageDimension >                 CacheImageType;
  typedef LevelSetDomainMapImageFilter< IdListImageType, CacheImageType >
                                                         DomainMapImageFilterType;
  typedef typename DomainMapImageFilterType::Pointer     DomainMapImageFilterPointer;
  typedef typename DomainMapImageFilterType::NounToBeDefined NounToBeDefined;
//   typedef typename DomainMapImageFilterType::DomainIteratorType DomainIteratorType;
typedef typename std::map< itk::IdentifierType, NounToBeDefined >::iterator DomainIteratorType;

  // create another class which contains all the equations
  // i.e. it is a container of term container :-):
  // set the i^th term container
  // This container should also hold the LevelSetContainer
//   void SetLevelSetEquations( EquationContainer );
  itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  void Update()
    {
    this->GenerateData();
    }

  // set the term container
  itkSetObjectMacro( TermContainer, TermContainerType );
  itkGetObjectMacro( TermContainer, TermContainerType );

  // set the number of iterations
  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetMacro( NumberOfIterations, unsigned int );

  // set the domain map image filter
  itkSetObjectMacro( DomainMapFilter, DomainMapImageFilterType );
  itkGetObjectMacro( DomainMapFilter, DomainMapImageFilterType );

protected:

  void ComputeIteration()
    {
    DomainIteratorType map_it = m_DomainMapFilter->m_LevelSetMap.begin();
    DomainIteratorType map_end = m_DomainMapFilter->m_LevelSetMap.end();
    std::cout << "Begin iteration" << std::endl;

    LevelSetPointer levelSet;
//     ChanAndVeseTermType::Pointer eqTerm;
    while( map_it != map_end )
      {
      std::cout << map_it->second.m_Region << std::endl;

      InputImageIteratorType it( m_InputImage, map_it->second.m_Region );
      it.GoToBegin();
      while( !it.IsAtEnd() )
        {
        std::cout << it.GetIndex() << std::endl;
        IdListType lout = map_it->second.m_List;

        if( lout.empty() )
          {
            itkGenericExceptionMacro( <<"No level set exists at voxel" );
          }

        for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
          {
            std::cout << *lIt << " ";
            levelSet = m_LevelSetContainer->GetLevelSet( *lIt - 1);
            std::cout << levelSet->Evaluate( it.GetIndex() ) << std::endl;

            // Store this update.
            //Run through all the terms associated with a given equation.
            // Another loop needs to come here.
            // TODO: Dynamic cast problems here
            // std::cout << m_TermContainer->GetTerm( *lIt - 1 )->Evaluate( it.GetIndex() ) << std::endl;
          }
        std::cout << std::endl;
        ++it;
        }
      ++map_it;
      }
    }

  void GenerateData()
    {
      m_InputImage = m_TermContainer->GetInput();

      // Get the LevelSetContainer from the EquationContainer
//       m_LevelSetContainer = m_EquationContainer->GetLevelSetContainer();
      for( unsigned int iter = 0; iter < m_NumberOfIterations; iter++ )
      {
        // one iteration over all container
        // update each level set based on the different equations provided
        ComputeIteration();

        //ComputeCFL();

        //ComputeDtForNextIteration();

        // Update

        // at first we would not reinitialize, then reinitialize at each iteration
        // then improve this later
        //ReInitialize();
      }
    }

  unsigned int                m_NumberOfIterations;
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  TermContainerPointer        m_TermContainer;
  LevelSetContainerPointer    m_LevelSetContainer;
  DomainMapImageFilterPointer m_DomainMapFilter;

//   EquationContainerPointer m_EquationContainer;

private:
};
}
#endif // __itkLevelSetEvolutionBase_h
