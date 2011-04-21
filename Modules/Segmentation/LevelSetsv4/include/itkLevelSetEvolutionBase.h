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
template< class TEquationContainer >
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

  typedef TEquationContainer                      EquationContainerType;
  typedef typename EquationContainerType::Pointer EquationContainerPointer;
  typedef typename EquationContainerType::TermContainerType
                                                  TermContainerType;
  typedef typename TermContainerType::Pointer     TermContainerPointer;

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
  itkSetObjectMacro( EquationContainer, EquationContainerType );
  itkGetObjectMacro( EquationContainer, EquationContainerType );

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
      m_InputImage = m_EquationContainer->GetInput();

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
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;
  DomainMapImageFilterPointer m_DomainMapFilter;

//   EquationContainerPointer m_EquationContainer;

private:
};
}
#endif // __itkLevelSetEvolutionBase_h
