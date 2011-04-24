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
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
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
  typedef typename LevelSetType::ImageType             LevelSetImageType;
  typedef typename LevelSetImageType::Pointer          LevelSetImagePointer;
  typedef typename LevelSetImageType::PixelType        LevelSetImagePixelType;

  typedef BinaryThresholdImageFilter< LevelSetImageType, LevelSetImageType >
                                                       ThresholdFilterType;
  typedef typename ThresholdFilterType::Pointer        ThresholdFilterPointer;
  typedef SignedMaurerDistanceMapImageFilter< LevelSetImageType, LevelSetImageType >
                                                       MaurerType;
  typedef typename MaurerType::Pointer                 MaurerPointer;

  typedef ImageRegionIteratorWithIndex< LevelSetImageType > LevelSetImageIteratorType;

  typedef ImageRegionConstIteratorWithIndex< LevelSetImageType > LevelSetImageConstIteratorType;

  typedef ImageRegionIteratorWithIndex< InputImageType > InputImageIteratorType;

  typedef ImageRegionConstIteratorWithIndex< InputImageType > InputImageConstIteratorType;

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
    // Allocate update buffer
    this->AllocateUpdateBuffer();
    //Run iteration
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

  void AllocateUpdateBuffer()
    {
    this->m_UpdateBuffer = LevelSetContainerType::New();
    this->m_UpdateBuffer->CopyInformationAndAllocate( m_LevelSetContainer, true );
    }

  void ComputeIteration()
    {
    DomainIteratorType map_it = m_DomainMapFilter->m_LevelSetMap.begin();
    DomainIteratorType map_end = m_DomainMapFilter->m_LevelSetMap.end();

    std::cout << "Begin iteration" << std::endl;

    LevelSetPointer levelSet, levelSetUpdate;
    InputPixelRealType p;
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

            levelSetUpdate = m_UpdateBuffer->GetLevelSet( *lIt - 1);

            // Store this update.
            //Run through all the terms associated with a given equation.
            // Another loop needs to come here.
            // TODO: Dynamic cast problems here
            // std::cout << m_TermContainer->GetTerm( *lIt - 1 )->Evaluate( it.GetIndex() ) << std::endl;
            // p = m_TermContainer->GetTerm( *lIt - 1 )->Evaluate( it.GetIndex() );
            levelSetUpdate->GetImage()->SetPixel( p );
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
        m_RMSChangeAccumulator = 0;

        // one iteration over all container
        // update each level set based on the different equations provided
        ComputeIteration();

        //ComputeCFL();

        ComputeDtForNextIteration();
        ApplyUpdate();
        ReInitialize();
        }
      }

  void ComputeDtForNextIteration()
    {
    m_Dt = 0.2;
    }

  void ApplyUpdate()
    {
    LevelSetContainerIteratorType it1 = m_LevelSetContainer->m_Container.begin();
    LevelSetContainerConstIteratorType it2 = m_UpdateBuffer->m_Container.begin();

    LevelSetImagePixelType p;

    while( it1 != m_LevelSetContainer->m_Container.end() )
    {
      LevelSetImagePointer image1 = it1->second->GetInput();
      LevelSetImagePointer image2 = it2->second->GetInput();

      LevelSetImageIteratorType It1( image1, image1->GetBufferedRegion() );
      LevelSetImageIteratorType It2( image2, image2->GetBufferedRegion() );
      It1.GoToBegin();
      It2.GoToBegin();
      while( !It1.IsAtEnd() )
      {
        p = m_Dt * It2.Get();
        It1.Set( It1.Get() + p );
        m_RMSChangeAccumulator += p;
        ++It1;
        ++It2;
      }

      ++it1;
      ++it2;
    }
  }

  void Reinitialize()
  {
    LevelSetContainerIteratorType it = m_LevelSetContainer->m_Container.begin();

    while( it != m_LevelSetContainer->m_Container.end() )
    {
      LevelSetImagePointer image = it->second->GetInput();

      ThresholdFilterPointer thresh = ThresholdFilterType::New();
      thresh->SetLowerThreshold( NumericTraits< LevelSetImagePixelType >::NonpositiveMin() );
      thresh->SetUpperThreshold( 0 );
      thresh->SetInsideValue( 1 );
      thresh->SetOutsideValue( 0 );
      thresh->SetInput( image );
      thresh->Update();

      MaurerPointer maurer = MaurerType::New();
      maurer->SetInput( thresh->GetOutput() );
      maurer->SetSquaredDistance( false );
      maurer->SetUseImageSpacing( true );
      maurer->SetInsideIsPositive( true );
      maurer->Update();

      LevelSetImageIteratorType It1( image, image->GetBufferedRegion() );
      LevelSetImageIteratorType It2( maurer->GetOutput(), image->GetBufferedRegion() );
      It1.GoToBegin();
      It2.GoToBegin();
      while( !It1.IsAtEnd() )
      {
        It1.Set( It2.Get() );
        ++It1;
        ++It2;
      }
      ++it;
    }
  }

  unsigned int                m_NumberOfIterations;
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;
  LevelSetContainerPointer    m_UpdateBuffer;
  DomainMapImageFilterPointer m_DomainMapFilter;
  InputPixelRealType          m_Dt;
  InputPixelRealType          m_RMSChangeAccumulator;

private:
};
}
#endif // __itkLevelSetEvolutionBase_h
