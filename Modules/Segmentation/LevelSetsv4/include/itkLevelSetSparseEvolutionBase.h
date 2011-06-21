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


#ifndef __itkLevelSetSparseEvolutionBase_h
#define __itkLevelSetSparseEvolutionBase_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkUpdateWhitakerSparseLevelSet.h"
#include <list>
#include "itkObject.h"
#include "itkImageFileWriter.h"

namespace itk
{
template< class TEquationContainer >
class LevelSetSparseEvolutionBase : public Object
{
public:
  typedef LevelSetSparseEvolutionBase      Self;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;
  typedef Object                           Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetSparseEvolutionBase, Object );

  typedef TEquationContainer                      EquationContainerType;
  typedef typename EquationContainerType::Pointer EquationContainerPointer;
  typedef typename EquationContainerType::TermContainerType
                                                  TermContainerType;
  typedef typename TermContainerType::Pointer     TermContainerPointer;

  typedef typename TermContainerType::TermType TermType;
  typedef typename TermType::Pointer           TermPointer;

  typedef typename TermContainerType::InputImageType InputImageType;
  typedef typename InputImageType::PixelType         InputImagePixelType;
  typedef typename InputImageType::Pointer           InputImagePointer;
  typedef typename InputImageType::RegionType        InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                     InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int, InputImageType::ImageDimension );

  typedef typename TermContainerType::LevelSetContainerType LevelSetContainerType;
  typedef typename LevelSetContainerType::IdentifierType    IdentifierType;
  typedef typename LevelSetContainerType::Pointer           LevelSetContainerPointer;
  typedef typename LevelSetContainerType::LevelSetContainerConstIteratorType
                                                            LevelSetContainerConstIteratorType;
  typedef typename LevelSetContainerType::LevelSetContainerIteratorType
                                                            LevelSetContainerIteratorType;

  typedef typename LevelSetContainerType::LevelSetType LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::ImageType             LevelSetImageType;
  typedef typename LevelSetType::OutputRealType        LevelSetOutputRealType;
  typedef typename LevelSetType::OutputType            LevelSetOutputType;
  typedef typename LevelSetImageType::Pointer          LevelSetImagePointer;

  typedef typename LevelSetType::NodePairType     NodePairType;
  typedef typename LevelSetType::NodeListIterator NodeListIterator;

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

  typedef UpdateWhitakerSparseLevelSet< ImageDimension, LevelSetOutputType > UpdateLevelSetFilterType;

  typedef typename LevelSetType::StatusImageType  StatusImageType;
  typedef typename LevelSetType::OutputImageType  OutputImageType;

  typedef typename UpdateLevelSetFilterType::Pointer                         UpdateLevelSetFilterPointer;
  typedef typename UpdateLevelSetFilterType::UpdateListType                  UpdateListType;

  // create another class which contains all the equations
  // i.e. it is a container of term container :-):
  // set the i^th term container
  // This container should also hold the LevelSetContainer
  //   void SetLevelSetEquations( EquationContainer );

  itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  void Update()
    {
    m_DomainMapFilter = m_LevelSetContainer->GetDomainMapFilter();

    // Get the image to be segmented
    m_InputImage = m_EquationContainer->GetInput();

    // Get the LevelSetContainer from the EquationContainer
    m_LevelSetContainer = m_EquationContainer->GetEquation( 0 )->GetTerm( 0 )->GetLevelSetContainer();

    //Run iteration
    this->GenerateData();
    }

  itkSetMacro( Alpha, LevelSetOutputRealType );
  itkGetMacro( Alpha, LevelSetOutputRealType );

  void SetTimeStep( const LevelSetOutputRealType& iDt )
    {
    if( iDt > NumericTraits< LevelSetOutputRealType >::epsilon() )
      {
      m_UserDefinedDt = true;
      m_Dt = iDt;
      this->Modified();
      }
    else
      {
      itkGenericExceptionMacro( <<"iDt should be > epsilon")
      }
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
  LevelSetSparseEvolutionBase() : m_NumberOfIterations( 0 ), m_NumberOfLevelSets( 0 ),
    m_InputImage( NULL ), m_EquationContainer( NULL ), m_LevelSetContainer( NULL ),
    m_UpdateBuffer( NULL ), m_DomainMapFilter( NULL ), m_Alpha( 0.9 ),
    m_Dt( 1. ), m_RMSChangeAccumulator( -1. ), m_UserDefinedDt( false )
  {
    m_UpdateBuffer =  new UpdateListType;
  }
  ~LevelSetSparseEvolutionBase()
  {
    delete m_UpdateBuffer;
  }

  unsigned int                m_NumberOfIterations;
  /// \todo is it useful?
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;

  // For sparse case, the update buffer needs to be the size of the active layer
  UpdateListType*             m_UpdateBuffer;
  DomainMapImageFilterPointer m_DomainMapFilter;

  LevelSetOutputRealType      m_Alpha;
  LevelSetOutputRealType      m_Dt;
  LevelSetOutputRealType      m_RMSChangeAccumulator;
  bool                        m_UserDefinedDt;

  void AllocateUpdateBuffer()
    {
    m_UpdateBuffer->clear();
    }


  void GenerateData()
    {
    AllocateUpdateBuffer();

    InitializeIteration();

    for( unsigned int iter = 0; iter < m_NumberOfIterations; iter++ )
      {
      m_RMSChangeAccumulator = 0;

      // one iteration over all container
      // update each level set based on the different equations provided
      ComputeIteration();

      ComputeDtForNextIteration();

      UpdateLevelSets();
      UpdateEquations();
      }
    }


  void InitializeIteration()
  {
    std::cout << "Initialize iteration" << std::endl;
    DomainIteratorType map_it = m_DomainMapFilter->m_LevelSetMap.begin();
    DomainIteratorType map_end = m_DomainMapFilter->m_LevelSetMap.end();

    while( map_it != map_end )
    {
      // std::cout << map_it->second.m_Region << std::endl;
      InputImageIteratorType it( m_InputImage, map_it->second.m_Region );
      it.GoToBegin();

      while( !it.IsAtEnd() )
      {
        // std::cout << it.GetIndex() << std::endl;
        IdListType lout = map_it->second.m_List;

        if( lout.empty() )
        {
          itkGenericExceptionMacro( <<"No level set exists at voxel" );
        }

        for( IdListIterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
        {
          m_EquationContainer->GetEquation( *lIt - 1 )->Initialize( it.GetIndex() );
        }
        ++it;
      }
      ++map_it;
    }
    m_EquationContainer->Update();
  }

  void ComputeIteration()
  {
    std::cout << "Compute iteration" << std::endl;
    LevelSetContainerIteratorType it = m_LevelSetContainer->Begin();
    while( it != m_LevelSetContainer->End() )
    {
      LevelSetPointer levelSet = it->second;
      NodeListIterator list_it = levelSet->GetListNode( 0 )->begin();
      NodeListIterator list_end = levelSet->GetListNode( 0 )->end();
      NodePairType p;
      while( list_it != list_end )
      {
        p = (*list_it);

        // TODO: Terms should update their values here dynamically
        // no need to call Update() later on
//         std::cout << p.first << std::endl;
        InputPixelRealType temp_update = m_EquationContainer->GetEquation( it->first )->Evaluate( p.first );

        // TODO: Need to index the correct levelset
        m_UpdateBuffer->push_back( temp_update );
//         std::cout << temp_update << std::endl;
        ++list_it;
      }
    ++it;
    }
  }


  void ComputeDtForNextIteration()
    {
    std::cout << "ComputeDtForNextIteration" << std::endl;
    if( !m_UserDefinedDt )
      {
      if( ( m_Alpha > NumericTraits< LevelSetOutputRealType >::Zero ) &&
          ( m_Alpha < NumericTraits< LevelSetOutputRealType >::One ) )
        {
        LevelSetOutputRealType contribution = m_EquationContainer->GetCFLContribution();

        if( contribution > NumericTraits< LevelSetOutputRealType >::epsilon() )
          {
          m_Dt = m_Alpha / contribution;
          }
        else
          {
          itkGenericExceptionMacro( << "contribution is too low" );
          }
        }
      else
        {
        itkGenericExceptionMacro( <<"m_Alpha should be in [0,1]" );
        }
      }

      std::cout << "Dt = " << m_Dt << std::endl;
//       m_Dt = 0.08;
    }

  virtual void UpdateLevelSets()
    {
      std::cout << "Update levelsets" << std::endl;
      LevelSetPointer levelSet = m_LevelSetContainer->GetLevelSet( 0 );

      typedef ImageFileWriter< StatusImageType > WriterType;
      typedef typename WriterType::Pointer       WriterPointer;

//       WriterPointer writer1 = WriterType::New();
//       writer1->SetInput( levelSet->GetOutputImage() );
//       writer1->SetFileName("/home/krm15/1.mha");
//       writer1->Update();

      UpdateLevelSetFilterPointer update_levelset = UpdateLevelSetFilterType::New();
      update_levelset->SetSparseLevelSet( levelSet );
      update_levelset->SetUpdate( m_UpdateBuffer );
      update_levelset->SetDt( m_Dt );
      update_levelset->Update();

//             NodeListIterator list_it = levelSet->GetListNode( 0 )->begin();
//             NodeListIterator list_end = levelSet->GetListNode( 0 )->end();
//             NodePairType p;
//             while( list_it != list_end )
//             {
//               p = (*list_it);
//               std::cout << p.first << std::endl;
//               ++list_it;
//             }

      WriterPointer writer2 = WriterType::New();
      writer2->SetInput( levelSet->GetStatusImage() );
      writer2->SetFileName("/home/krm15/2.mha");
      writer2->Update();

      m_RMSChangeAccumulator = update_levelset->GetRMSChangeAccumulator();

      m_UpdateBuffer->clear();
    }

  void UpdateEquations()
    {
    std::cout << "Update equations" << std::endl;
    InitializeIteration();
//     m_EquationContainer->Update();
    }

private:
  LevelSetSparseEvolutionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkLevelSetSparseEvolutionBase_h
