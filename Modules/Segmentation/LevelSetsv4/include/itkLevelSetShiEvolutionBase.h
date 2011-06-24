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


#ifndef __itkLevelSetShiEvolutionBase_h
#define __itkLevelSetShiEvolutionBase_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkUpdateShiSparseLevelSet.h"
#include <list>
#include "itkObject.h"
#include "itkImageFileWriter.h"

namespace itk
{
template< class TEquationContainer >
class LevelSetShiEvolutionBase : public Object
{
public:
  typedef LevelSetShiEvolutionBase      Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetShiEvolutionBase, Object );

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

  typedef UpdateShiSparseLevelSet< ImageDimension, EquationContainerType > UpdateLevelSetFilterType;

  typedef typename LevelSetType::ImageType  OutputImageType;

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
    m_LevelSetContainer =
        m_EquationContainer->GetEquation( 0 )->GetTerm( 0 )->GetLevelSetContainer();

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
  LevelSetShiEvolutionBase() : m_NumberOfIterations( 0 ), m_NumberOfLevelSets( 0 ),
    m_InputImage( NULL ), m_EquationContainer( NULL ), m_LevelSetContainer( NULL ),
    m_DomainMapFilter( NULL ), m_Alpha( 0.9 ),
    m_Dt( 1. ), m_RMSChangeAccumulator( -1. ), m_UserDefinedDt( false )
  {
    m_UpdateBuffer[-1] = new UpdateListType;
    m_UpdateBuffer[1] = new UpdateListType;
  }
  ~LevelSetShiEvolutionBase()
  {
    delete m_UpdateBuffer[-1];
    delete m_UpdateBuffer[1];
  }

  unsigned int                m_NumberOfIterations;
  /// \todo is it useful?
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;

  // For sparse case, the update buffer needs to be the size of the active layer
  std::map< char, UpdateListType* >   m_UpdateBuffer;
  DomainMapImageFilterPointer         m_DomainMapFilter;

  LevelSetOutputRealType          m_Alpha;
  LevelSetOutputRealType          m_Dt;
  LevelSetOutputRealType          m_RMSChangeAccumulator;
  bool                            m_UserDefinedDt;

  void AllocateUpdateBuffer()
    {
//     m_UpdateBuffer[-1]->clear();
//     m_UpdateBuffer[1]->clear();
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

      this->InvokeEvent( IterationEvent() );
      }
    }


  void InitializeIteration()
  {
    std::cout << "Initialize iteration" << std::endl;
    DomainIteratorType map_it = m_DomainMapFilter->m_LevelSetMap.begin();
    DomainIteratorType map_end = m_DomainMapFilter->m_LevelSetMap.end();

    // Initialize parameters here
    m_EquationContainer->InitializeParameters();

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
  {}


  void ComputeDtForNextIteration()
    {}

  virtual void UpdateLevelSets()
    {
      std::cout << "Update levelsets" << std::endl;
      LevelSetPointer levelSet = m_LevelSetContainer->GetLevelSet( 0 );

      UpdateLevelSetFilterPointer update_levelset = UpdateLevelSetFilterType::New();
      update_levelset->SetSparseLevelSet( levelSet );
      update_levelset->SetEquationContainer( m_EquationContainer );
      update_levelset->Update();

      typedef ImageFileWriter< OutputImageType > WriterType;
      typedef typename WriterType::Pointer       WriterPointer;

      WriterPointer writer2 = WriterType::New();
      writer2->SetInput( levelSet->GetImage() );
      writer2->SetFileName("/home/krm15/2.mha");
      writer2->Update();

      m_RMSChangeAccumulator = update_levelset->GetRMSChangeAccumulator();

//       m_UpdateBuffer[-1]->clear();
//       m_UpdateBuffer[1]->clear();
    }

  void UpdateEquations()
    {
    std::cout << "Update equations" << std::endl << std::endl;
    m_EquationContainer->Update();
//     InitializeIteration();
    }

private:
  LevelSetShiEvolutionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkLevelSetShiEvolutionBase_h
