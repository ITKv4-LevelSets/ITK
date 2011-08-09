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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "itkLabelMapToLabelImageFilter.h"

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
  typedef typename LevelSetType::InputType             LevelSetInputType;
  typedef typename LevelSetType::OutputRealType        LevelSetOutputRealType;
  typedef typename LevelSetType::OutputType            LevelSetOutputType;

  typedef typename LevelSetType::LayerType             LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator         LevelSetLayerIterator;

  typedef typename LevelSetType::LabelMapType          LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer       LevelSetLabelMapPointer;

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

  typedef typename std::map< itk::IdentifierType, NounToBeDefined >::iterator DomainIteratorType;

  typedef UpdateShiSparseLevelSet< ImageDimension, EquationContainerType > UpdateLevelSetFilterType;

  typedef typename UpdateLevelSetFilterType::Pointer     UpdateLevelSetFilterPointer;

  itkSetObjectMacro( LevelSetContainer, LevelSetContainerType );
  itkGetObjectMacro( LevelSetContainer, LevelSetContainerType );

  void Update()
    {
    if( m_EquationContainer.IsNull() )
      {
      itkGenericExceptionMacro( << "m_EquationContainer is NULL" );
      }

    if( !m_EquationContainer->GetEquation( 0 ) )
      {
      itkGenericExceptionMacro( << "m_EquationContainer->GetEquation( 0 ) is NULL" );
      }

    m_DomainMapFilter = m_LevelSetContainer->GetDomainMapFilter();

    // Get the image to be segmented
    m_InputImage = m_EquationContainer->GetInput();

    if( m_InputImage.IsNull() )
      {
      itkGenericExceptionMacro( << "m_InputImage is NULL" );
      }

    TermContainerPointer Equation0 = m_EquationContainer->GetEquation( 0 );
    TermPointer term0 = Equation0->GetTerm( 0 );

    // Get the LevelSetContainer from the EquationContainer
    m_LevelSetContainer = term0->GetLevelSetContainer();

    if( term0.IsNull() )
      {
      itkGenericExceptionMacro( << "m_EquationContainer->GetEquation( 0 ) is NULL" );
      }

    if( !term0->GetLevelSetContainer() )
      {
      itkGenericExceptionMacro( << "m_LevelSetContainer is NULL" );
      }

    //Run iteration
    this->GenerateData();
    }

  itkSetMacro( Alpha, LevelSetOutputRealType );
  itkGetMacro( Alpha, LevelSetOutputRealType );

  void SetTimeStep( const LevelSetOutputRealType& )
    {
    // in this method, the time step is not used at all, only sign matters!
    // so m_Dt is constant and equal to 1 through all iteration
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
    m_Dt( 1. ), m_RMSChangeAccumulator( -1. )
  {
  }
  ~LevelSetShiEvolutionBase()
  {}

  unsigned int                m_NumberOfIterations;
  /// \todo is it useful?
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;

  DomainMapImageFilterPointer     m_DomainMapFilter;

  LevelSetOutputRealType          m_Alpha;
  LevelSetOutputRealType          m_Dt;
  LevelSetOutputRealType          m_RMSChangeAccumulator;

  void AllocateUpdateBuffer()
    {
    }


  void GenerateData()
    {
    AllocateUpdateBuffer();

    InitializeIteration();

    for( unsigned int iter = 0; iter < m_NumberOfIterations; iter++ )
      {
      std::cout <<"Iteration " <<iter << std::endl;
      m_RMSChangeAccumulator = NumericTraits< LevelSetOutputRealType >::Zero;

      // one iteration over all container
      // update each level set based on the different equations provided
      ComputeIteration();

      ComputeDtForNextIteration();

      UpdateLevelSets();
      UpdateEquations();

      // DEBUGGING
//       typedef Image< char, ImageDimension >     LabelImageType;
//       typedef typename LabelImageType::Pointer  LabelImagePointer;
//       typedef LabelMapToLabelImageFilter<LevelSetLabelMapType, LabelImageType> LabelMapToLabelImageFilterType;
//       typedef ImageFileWriter< LabelImageType > WriterType;
//
//       typename LevelSetContainerType::Iterator it = m_LevelSetContainer->Begin();
//       while( it != m_LevelSetContainer->End() )
//         {
//         std::ostringstream filename;
//         filename << "/home/krm15/temp/" << iter << "_" <<  it->GetIdentifier() << ".mha";
//
//         LevelSetPointer levelSet = it->GetLevelSet();
//
//         typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
//         labelMapToLabelImageFilter->SetInput( levelSet->GetLabelMap() );
//         labelMapToLabelImageFilter->Update();
//
//         typename WriterType::Pointer writer = WriterType::New();
//         writer->SetInput( labelMapToLabelImageFilter->GetOutput() );
//         writer->SetFileName( filename.str().c_str() );
//         writer->Update();
//
//         ++it;
//         }

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
    {
    std::cout << "Compute iteration" << std::endl;
    }


  void ComputeDtForNextIteration()
    {
    std::cout << "ComputeDtForNextIteration" << std::endl;
    }

  virtual void UpdateLevelSets()
    {
    std::cout << "Update levelsets" << std::endl;

    typename LevelSetContainerType::Iterator it = m_LevelSetContainer->Begin();
    bool singleLevelSet = ( m_LevelSetContainer->Size() == 1 );

    while( it != m_LevelSetContainer->End() )
      {
      LevelSetPointer levelSet = it->GetLevelSet();

      UpdateLevelSetFilterPointer update_levelset = UpdateLevelSetFilterType::New();
      update_levelset->SetInputLevelSet( levelSet );
      update_levelset->SetCurrentLevelSetId( it->GetIdentifier() );
      update_levelset->SetEquationContainer( m_EquationContainer );
      update_levelset->SetSingleLevelSet( singleLevelSet );
      update_levelset->Update();

      levelSet->Graft( update_levelset->GetOutputLevelSet() );

      m_RMSChangeAccumulator = update_levelset->GetRMSChangeAccumulator();

      ++it;
      }
    }

  void UpdateEquations()
    {
    std::cout << "Update equations" << std::endl << std::endl;
//     m_EquationContainer->Update();
    this->InitializeIteration();
    }

private:
  LevelSetShiEvolutionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkLevelSetShiEvolutionBase_h
