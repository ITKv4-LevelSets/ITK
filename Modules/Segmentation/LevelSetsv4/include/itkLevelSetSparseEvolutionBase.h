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
#include "itkLevelSetEvolutionStoppingCriterionBase.h"
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
  typedef typename LevelSetContainerType::LevelSetIdentifierType
                                                            LevelSetIdentifierType;
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

  //   typedef typename DomainMapImageFilterType::DomainIteratorType DomainIteratorType;
  typedef typename std::map< itk::IdentifierType, NounToBeDefined >::iterator DomainIteratorType;

  typedef UpdateWhitakerSparseLevelSet< ImageDimension, LevelSetOutputType, EquationContainerType >
                                                  UpdateLevelSetFilterType;
//  typedef typename LevelSetType::StatusImageType  StatusImageType;
//  typedef typename LevelSetType::OutputImageType  OutputImageType;

  typedef typename UpdateLevelSetFilterType::Pointer                         UpdateLevelSetFilterPointer;
//  typedef typename UpdateLevelSetFilterType::UpdateListType                  UpdateListType;

  typedef LevelSetEvolutionStoppingCriterionBase< LevelSetContainerType >
                                                  StoppingCriterionType;
  typedef typename StoppingCriterionType::Pointer StoppingCriterionPointer;

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

    // Get the LevelSetContainer from the EquationContainer
    m_LevelSetContainer = m_EquationContainer->GetEquation( 0 )->GetTerm( 0 )->GetLevelSetContainer();

    if( m_StoppingCriterion.IsNull() )
      {
      itkGenericExceptionMacro( << "m_StoppingCriterion is NULL" );
      }

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

  /** \brief Set/Get the Stopping Criterion */
  itkGetObjectMacro( StoppingCriterion, StoppingCriterionType );
  itkSetObjectMacro( StoppingCriterion, StoppingCriterionType );

  // set the domain map image filter
  itkSetObjectMacro( DomainMapFilter, DomainMapImageFilterType );
  itkGetObjectMacro( DomainMapFilter, DomainMapImageFilterType );

protected:
  LevelSetSparseEvolutionBase() : m_StoppingCriterion( NULL ), m_NumberOfLevelSets( 0 ),
    m_InputImage( NULL ), m_EquationContainer( NULL ), m_LevelSetContainer( NULL ),
    m_DomainMapFilter( NULL ), m_Alpha( 0.9 ),
    m_Dt( 1. ), m_RMSChangeAccumulator( -1. ), m_UserDefinedDt( false )
  {
  }
  ~LevelSetSparseEvolutionBase()
  {
    typename LevelSetContainerType::ConstIterator it = m_LevelSetContainer->Begin();
    while( it != m_LevelSetContainer->End() )
      {
      delete m_UpdateBuffer[ it->GetIdentifier() ];
      ++it;
      }
  }


  StoppingCriterionPointer  m_StoppingCriterion;

  /// \todo is it useful?
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;

  // For sparse case, the update buffer needs to be the size of the active layer
  std::map< IdentifierType, LevelSetLayerType* >  m_UpdateBuffer;
  DomainMapImageFilterPointer                     m_DomainMapFilter;

  LevelSetOutputRealType      m_Alpha;
  LevelSetOutputRealType      m_Dt;
  LevelSetOutputRealType      m_RMSChangeAccumulator;
  bool                        m_UserDefinedDt;

  void AllocateUpdateBuffer()
    {
    typename LevelSetContainerType::Iterator it = m_LevelSetContainer->Begin();
    while( it != m_LevelSetContainer->End() )
      {
      IdentifierType id = it->GetIdentifier();

      if( m_UpdateBuffer.find( id ) == m_UpdateBuffer.end() )
        {
        m_UpdateBuffer[ id ] = new LevelSetLayerType;
        }
      else
        {
        if( m_UpdateBuffer[ id ] )
          {
          m_UpdateBuffer[ id ]->clear();
          }
        else
          {
          m_UpdateBuffer[ id ] = new LevelSetLayerType;
          }
        }
      ++it;
      }
    }

  void GenerateData()
    {
    AllocateUpdateBuffer();

    InitializeIteration();

    typename StoppingCriterionType::IterationIdType iter = 0;
    m_StoppingCriterion->SetCurrentIteration( iter );
    m_StoppingCriterion->SetLevelSetContainer( m_LevelSetContainer );

//    for( unsigned int iter = 0; iter < m_NumberOfIterations; iter++ )
    while( !m_StoppingCriterion->IsSatisfied() )
      {
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

      ++iter;

      m_StoppingCriterion->SetRMSChangeAccumulator( m_RMSChangeAccumulator );
      m_StoppingCriterion->SetCurrentIteration( iter );
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
    typename LevelSetContainerType::Iterator it = m_LevelSetContainer->Begin();

    while( it != m_LevelSetContainer->End() )
      {
      LevelSetPointer levelSet = it->GetLevelSet();
      LevelSetLayerIterator list_it = levelSet->GetLayer( 0 ).begin();
      LevelSetLayerIterator list_end = levelSet->GetLayer( 0 ).end();

      while( list_it != list_end )
        {
        LevelSetInputType idx = list_it->first;

        InputPixelRealType temp_update =
            m_EquationContainer->GetEquation( it->GetIdentifier() )->Evaluate( idx );
        m_UpdateBuffer[ it->GetIdentifier() ]->insert(
              std::pair< LevelSetInputType, LevelSetOutputType >( idx, temp_update ) );
//         std::cout << idx << ' ' << temp_update << std::endl;
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

      typename LevelSetContainerType::Iterator it = m_LevelSetContainer->Begin();
      while( it != m_LevelSetContainer->End() )
      {
//         std::cout << "m_UpdateBuffer size=" << m_UpdateBuffer[it->GetIdentifier()]->size() << std::endl;
//         std::cout << "Layer -2 =" << it->GetLevelSet()->GetLayer( -2 ).size() << std::endl;
//         std::cout << "Layer -1 =" << it->GetLevelSet()->GetLayer( -1 ).size() << std::endl;
//         std::cout << "Layer  0 =" << it->GetLevelSet()->GetLayer( 0 ).size() << std::endl;
//         std::cout << "Layer  1 =" << it->GetLevelSet()->GetLayer( 1 ).size() << std::endl;
//         std::cout << "Layer  2 =" << it->GetLevelSet()->GetLayer( 2 ).size() << std::endl << std::endl;

        LevelSetPointer levelSet = it->GetLevelSet();

        UpdateLevelSetFilterPointer update_levelset = UpdateLevelSetFilterType::New();
        update_levelset->SetInputLevelSet( levelSet );
        update_levelset->SetUpdate( * m_UpdateBuffer[it->GetIdentifier()] );
        update_levelset->SetEquationContainer( m_EquationContainer );
        update_levelset->SetDt( m_Dt );
        update_levelset->Update();

//         std::cout << "Layer -2 =" << update_levelset->GetOutputLevelSet()->GetLayer( -2 ).size() << std::endl;
//         std::cout << "Layer -1 =" << update_levelset->GetOutputLevelSet()->GetLayer( -1 ).size() << std::endl;
//         std::cout << "Layer  0 =" << update_levelset->GetOutputLevelSet()->GetLayer( 0 ).size() << std::endl;
//         std::cout << "Layer  1 =" << update_levelset->GetOutputLevelSet()->GetLayer( 1 ).size() << std::endl;
//         std::cout << "Layer  2 =" << update_levelset->GetOutputLevelSet()->GetLayer( 2 ).size() << std::endl;

        levelSet->Graft( update_levelset->GetOutputLevelSet() );

        m_RMSChangeAccumulator = update_levelset->GetRMSChangeAccumulator();

        m_UpdateBuffer[it->GetIdentifier()]->clear();
        ++it;
      }
    }

  void UpdateEquations()
    {
    std::cout << "Update equations" << std::endl;
//     m_EquationContainer->Update();
    InitializeIteration();
    }

private:
  LevelSetSparseEvolutionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkLevelSetSparseEvolutionBase_h
