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


#ifndef __itkLevelSetWhitakerSparseEvolutionBase_h
#define __itkLevelSetWhitakerSparseEvolutionBase_h

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
class LevelSetWhitakerSparseEvolutionBase : public Object
{
public:
  typedef LevelSetWhitakerSparseEvolutionBase Self;
  typedef SmartPointer< Self >                Pointer;
  typedef SmartPointer< const Self >          ConstPointer;
  typedef Object                              Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetWhitakerSparseEvolutionBase, Object );

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

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

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
  typedef typename LevelSetType::InputType             LevelSetInputType;
  typedef typename LevelSetType::OutputType            LevelSetOutputType;
  typedef typename LevelSetImageType::Pointer          LevelSetImagePointer;

  typedef typename LevelSetType::NodeStatusType        LevelSetNodeStatusType;
  typedef typename LevelSetType::NodePairType          LevelSetNodePairType;
  typedef typename LevelSetType::NodeListType          LevelSetNodeListType;
  typedef typename LevelSetType::NodeListIterator      LevelSetNodeListIterator;
  typedef typename LevelSetType::NodeListConstIterator LevelSetNodeListConstIterator;

  typedef typename LevelSetType::SparseLayerMapType           SparseLayerMapType;
  typedef typename LevelSetType::SparseLayerMapIterator       SparseLayerMapIterator;
  typedef typename LevelSetType::SparseLayerMapConstIterator  SparseLayerMapConstIterator;


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

  itkSetMacro( Alpha, LevelSetOutputType );
  itkGetMacro( Alpha, LevelSetOutputType );

  void Update()
    {
    //Run iteration
    this->GenerateData();
    }

  void SetTimeStep( const LevelSetOutputType& iDt )
    {
    if( iDt > NumericTraits< LevelSetOutputType >::epsilon() )
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
  LevelSetWhitakerSparseEvolutionBase() : m_NumberOfIterations( 0 ), m_NumberOfLevelSets( 0 ),
    m_InputImage( NULL ), m_EquationContainer( NULL ), m_LevelSetContainer( NULL ),
    m_UpdateBuffer( NULL ), m_DomainMapFilter( NULL ), m_Alpha( 0.9 ),
    m_Dt( 1. ), m_RMSChangeAccumulator( -1. ), m_UserDefinedDt( false ),
    m_ConstantGradientValue( 1. )
  {
  // let's create here a city block neighborhood
  }

  ~LevelSetWhitakerSparseEvolutionBase() {}

  unsigned int                m_NumberOfIterations;
  /// \todo is it useful?
  unsigned int                m_NumberOfLevelSets;
  InputImagePointer           m_InputImage;
  EquationContainerPointer    m_EquationContainer;
  LevelSetContainerPointer    m_LevelSetContainer;
  DomainMapImageFilterPointer m_DomainMapFilter;

  LevelSetOutputType          m_Alpha;
  LevelSetOutputType          m_Dt;
  LevelSetOutputType          m_RMSChangeAccumulator;
  bool                        m_UserDefinedDt;

  LevelSetOutputType          m_ConstantGradientValue;

  std::map< unsigned int, std::list< InputPixelRealType > > m_UpdateBuffer;

  /** \brief Creates an empty list for each level-set, and stores it in
  m_UpdateBuffer */
  void AllocateUpdateBuffer()
    {
    // Kishore:
    // Make sure that the LevelSetContainer class works for SparseLevelSetBase
    // and LevelSetImageBase classes.
    // Write a small test to verify that!!

    // Arnaud: Since the update values (computed from each terms) are only
    // computed on the zero level-set, we can only keep a list of updates for
    // each level-set, right?
    LevelSetContainerIteratorType ls_it = m_LevelSetContainer->Begin();
    LevelSetContainerIteratorType ls_end = m_LevelSetContainer->End();

    while( ls_it != ls_end )
      {
      m_UpdateBuffer[ ls_it->first ] = std::list< InputPixelRealType >();
      ++ls_it;
      }
    }


  /** \brief Compute*/
  void ComputeIteration()
    {
    // first let's get the 0-list
    LevelSetContainerIteratorType ls_it = m_LevelSetContainer->Begin();
    LevelSetContainerIteratorType ls_end = m_LevelSetContainer->End();

    typename std::map< unsigned int, std::list< InputPixelRealType > >::iterator
        update_it = m_UpdateBuffer.begin();

    while( ls_it != ls_end )
      {
      LevelSetNodeListType* list_of_nodes = ( ls_it->second )->GetListNode( 0 );

      LevelSetNodeListIterator node_it = list_of_nodes->begin();
      LevelSetNodeListIterator node_end = list_of_nodes->end();

      LevelSetPointer LevelSetUpdate = m_UpdateBuffer->GetLevelSet( ls_it->first );

      update_it->second.clear();

      while( node_it != node_end )
        {
        LevelSetInputType idx = ( *node_it )->first;

        // get the contribution from all terms at idx
        InputPixelRealType temp_update =
            m_EquationContainer->GetEquation( ls_it->first )->Evaluate( idx );

        update_it->second.push_back( temp_update );

        ++node_it;
        }
      ++update_it;
      ++ls_it;
      }
    }


  // ---------------------------------------------------------------------------
  /** \brief Iterate over all equations/terms, and initialize each internal
  parameter.
  \note In the sparse case, parameters are computed only once; they are then
  updated through the evolution of the level set.
  */
  void InitializeIteration()
  {
    // Here we are setting the LevelSetFunction specific constants (lambda_1 etc) to 0
    // In the sparse case, the numerator and denominator always remain computed
    // So, we need to write member functions in LevelSetFunction specific to sparse cases
    // where we are only computing the ratio of numerator to denominator
    // In dense, we have to go back and compute the numerator and denominator
    // and then compute the ratio
    // In sparse, you know exactly the changes to be made. Let me see what was implemented in itkv3 for this function

  }


  // ---------------------------------------------------------------------------
  void GenerateData()
    {
    m_InputImage = m_EquationContainer->GetInput();

    // Get the LevelSetContainer from the EquationContainer
    m_LevelSetContainer =
        m_EquationContainer->GetEquation( 0 )->GetTerm( 0 )->GetLevelSetContainer();

    // allocate an update buffer to store the update values
    AllocateUpdateBuffer();

    InitializeIteration();

    for( unsigned int iter = 0; iter < m_NumberOfIterations; ++iter )
      {
      m_RMSChangeAccumulator = NumericTraits< LevelSetOutputType >::Zero;

      // one iteration over all container
      // update each level set based on the different equations provided
      ComputeIteration();

      ComputeDtForNextIteration();

      UpdateLevelSets();

      Reinitialize();

      UpdateEquations();

      this->InvokeEvent( IterationEvent() );
      }
    }

  // ---------------------------------------------------------------------------
  void ComputeDtForNextIteration()
    {
    if( !m_UserDefinedDt )
      {
      if( ( m_Alpha > NumericTraits< LevelSetOutputType >::Zero ) &&
          ( m_Alpha < NumericTraits< LevelSetOutputType >::One ) )
        {
        LevelSetOutputType contribution = m_EquationContainer->GetCFLContribution();

        if( contribution > NumericTraits< LevelSetOutputType >::epsilon() )
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
    }

  // ---------------------------------------------------------------------------
  void UpdateAllLayers()
    {
    LevelSetContainerIteratorType ls_it = m_LevelSetContainer->Begin();
    LevelSetContainerIteratorType ls_end = m_LevelSetContainer->End();

    // iterate on all level set
    while( ls_it != ls_end )
      {
      // get the active layer (zero level-set)
      LevelSetNodeListType* list_of_nodes = ( ls_it->second )->GetListNode( 0 );

      // create a 2 temporary lists to keep track of what is going up, and what
      // is going down
      LevelSetNodeListType up_list, down_list;

      // update the active layer
      UpdateActiveLayer( list_of_nodes, &up_list, &down_list );


      ++ls_it;
      }
    }

  void UpdateActiveLayer( const std::list< InputPixelRealType >& iUpdates,
    LevelSetNodeListType* iCurrentList,
    LevelSetNodeListType* ioUpList,
    LevelSetNodeListType* ioDownList )
    {
    const LevelSetOutputType LowerActiveThreshold = - 0.5 * m_ConstantGradientValue;
    const LevelSetOutputType UpperActiveThreshold = 0.5 * m_ConstantGradientValue;

    LevelSetNodeListIterator node_it = iCurrentList->begin();
    LevelSetNodeListIterator node_end = iCurrentList->end();

    typename std::list< InputPixelRealType >::const_iterator
        update_it = iUpdates.begin();

    NeighborhoodIterator< LevelSetImageType >
      n_It( m_NeighborList.GetRadius(),
            ( ls_it->second )->GetImage(), // get the resulting level set image
            this->GetOutput()->GetRequestedRegion() );

    while( node_it != node_end )
      {
      current_index = ( *node_it )->first;
      current_attributes = ( *node_it )->second;

      n_It.SetLocation( current_index );

      // get the update for the given level set and the given pixel
      InputPixelRealType temp_update = *update_it;

      // initialize new_value to the current level set value
      LevelSetOutputType new_value = attributes.m_Value;
        new_value += m_Dt * temp_update.m_Value;

      if( new_value >= UpperActiveThreshold )
        {
        // This index will move UP into a positive (outside) layer.
        // First check for active layer neighbors moving in the opposite
        // direction.
        flag = false;

        for( unsigned int i = 0; i < m_NeighborList.GetSize(); i++ )
          {
          neigh_attributes = n_It.GetPixel( m_NeighborList.GetArrayIndex(i) );

          if( neigh_attributes.m_Status == m_StatusChangingDown )
            {
            flag = true;
            break;
            }
          }
        if( flag )
          {
          ++node_it;
          ++update_it;
          }
        else
          {
          rms_change_accumulator +=
              vnl_math_sqr( new_value - outputIt.GetCenterPixel() );

          // Search the neighborhood for inside indicies.
          LevelSetOutputType temp_value = new_value - m_ConstantGradientValue;

          for ( unsigned int i = 0; i < m_NeighborList.GetSize(); ++i )
            {
            internal_idx = m_NeighborList.GetArrayIndex(i);
            neigh_attributes = n_It.GetPixel(idx);

            if ( neigh_attributes.m_Status == 1 )
              {
              // Keep the smallest possible value for the new active node.  This
              // places the new active layer node closest to the zero level-set.
              LevelSetOutputType output_value = neigh_attributes.m_Value;

              if( ( output_value < LowerActiveThreshold ) ||
                    ( vnl_math_abs(temp_value) < vnl_math_abs( output_value ) ) )
                {
                neigh_attributes.m_Value = temp_value;

                n_It.SetPixel( idx, neigh_attributes, bounds_status );
                }
              }
            }

          current_attributes.m_Status = m_StatusActiveChangingUp;
          current_attributes.m_Value = new_value;

          NodePairType node( current_index, attributes );

          ioUpList->push_front( node );

          release_node = layerIt.GetPointer();
          node_it = iCurrentList->erase( node_it );
          }
        }
      ++node_it;
      }
    }


  void
  UpdateStatus(
    LevelSetPointer iLevelSet,
    LevelSetNodeListType *InputList, LevelSetNodeListType *OutputList,
    LevelSetNodeStatusType iChangeToStatus, LevelSetNodeStatusType iSearchForStatus)
  {
  NeighborhoodIterator< LevelSetImageType >
  n_It( m_NeighborList.GetRadius(),
           ( ls_it->second )->GetImage(), // get the resulting level set image
           this->GetOutput()->GetRequestedRegion() );

  if ( !m_BoundsCheckingActive )
    {
    statusIt.NeedToUseBoundaryConditionOff();
    }

  // Push each index in the input list into its appropriate status layer
  // (ChangeToStatus) and update the status image value at that index.
  // Also examine the neighbors of the index to determine which need to go onto
  // the output list (search for SearchForStatus).
  while ( !InputList->empty() )
    {
    NodePairType node_pair = InputList->front();

    IndexType current_index = node_pair.first;
    n_It.SetLocation( current_index );

    node_pair.m_Status = iChangeToStatus;

    InputList->pop_front();      // _before_ transferring to another list.

    iLevelSet->GetListNode( iChangeToStatus )->push_front( node_pair );

    for ( unsigned int i = 0; i < m_NeighborList.GetSize(); ++i )
      {
      // get the status of the neighbor
      idx = m_NeighborList.GetArrayIndex(i);
      StatusType neighbor_status = statusIt.GetPixel( idx );

      // Have we bumped up against the boundary?  If so, turn on bounds
      // checking.
      if ( neighbor_status == m_StatusBoundaryPixel )
        {
        m_BoundsCheckingActive = true;
        }

      if ( neighbor_status == iSearchForStatus )
        {
        bool bounds_status = false;

        // mark this pixel so we don't add it twice.
        statusIt.SetPixel( idx, m_StatusChanging, bounds_status );
        if ( bounds_status == true )
          {
          node = m_LayerNodeStore->Borrow();
          node->m_Value = statusIt.GetIndex()
                          + m_NeighborList.GetNeighborhoodOffset(i);
          OutputList->PushFront(node);
          } // else this index was out of bounds.
        }
      }
    }
  }

  virtual void UpdateLevelSets()
    {
//    LevelSetContainerIteratorType it1 = m_LevelSetContainer->Begin();
//    LevelSetContainerConstIteratorType it2 = m_UpdateBuffer->Begin();

//    LevelSetOutputType p;

//    while( it1 != m_LevelSetContainer->End() )
//      {
//      LevelSetImagePointer image1 = it1->second->GetImage();
//      LevelSetImagePointer image2 = it2->second->GetImage();

//      LevelSetImageIteratorType It1( image1, image1->GetBufferedRegion() );
//      LevelSetImageIteratorType It2( image2, image2->GetBufferedRegion() );
//      It1.GoToBegin();
//      It2.GoToBegin();

//      while( !It1.IsAtEnd() )
//        {
//        p = m_Dt * It2.Get();
//        It1.Set( It1.Get() + p );

//        m_RMSChangeAccumulator += p*p;

//        ++It1;
//        ++It2;
//        }

//      ++it1;
//      ++it2;
//      }
  }

  void UpdateEquations()
    {
//    InitializeIteration();
//     m_EquationContainer->Update();
    }

  void Reinitialize()
  {
//    LevelSetContainerIteratorType it = m_LevelSetContainer->Begin();

//    while( it != m_LevelSetContainer->End() )
//      {
//      LevelSetImagePointer image = it->second->GetImage();

//      ThresholdFilterPointer thresh = ThresholdFilterType::New();
//      thresh->SetLowerThreshold(
//            NumericTraits< LevelSetOutputType >::NonpositiveMin() );
//      thresh->SetUpperThreshold( 0 );
//      thresh->SetInsideValue( 1 );
//      thresh->SetOutsideValue( 0 );
//      thresh->SetInput( image );
//      thresh->Update();

//      MaurerPointer maurer = MaurerType::New();
//      maurer->SetInput( thresh->GetOutput() );
//      maurer->SetSquaredDistance( false );
//      maurer->SetUseImageSpacing( true );
//      maurer->SetInsideIsPositive( false );
//      maurer->Update();

//      LevelSetImageIteratorType It1( image, image->GetBufferedRegion() );
//      LevelSetImageIteratorType It2( maurer->GetOutput(), image->GetBufferedRegion() );
//      It1.GoToBegin();
//      It2.GoToBegin();
//      while( !It1.IsAtEnd() )
//        {
//        It1.Set( It2.Get() );
//        ++It1;
//        ++It2;
//        }
//      ++it;
//      }
  }

private:
  LevelSetWhitakerSparseEvolutionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkLevelSetWhitakerSparseEvolutionBase_h
