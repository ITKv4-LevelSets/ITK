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

#ifndef __itkUpdateWhitakerSparseLevelSet_h
#define __itkUpdateWhitakerSparseLevelSet_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkWhitakerSparseLevelSetBase.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include <list>
#include "itkObject.h"


#include "itkImageFileWriter.h"

namespace itk
{
template< unsigned int VDimension,
          typename TLevelSetValueType,
          class TEquationContainer >
class UpdateWhitakerSparseLevelSet : public Object
{
public:
  typedef UpdateWhitakerSparseLevelSet  Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( UpdateWhitakerSparseLevelSet, Object );

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef TLevelSetValueType  LevelSetOutputType;

  typedef WhitakerSparseLevelSetBase< LevelSetOutputType, ImageDimension >
                                                       LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::InputType             LevelSetInputType;

  typedef typename LevelSetType::LabelMapType          LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer       LevelSetLabelMapPointer;

  typedef typename LevelSetType::LabelObjectType       LevelSetLabelObjectType;
  typedef typename LevelSetType::LabelObjectPointer    LevelSetLabelObjectPointer;
  typedef typename LevelSetType::LabelObjectLengthType LevelSetLabelObjectLengthType;
  typedef typename LevelSetType::LabelObjectLineType   LevelSetLabelObjectLineType;
  typedef typename LevelSetType::LabelObjectLineContainerType
                                                       LevelSetLabelObjectLineContainerType;

  typedef typename LevelSetType::LayerType             LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator         LevelSetLayerIterator;
  typedef typename LevelSetType::LayerConstIterator    LevelSetLayerConstIterator;
  typedef typename LevelSetType::OutputRealType        LevelSetOutputRealType;

  typedef typename LevelSetType::LayerMapType           LevelSetLayerMapType;
  typedef typename LevelSetType::LayerMapIterator       LevelSetLayerMapIterator;
  typedef typename LevelSetType::LayerMapConstIterator  LevelSetLayerMapConstIterator;

  typedef TEquationContainer                      EquationContainerType;
  typedef typename EquationContainerType::Pointer EquationContainerPointer;

  itkGetObjectMacro( OutputLevelSet, LevelSetType );

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void Update()
  {
    if( m_SparseLevelSet.IsNull() )
      {
      itkGenericExceptionMacro( <<"m_SparseLevelSet is NULL" );
      }
    if( m_Update.empty() )
      {
      itkGenericExceptionMacro( <<"m_Update is empty" );
      }

    m_OutputLevelSet->SetLayer( -2, m_SparseLevelSet->GetLayer( -2 ) );
    m_OutputLevelSet->SetLayer( -1, m_SparseLevelSet->GetLayer( -1 ) );
    m_OutputLevelSet->SetLayer(  0, m_SparseLevelSet->GetLayer(  0 ) );
    m_OutputLevelSet->SetLayer(  1, m_SparseLevelSet->GetLayer(  1 ) );
    m_OutputLevelSet->SetLayer(  2, m_SparseLevelSet->GetLayer(  2 ) );

    m_OutputLevelSet->SetLabelMap( m_SparseLevelSet->GetLabelMap() );

    typedef LabelMapToLabelImageFilter<LevelSetLabelMapType, LabelImageType> LabelMapToLabelImageFilterType;
    typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput( m_SparseLevelSet->GetLabelMap() );
    labelMapToLabelImageFilter->Update();

    m_InternalImage = labelMapToLabelImageFilter->GetOutput();
    m_InternalImage->DisconnectPipeline();

    m_TempPhi.clear();

    for( char status = -1; status < 2; status++ )
      {
      LevelSetLayerType layer = m_SparseLevelSet->GetLayer( status );

      LevelSetLayerConstIterator it = layer.begin();
      while( it != layer.end() )
        {
        m_TempPhi[ it->first ] = it->second;
        ++it;
        }
      }

    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType layerMinus2 = m_SparseLevelSet->GetLayer( -2 );

    LevelSetLayerConstIterator it = layerMinus2.begin();
    while( it != layerMinus2.end() )
      {
      LevelSetInputType currentIndex = it->first;
      m_TempPhi[ currentIndex ] = -2;
      neighIt.SetLocation( currentIndex );

      for( typename NeighborhoodIteratorType::Iterator nIt = neighIt.Begin();
           !nIt.IsAtEnd();
           ++nIt )
        {
        if( nIt.Get() == -3 )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( nIt.GetNeighborhoodOffset() );

          m_TempPhi[ tempIndex ] = -3;
          }
        }

      ++it;
      }

    LevelSetLayerType layerPlus2 = m_SparseLevelSet->GetLayer( 2 );

    it = layerPlus2.begin();
    while( it != layerPlus2.end() )
      {
      LevelSetInputType currentIndex = it->first;
      m_TempPhi[ currentIndex ] = 2;
      neighIt.SetLocation( currentIndex );

      for( typename NeighborhoodIteratorType::Iterator nIt = neighIt.Begin();
           !nIt.IsAtEnd();
           ++nIt )
        {
        if( nIt.Get() == 3 )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( nIt.GetNeighborhoodOffset() );

          m_TempPhi[ tempIndex ] = 3;
          }
        }

      ++it;
      }

    this->UpdateZeroLevelSet();

    this->UpdateLminus1();

    this->UpdateLplus1();

    this->UpdateLminus2();

    this->UpdateLplus2();

    this->MovePointIntoZeroLevelSet();

    this->MovePointFromMinus1();

    this->MovePointFromPlus1();

    this->MovePointFromMinus2();

    this->MovePointFromPlus2();

    typedef ImageFileWriter< LabelImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( m_InternalImage );
    writer->SetFileName( "internal_image.mha" );
    writer->Update();

    typedef LabelImageToLabelMapFilter< LabelImageType, LevelSetLabelMapType> LabelImageToLabelMapFilterType;
    typename LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();
    labelImageToLabelMapFilter->SetInput( m_InternalImage );
    labelImageToLabelMapFilter->SetBackgroundValue( 3 );
    labelImageToLabelMapFilter->Update();

    m_OutputLevelSet->GetLabelMap( )->Graft( labelImageToLabelMapFilter->GetOutput() );
    m_TempPhi.clear();
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set/Get the Dt for the update
  itkSetMacro( Dt, LevelSetOutputType );
  itkGetMacro( Dt, LevelSetOutputType );

  // Set/Get the RMS change for the update
  itkGetMacro( RMSChangeAccumulator, LevelSetOutputType );

  itkSetObjectMacro( EquationContainer, EquationContainerType );
  itkGetObjectMacro( EquationContainer, EquationContainerType );

  void SetUpdate( const LevelSetLayerType& iUpdate )
    {
    m_Update = iUpdate;
    }

protected:
  UpdateWhitakerSparseLevelSet() : m_Dt( NumericTraits< LevelSetOutputType >::One ),
    m_RMSChangeAccumulator( NumericTraits< LevelSetOutputType >::Zero ),
    m_MinStatus( -3 ),  m_MaxStatus( 3 )
    {
    m_TempLevelSet = LevelSetType::New();
    m_OutputLevelSet = LevelSetType::New();
    }
  ~UpdateWhitakerSparseLevelSet() {}

  LevelSetOutputType m_Dt;
  LevelSetOutputType m_RMSChangeAccumulator;

  EquationContainerPointer m_EquationContainer;

  LevelSetLayerType  m_Update;
  LevelSetPointer    m_SparseLevelSet;
  LevelSetPointer    m_OutputLevelSet;

  LevelSetLayerType m_TempPhi;

  LevelSetPointer   m_TempLevelSet;
  char              m_MinStatus;
  char              m_MaxStatus;

  typedef Image< char, ImageDimension >     LabelImageType;
  typedef typename LabelImageType::Pointer  LabelImagePointer;

  LabelImagePointer m_InternalImage;

  typedef ShapedNeighborhoodIterator< LabelImageType > NeighborhoodIteratorType;

  void UpdateZeroLevelSet()
  {
    LevelSetLayerType& layer0 = m_OutputLevelSet->GetLayer( static_cast< char >( 0 ) );

    assert( m_Update.size() == layer0.size() );

    LevelSetLayerIterator nodeIt   = layer0.begin();
    LevelSetLayerIterator nodeEnd  = layer0.end();

    LevelSetLayerIterator upIt     = m_Update.begin();

    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    while( nodeIt != nodeEnd )
      {
      assert( nodeIt->first == upIt->first );

      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;
      LevelSetOutputType  tempUpdate =
          m_Dt * static_cast< LevelSetOutputType >( upIt->second );

//      if( tempUpdate > 0.5 )
//        {
//        tempUpdate = 0.5;
//        }
//      else if( tempUpdate < - 0.5 )
//        {
//        tempUpdate = - 0.5;
//        }

      LevelSetOutputType tempValue = currentValue + tempUpdate;
      m_RMSChangeAccumulator += tempUpdate*tempUpdate;

      if( tempValue > 0.5 )
        {
        // is there any point moving in the opposite direction?
        bool samedirection = true;

        neighIt.SetLocation( currentIndex );

        for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
             !it.IsAtEnd();
             ++it )
          {
          if( it.Get() == 0 )
            {
            LevelSetInputType tempIndex =
                neighIt.GetIndex( it.GetNeighborhoodOffset() );

            LevelSetLayerIterator it = m_TempPhi.find( tempIndex );
            if( it != m_TempPhi.end() )
              {
              if( it->second < -0.5 )
                {
                samedirection = false;
                }
              }
            }
          }

        if( samedirection )
          {
          LevelSetLayerIterator it = m_TempPhi.find( currentIndex );

          if( it != m_TempPhi.end() )
            {
            it->second = tempValue;
            }
          else
            {
            m_TempPhi.insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                      tempValue ) );
            }

          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          ++upIt;
          layer0.erase( tempIt );

          m_TempLevelSet->GetLayer( 1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                    tempValue ) );

//          m_InternalImage->SetPixel( currentIndex, static_cast< char >( 1 ) );
          }
        else // samedirection == false
          {
//          m_InternalImage->SetPixel( currentIndex, static_cast< char >( 0 ) );
          ++nodeIt;
          ++upIt;
          }
        } // end of if( tempValue > 0.5 )

      else // ( tempValue <= 0.5 )
        {
        if( tempValue < -0.5 )
          {
          bool samedirection = true;

          neighIt.SetLocation( currentIndex );

          for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
               !it.IsAtEnd();
               ++it )
            {
            if( it.Get() == 0 )
              {
              LevelSetInputType tempIndex =
                  neighIt.GetIndex( it.GetNeighborhoodOffset() );

              LevelSetLayerIterator it = m_TempPhi.find( tempIndex );
              if( it != m_TempPhi.end() )
                {
                if( it->second > 0.5 )
                  {
                  samedirection = false;
                  }
                }
              }
            }

          if( samedirection )
            {
            LevelSetLayerIterator it = m_TempPhi.find( currentIndex );

            if( it != m_TempPhi.end() )
              {
              it->second = tempValue;
              }
            else
              {
              m_TempPhi.insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                        tempValue ) );
              }

            LevelSetLayerIterator tempIt = nodeIt;
            ++nodeIt;
            ++upIt;

            layer0.erase( tempIt );

            m_TempLevelSet->GetLayer( -1 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                      tempValue ) );

//            m_InternalImage->SetPixel( currentIndex, static_cast< char >( -1 ) );
            }
          else // samedirection == false
            {
//            m_InternalImage->SetPixel( currentIndex, static_cast< char >( 0 ) );
            ++nodeIt;
            ++upIt;
            }
          }
        else // -0.5 <= temp <= 0.5
          {
//          m_InternalImage->SetPixel( currentIndex, static_cast< char >( 0 ) );
          ++nodeIt;
          ++upIt;
          }
        }
      } // while( nodeIt != nodeEnd )
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void UpdateLminus1()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType& LayerMinus1 = m_OutputLevelSet->GetLayer( static_cast< char >( -1 ) );

    LevelSetLayerIterator nodeIt   = LayerMinus1.begin();
    LevelSetLayerIterator nodeEnd  = LayerMinus1.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool IsThereAPointWithLabelEqualTo0 = false;

      LevelSetOutputType M = NumericTraits< LevelSetOutputType >::NonpositiveMin();

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        char label = it.Get();

        if( label >= 0 )
          {
          if( label == 0 )
            {
            IsThereAPointWithLabelEqualTo0 = true;
            }
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
          assert( phiIt != m_TempPhi.end() );

          M = vnl_math_max( M, phiIt->second );
          }
        } // end for

      if( IsThereAPointWithLabelEqualTo0 )
        {
        LevelSetLayerIterator phiIt = m_TempPhi.find( currentIndex );

        M = M - 1.;

        if( phiIt != m_TempPhi.end() )
          {
          phiIt->second = M;
          }
        else
          {
          m_TempPhi.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M ) );
          }

        if( M >= -0.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerMinus1.erase( tempIt );
          m_TempLevelSet->GetLayer( 0 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else if( M < -1.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerMinus1.erase( tempIt );
          m_TempLevelSet->GetLayer( -2 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else
          {
          ++nodeIt;
          }
        }
      else
        {
        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        LayerMinus1.erase( tempIt );
        m_TempLevelSet->GetLayer( -2 ).insert(
              std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                  nodeIt->second) );
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void UpdateLplus1()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType& LayerPlus1 = m_OutputLevelSet->GetLayer( static_cast< char >( 1 ) );

    LevelSetLayerIterator nodeIt   = LayerPlus1.begin();
    LevelSetLayerIterator nodeEnd  = LayerPlus1.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool IsThereAPointWithLabelEqualTo0 = false;

      LevelSetOutputType M = NumericTraits< LevelSetOutputType >::max();

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        char label = it.Get();

        if( label <= 0 )
          {
          if( label == 0 )
            {
            IsThereAPointWithLabelEqualTo0 = true;
            }
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
          if( phiIt != m_TempPhi.end() )
            {
            M = vnl_math_min( M, phiIt->second );
            }
          else
            {
            std::cout << tempIndex << "is not in m_TempPhi"<< std::endl;
            }
          }
        } // end for

      if( IsThereAPointWithLabelEqualTo0 )
        {
        LevelSetLayerIterator phiIt = m_TempPhi.find( currentIndex );

        M = M + 1.;

        if( phiIt != m_TempPhi.end() )
          {
          phiIt->second = M;
          }
        else
          {
          m_TempPhi.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M ) );
          }

        if( M <= 0.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerPlus1.erase( tempIt );
          m_TempLevelSet->GetLayer( 0 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else if( M > 1.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerPlus1.erase( tempIt );
          m_TempLevelSet->GetLayer( 2 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else
          {
          ++nodeIt;
          }
        }
      else
        {
        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        LayerPlus1.erase( tempIt );
        m_TempLevelSet->GetLayer( 2 ).insert(
              std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                  nodeIt->second) );
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void UpdateLminus2()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType& LayerMinus2 = m_OutputLevelSet->GetLayer( static_cast< char >( -2 ) );

    LevelSetLayerIterator nodeIt   = LayerMinus2.begin();
    LevelSetLayerIterator nodeEnd  = LayerMinus2.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool IsThereAPointWithLabelEqualToMinus1 = false;

      LevelSetOutputType M = NumericTraits< LevelSetOutputType >::NonpositiveMin();

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        char label = it.Get();

        if( label >= -1 )
          {
          if( label == -1 )
            {
            IsThereAPointWithLabelEqualToMinus1 = true;
            }
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
          assert( phiIt != m_TempPhi.end() );

          M = vnl_math_max( M, phiIt->second );
          }
        } // end for

      if( IsThereAPointWithLabelEqualToMinus1 )
        {
        LevelSetLayerIterator phiIt = m_TempPhi.find( currentIndex );

        M = M - 1.;

        if( phiIt != m_TempPhi.end() )
          {
          phiIt->second = M;
          }
        else
          {
          m_TempPhi.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M ) );
          }

        if( M >= -1.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerMinus2.erase( tempIt );
          m_TempLevelSet->GetLayer( -1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else if( M < -2.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerMinus2.erase( tempIt );
          m_InternalImage->SetPixel( currentIndex, -3 );
          m_TempPhi.erase( currentIndex );
          }
        else
          {
          ++nodeIt;
          }
        }
      else
        {
        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        LayerMinus2.erase( tempIt );
        m_InternalImage->SetPixel( currentIndex, -3 );
        m_TempPhi.erase( currentIndex );
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void UpdateLplus2()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType& LayerPlus2 = m_OutputLevelSet->GetLayer( static_cast< char >( 2 ) );

    LevelSetLayerIterator nodeIt   = LayerPlus2.begin();
    LevelSetLayerIterator nodeEnd  = LayerPlus2.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool IsThereAPointWithLabelEqualToPlus1 = false;

      LevelSetOutputType M = NumericTraits< LevelSetOutputType >::max();

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        char label = it.Get();

        if( label <= 1 )
          {
          if( label == 1 )
            {
            IsThereAPointWithLabelEqualToPlus1 = true;
            }
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
          if( phiIt != m_TempPhi.end() )
            {
            M = vnl_math_min( M, phiIt->second );
            }
          else
            {
            std::cout << tempIndex << " is not in m_TempPhi" <<std::endl;
            }
          }
        } // end for

      if( IsThereAPointWithLabelEqualToPlus1 )
        {
        LevelSetLayerIterator phiIt = m_TempPhi.find( currentIndex );

        M = M + 1.;

        if( phiIt != m_TempPhi.end() )
          {
          phiIt->second = M;
          }
        else
          {
          m_TempPhi.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M ) );
          }

        if( M <= 1.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerPlus2.erase( tempIt );
          m_TempLevelSet->GetLayer( 1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, M) );
          }
        else if( M > 2.5 )
          {
          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          LayerPlus2.erase( tempIt );
          m_InternalImage->SetPixel( currentIndex, 3 );
          m_TempPhi.erase( currentIndex );
          }
        else
          {
          ++nodeIt;
          }
        }
      else
        {
        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        LayerPlus2.erase( tempIt );
        m_InternalImage->SetPixel( currentIndex, 3 );
        m_TempPhi.erase( currentIndex );
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void MovePointIntoZeroLevelSet()
    {
    LevelSetLayerType layer0 = m_TempLevelSet->GetLayer( static_cast< char >( 0 ) );

    LevelSetLayerIterator nodeIt = layer0.begin();
    LevelSetLayerIterator nodeEnd = layer0.end();

    while( nodeIt != nodeEnd )
      {
      m_OutputLevelSet->GetLayer( 0 ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( nodeIt->first, nodeIt->second ) );
      m_InternalImage->SetPixel( nodeIt->first, 0 );
      LevelSetLayerIterator tempIt = nodeIt;
      ++nodeIt;
      layer0.erase( tempIt );
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void MovePointFromMinus1()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType layerMinus1 = m_TempLevelSet->GetLayer( static_cast< char >( -1 ) );

    LevelSetLayerIterator nodeIt = layerMinus1.begin();
    LevelSetLayerIterator nodeEnd = layerMinus1.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;

      m_OutputLevelSet->GetLayer( -1 ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( nodeIt->first, nodeIt->second ) );
      m_InternalImage->SetPixel( currentIndex, -1 );
      LevelSetLayerIterator tempIt = nodeIt;
      ++nodeIt;
      layerMinus1.erase( tempIt );

      neighIt.SetLocation( currentIndex );

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        LevelSetInputType tempIndex =
            neighIt.GetIndex( it.GetNeighborhoodOffset() );

        LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
        if( phiIt != m_TempPhi.end() )
          {
          if( phiIt->second == -3. )
            {
            phiIt->second = currentValue - 1;
            m_TempLevelSet->GetLayer( -2 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, currentValue - 1 ) );
            }
          }
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void MovePointFromPlus1()
    {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType neighOffset;
    neighOffset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      neighOffset[dim] = -1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 1;
      neighIt.ActivateOffset( neighOffset );
      neighOffset[dim] = 0;
      }

    LevelSetLayerType layerPlus1 = m_TempLevelSet->GetLayer( static_cast< char >( 1 ) );

    LevelSetLayerIterator nodeIt = layerPlus1.begin();
    LevelSetLayerIterator nodeEnd = layerPlus1.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;

      m_OutputLevelSet->GetLayer( 1 ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( nodeIt->first, nodeIt->second ) );
      m_InternalImage->SetPixel( currentIndex, 1 );
      LevelSetLayerIterator tempIt = nodeIt;
      ++nodeIt;
      layerPlus1.erase( tempIt );

      neighIt.SetLocation( currentIndex );

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        LevelSetInputType tempIndex =
            neighIt.GetIndex( it.GetNeighborhoodOffset() );

        LevelSetLayerIterator phiIt = m_TempPhi.find( tempIndex );
        if( phiIt != m_TempPhi.end() )
          {
          if( phiIt->second == 3. )
            {
            phiIt->second = currentValue + 1;
            m_TempLevelSet->GetLayer( 2 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, currentValue + 1 ) );
            }
          }
//        else
//          {
//          char status = m_InternalImage->GetPixel( tempIndex );
//          if( status == 3 )
//            {
//            m_TempPhi.insert(
//                  std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, currentValue + 1 ) );
//            m_TempLevelSet->GetLayer( 2 ).insert(
//                  std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, currentValue + 1 ) );
//            }
//          }
        }
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void MovePointFromMinus2()
    {
    LevelSetLayerType layerMinus2 = m_TempLevelSet->GetLayer( static_cast< char >( -2 ) );

    LevelSetLayerIterator nodeIt = layerMinus2.begin();
    LevelSetLayerIterator nodeEnd = layerMinus2.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType   currentIndex = nodeIt->first;

      m_OutputLevelSet->GetLayer( -2 ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( nodeIt->first, nodeIt->second ) );
      m_InternalImage->SetPixel( currentIndex, -2 );
      LevelSetLayerIterator tempIt = nodeIt;
      ++nodeIt;
      layerMinus2.erase( tempIt );
      }
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void MovePointFromPlus2()
    {
    LevelSetLayerType layerPlus2 = m_TempLevelSet->GetLayer( static_cast< char >( 2 ) );

    LevelSetLayerIterator nodeIt = layerPlus2.begin();
    LevelSetLayerIterator nodeEnd = layerPlus2.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;

      m_OutputLevelSet->GetLayer( 2 ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( nodeIt->first, nodeIt->second ) );
      m_InternalImage->SetPixel( currentIndex, 2 );
      LevelSetLayerIterator tempIt = nodeIt;
      ++nodeIt;
      layerPlus2.erase( tempIt );
      }
    }

private:
  UpdateWhitakerSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateWhitakerSparseLevelSet_h
