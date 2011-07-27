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
#include <list>
#include "itkObject.h"

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

  typedef typename LevelSetType::LabelObjectType       LevelSetLabelObjectType;
  typedef typename LevelSetType::LabelObjectPointer    LevelSetLabelObjectPointer;
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

  // this is the same as Procedure 2
  // Input is a update image point m_UpdateImage
  // Input is also WhitakerSparseLevelSetBasePointer
  void UpdateZeroLevelSet()
  {
    LevelSetLayerType& layer0 = m_OutputLevelSet->GetLayer( 0 );

    m_TempLevelSet->GetLayer( 1 ).clear();
    m_TempLevelSet->GetLayer( -1 ).clear();

    // for each point in Lz
    LevelSetLayerIterator upIt = m_Update.begin();
    LevelSetLayerIterator zeroSetIt = layer0.begin();

    while( upIt != m_Update.end() )
      {
      LevelSetInputType   currentIndex = zeroSetIt->first;
      LevelSetOutputType  currentValue = zeroSetIt->second;

      LevelSetInputType updateIndex = upIt->first;
      assert( updateIndex == currentIndex );

      // update the level set
      LevelSetOutputType tempUpdate =
          m_Dt * static_cast< LevelSetOutputType >( upIt->second );

      if ( tempUpdate > 1 )//0.5 )
      {
        tempUpdate = 1.;//0.5;
      }
      if ( tempUpdate < -1 )//0.5 )
      {
        tempUpdate = -1; //0.5;
      }

      LevelSetOutputType tempValue = currentValue + tempUpdate;
      m_RMSChangeAccumulator += tempUpdate*tempUpdate;

      // if(phi(p)> 0.5)
      if( tempValue > static_cast<LevelSetOutputType>( 0.5 ) )
        {
        // is there any point around moving in the opposite direction?
        bool ok = true;

        LevelSetInputType tempIndex = currentIndex;

        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
          for( int kk = -1; kk < 2; kk +=2 )
          {
//            if( ( tempIndex[dim] > 0 ) && ( tempIndex[dim] < imageSize[dim] - 1 ) )
            {
              tempIndex[dim] = currentIndex[dim] + kk;
              char tempStatus = m_OutputLevelSet->Status( tempIndex);

              if( tempStatus == 0 )
              {
                if( m_TempPhi[ tempIndex ] < -0.5 )
                {
                  ok = false;
                }
              }
            }
            tempIndex[dim] = currentIndex[dim];
          }
          tempIndex[dim] = currentIndex[dim];
        }

        // no point move in the opposite direction
        if( ok )
          {
          LevelSetLayerIterator tempIt = zeroSetIt;
          ++zeroSetIt;
          // remove p from Lz
          layer0.erase( tempIt );

          // update m_TempPhi
          currentValue = tempValue;
          m_TempPhi[ currentIndex ] = currentValue;

          // add p to Sp1
          m_TempLevelSet->GetLayer( 1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
        else // at least one point is moving in the opposite direction
          {
          // keep it in Lz
          ++zeroSetIt;
          }
        }
      else
        {
        // if(phi(p)<-0.5)
        if( tempValue < static_cast<LevelSetOutputType>( -0.5 ) )
          {
          // is there any point around moving in the opposite direction?
          bool ok = true;
          LevelSetInputType tempIndex = currentIndex;

          for( unsigned int dim = 0; dim < ImageDimension; dim++ )
          {
            for( int kk = -1; kk < 2; kk += 2 )
            {
//              if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
              {
                tempIndex[dim] = currentIndex[dim] + kk;
                char tempStatus = m_OutputLevelSet->Status( tempIndex );

                if( tempStatus == 0 )
                {
                  if( m_TempPhi[ tempIndex ] > 0.5 )
                  {
                    ok = false;
                  }
                }
              }
              tempIndex[dim] = currentIndex[dim];
            }
          }
          if( ok )
            {
            LevelSetLayerIterator tempIt = zeroSetIt;
            ++zeroSetIt;
            // remove p from Lz
            layer0.erase( tempIt );

            currentValue = tempValue;

            // update m_TempPhi
            m_TempPhi[ currentIndex ] = currentValue;

            // add p to Sn1
            m_TempLevelSet->GetLayer( -1 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
            }
          else // at least one point is moving in the opposite direction
            {
            // keep it in Lz
            ++zeroSetIt;
            }
          }
        else // else keep it in Lz but update m_TempPhi
          {
          zeroSetIt->second = tempValue;
          m_TempPhi[ currentIndex ] = tempValue;
          ++zeroSetIt;
          }
        }
      ++upIt;
      }
    }

  void UpdateMinusLevelSet( const char& status )
  {
    const LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) + 0.5;
    const LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) - 0.5;

    const char status_plus_1 = status + 1;
    const char status_minus_1 = status - 1;

    // for each point p in L_{status}
    LevelSetLayerType& list = m_OutputLevelSet->GetLayer( status );
    LevelSetLayerIterator layerIt = list.begin();

    while( layerIt != list.end() )
      {
      LevelSetInputType   currentIndex = layerIt->first;
      LevelSetOutputType  currentValue = layerIt->second;

      LevelSetOutputType M =
        NumericTraits<LevelSetOutputType>::NonpositiveMin();

      // flag = is there at least one neighbor q s.t. ( q.m_Status == status + 1 ) ?
      bool IsThereNeighborEqualToStatusPlus1 = false;

      LevelSetInputType tempIndex = currentIndex;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk += 2 )
        {
//          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            tempIndex[dim] = currentIndex[dim] + kk;

            const char tempStatus = m_OutputLevelSet->Status( tempIndex );

            if( tempStatus == status_plus_1 )
            {
              IsThereNeighborEqualToStatusPlus1 = true;
            }
            if( tempStatus >= status_plus_1 )
            {
              LevelSetOutputType tempValue = m_TempPhi[ tempIndex ];

              if ( tempValue > M )
              {
                // M = \max_{q \in N(p)} m_TempPhi[ q ]
                M = tempValue;
              }
            }
            tempIndex[dim] = currentIndex[dim];
          }
        }
        tempIndex[dim] = currentIndex[dim];
      }


      if ( !IsThereNeighborEqualToStatusPlus1 )
        {
          // let's make sure the layer at status_minus_1 is in the active layers
          // before pushing back the current node
          if ( status_minus_1 > m_MinStatus )
          {
            LevelSetLayerIterator tempIt = layerIt;
            ++layerIt;

            // remove p from the current layer L_{status}
            list.erase( tempIt );

            // add it to S_{status - 1}
            m_TempLevelSet->GetLayer( status_minus_1 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
          else // it is out of the layers (inside)
          {
            // we only need to remove it from the layers.
            LevelSetLayerIterator tempIt = layerIt;
            ++layerIt;

            m_TempPhi.erase( currentIndex );
            list.erase( tempIt );

            // set label( p ) = m_MinStatus
            m_OutputLevelSet->GetLabelMap()->SetPixel( currentIndex, m_MinStatus );
          }
        }
      else
        {
        currentValue = M-1.;

        if ( currentValue >= o1 )
          {
          LevelSetLayerIterator tempIt = layerIt;
          ++layerIt;

          // remove from the current layer L_{status}
          list.erase( tempIt );

          // update m_TempPhi
          m_TempPhi[ currentIndex ] = currentValue;

          // add it to S_{status + 1}
          m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
        else
          {
          if ( currentValue < o2 )
            {
            if ( status_minus_1 > m_MinStatus )
              {
                LevelSetLayerIterator tempIt = layerIt;
                ++layerIt;

                // remove from the current layer L_{status}
                list.erase( tempIt );

                // update m_TempPhi
                m_TempPhi[ currentIndex ] = currentValue;

                // add it to S_{status - 1}
                m_TempLevelSet->GetLayer( status_minus_1 ).insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
              }
            else // it is out of the layers (inside)
            {
              // remove it from the layers
              LevelSetLayerIterator tempIt = layerIt;
              ++layerIt;

              m_TempPhi.erase( currentIndex );
              list.erase( tempIt );

              // update the corresponding label to m_MinStatus
              m_OutputLevelSet->GetLabelMap()->SetPixel( currentIndex, m_MinStatus );
            }
            }
          else // we keep in the layer
            {
            // update m_TempPhi
            m_TempPhi[currentIndex] = currentValue;
            ++layerIt;
            }
          }
        }
     }
  }


  void UpdatePlusLevelSet( const char& status )
  {
    const LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) - 0.5;
    const LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) + 0.5;

    const char status_minus_1 = status - 1;
    const char status_plus_1 = status + 1;

    // for each point p in L_{status}
    LevelSetLayerType& list = m_OutputLevelSet->GetLayer( status );
    LevelSetLayerIterator layerIt = list.begin();

    while( layerIt != list.end() )
      {
      LevelSetInputType   currentIndex = layerIt->first;
      LevelSetOutputType  currentValue = layerIt->second;

      LevelSetOutputType M = NumericTraits<LevelSetOutputType>::max();

      // flag = is there at least one neighbor q s.t. ( q.m_Status == status - 1 ) ?
      bool IsThereNeighborEqualToStatusMinus1 = false;

      LevelSetInputType tempIndex = currentIndex;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk += 2 )
        {
//          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            tempIndex[dim] = currentIndex[dim] + kk;

            const char tempStatus = m_OutputLevelSet->Status( tempIndex );

            if( tempStatus == status_minus_1 )
            {
              IsThereNeighborEqualToStatusMinus1 = true;
            }
            if ( tempStatus <= status_minus_1 )
            {
              LevelSetOutputType tempValue = m_TempPhi[ tempIndex ];

              if( tempValue < M )
              {
                // M = M = \min_{q \in N(p)} m_TempPhi[ q ]
                M = tempValue;
              }
            }
            tempIndex[dim] = currentIndex[dim];
          }
        }
        tempIndex[dim] = currentIndex[dim];
      }

      if ( !IsThereNeighborEqualToStatusMinus1 )
        {
        // let's make sure the layer at status_plus_1 is in the active layers
        // before pushing back the current node
        if ( status_plus_1 < m_MaxStatus )
        {
          LevelSetLayerIterator tempIt = layerIt;
          ++layerIt;

          // remove p from the current layer L_{status}
          list.erase( tempIt );

          // add it to S_{status + 1}
          m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
        }
        else // it is out of the layers (outside)
          {
            // we only need to remove it from the layers.
            LevelSetLayerIterator tempIt = layerIt;
            ++layerIt;

            m_TempPhi.erase( currentIndex );
            list.erase( tempIt );

            // set label( p ) = m_MaxStatus
            m_OutputLevelSet->GetLabelMap()->SetPixel( currentIndex, m_MaxStatus );
          }
        }
      else
        {
        currentValue = M+1;

        if ( currentValue <= o1 )
        {
          LevelSetLayerIterator tempIt = layerIt;
          ++layerIt;

          // remove from the current layer L_{status}
          list.erase( tempIt );

          // update m_TempPhi
          m_TempPhi[ currentIndex ] = currentValue;

          // add it to S_{status - 1}
          m_TempLevelSet->GetLayer( status_minus_1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                    currentValue ) );
        }
        else
          {
          if ( currentValue > o2 )
            {
            if ( status_plus_1 < m_MaxStatus )
            {
              LevelSetLayerIterator tempIt = layerIt;
              ++layerIt;

              // remove from the current layer L_{status}
              list.erase( tempIt );

              // update m_TempPhi
              m_TempPhi[ currentIndex ] = currentValue;

              // add it to S_{status + 1}
              m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                        currentValue ) );
            }
            else // it is out of the layers (outside)
              {
                // remove it from the layers
                LevelSetLayerIterator tempIt = layerIt;
                ++layerIt;

                m_TempPhi.erase( currentIndex );
                list.erase( tempIt );

                // update the corresponding label to m_MaxStatus
                m_OutputLevelSet->GetLabelMap()->SetPixel( currentIndex, m_MaxStatus );
              }
            }
          else // we keep in the layer
          {
            // update m_TempPhi
            m_TempPhi[currentIndex] = currentValue;
            ++layerIt;
          }
          }
        }
      }
  }

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

    m_TempPhi.clear();

    for( char status = -2; status < 3; status++ )
      {
      LevelSetLayerType layer = m_SparseLevelSet->GetLayer( status );

      LevelSetLayerConstIterator it = layer.begin();
      while( it != layer.end() )
        {
        m_TempPhi[ it->first ] = it->second;
        ++it;
        }
      }

    UpdateZeroLevelSet();
    UpdateMinusLevelSet( -1 );
    UpdatePlusLevelSet( 1 );
    UpdateMinusLevelSet( -2 );
    UpdatePlusLevelSet( 2 );

    UpdatePointsChangingStatus();

    m_TempPhi.clear();
  }

  void UpdatePointsChangingStatus()
  {
    // Move points into the zero levelset
    LevelSetLayerType list0 = m_TempLevelSet->GetLayer( 0 );

    LevelSetLayerIterator it = list0.begin();

    while( it != list0.end() )
    {
      m_OutputLevelSet->GetLabelMap()->SetPixel( it->first, 0 );
      m_OutputLevelSet->GetLayer(0).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( it->first, it->second ) );
      ++it;
    }

    UpdatePointsChangingStatus( -1 );

    //for each point in Sp1
    //label(p) = 1, add p to Lp1, remove p from Sp1
    //for each point q in N(p)
    //if(phi(q)== 3), phi(q)=phi(p)+1, add q to Sp2
    UpdatePointsChangingStatus( 1 );

    // Move points into -2 and +2 level sets
    //for each point p in Sn2
    //label(p) = -2, add p to Ln2, remove p from Sn2
    UpdatePointsChangingStatus( -2 );

    //for each point p in Sp2
    //label(p) = 2, add p to Lp2, remove p from Sp2
    UpdatePointsChangingStatus( 2 );
  }

  void UpdatePointsChangingStatus( char iStatus )
  {
    int iSign = ( iStatus > 0 ) ? 1 : -1;

    // Move points into -1 and +1 level sets
    // and ensure -2, +2 neighbors
    LevelSetLayerType& listMinus1 = m_TempLevelSet->GetLayer( iStatus );

    LevelSetLayerIterator layerIt = listMinus1.begin();

    while( layerIt != listMinus1.end() )
      {
      LevelSetInputType   currentIndex = layerIt->first;
      LevelSetOutputType  currentValue = layerIt->second;

      // label(p) = iStatus
      m_OutputLevelSet->GetLabelMap()->SetPixel( currentIndex, iStatus );

      // add p to L_{iStatus}
      m_OutputLevelSet->GetLayer( iStatus ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                currentValue ) );

      LevelSetLayerIterator tempIt = layerIt;
      ++layerIt;

      // remove p from S_{iStatus}
      listMinus1.erase( tempIt );

      if( m_MaxStatus - static_cast< char >( vnl_math_abs( iStatus ) ) > 1 )
      {
        LevelSetInputType tempIndex = currentIndex;

        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
          for( int kk = -1; kk < 2; kk += 2 )
          {
//            if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
            {
              tempIndex[dim] = currentIndex[dim] + kk;

              char tempStatus = m_OutputLevelSet->Status( tempIndex );
              LevelSetOutputType tempValue = m_TempPhi[ tempIndex ];

              if( iStatus == -1 )
              {
                assert( iStatus + 2 * iSign == -3 );
              }
              if( iStatus == 1 )
              {
                assert( iStatus + 2 * iSign == 3 );
              }

              if( tempValue == static_cast< LevelSetOutputType >( iStatus + 2 * iSign ) )
              {
                tempStatus += iSign;
                tempValue += static_cast< LevelSetOutputType >( iSign );
                m_TempLevelSet->GetLayer( tempStatus ).insert(
                      std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, tempValue ) );
              }
            }
            tempIndex[dim] = currentIndex[dim];
          }
          tempIndex[dim] = currentIndex[dim];
        }
      }
    }
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

private:
  UpdateWhitakerSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateWhitakerSparseLevelSet_h
