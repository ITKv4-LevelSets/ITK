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

      if ( tempUpdate > 0.5 )
      {
        tempUpdate = 0.5;
      }
      if ( tempUpdate < -0.5 )
      {
        tempUpdate = -0.5;
      }

      LevelSetOutputType tempValue = currentValue + tempUpdate;
      m_RMSChangeAccumulator += tempUpdate*tempUpdate;

      // if(phi(p)> .5), remove p from Lz, add p to Sp1
      if( tempValue > static_cast<LevelSetOutputType>( 0.5 ) )
        {
        bool ok = true;
        LevelSetInputType tempIndex = currentIndex;

        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
          for( int kk = -1; kk < 2; kk +=2 )
          {
//            if( ( tempIndex[dim] > 0 ) && ( tempIndex[dim] < imageSize[dim] - 1 ) )
            {
              tempIndex[dim] = currentIndex[dim] + kk;
              LevelSetLayerIterator tempIt = layer0.find( tempIndex );

              if( tempIt != layer0.end() )
              {
                if( tempIt->second < -0.5 )
                {
                  ok = false;
                }
              }
            }
          }
          tempIndex[dim] = currentIndex[dim];
        }

        if( ok )
          {
          currentValue = tempValue;
          m_TempLevelSet->GetLayer( 1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
        else
          {
          m_TempLevelSet->GetLayer( 0 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue) );
          }
        }
      else
        {
        // if(phi(p)<-.5), remove p from Lz, add p to Sn1
        if( tempValue < static_cast<LevelSetOutputType>( -0.5 ) )
          {
          bool ok = true;
          LevelSetInputType tempIndex = currentIndex;

          for( unsigned int dim = 0; dim < ImageDimension; dim++ )
          {
            for( int kk = -1; kk < 2; kk += 2 )
            {
//              if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
              {
                tempIndex[dim] = currentIndex[dim] + kk;
                LevelSetLayerIterator tempIt = layer0.find( tempIndex );

                if( tempIt != layer0.end() )
                {
                  if( tempIt->second > 0.5 )
                  {
                    ok = false;
                  }
                }
              }
            }
          }
          if( ok )
            {
            layer0.erase( zeroSetIt );

            currentValue = tempValue;
            m_TempLevelSet->GetLayer( -1 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
            }
//          else
//            {
//            m_TempLevelSet->GetLayer( 0 ).insert(
//                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue) );
//            }
          }
        // else keep it in Lz
//        else
//          {
//          m_TempLevelSet->GetLayer( 0 ).insert(
//                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue) );
//          }
        }
      ++zeroSetIt;
      ++upIt;
      }
    }

  void UpdateMinusLevelSet( const char& status )
  {
    const LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) + 0.5;
    const LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) - 0.5;

    const char status_plus_1 = status + 1;
    const char status_minus_1 = status - 1;

    // for each point p in Ln1 -- status = -1
    LevelSetLayerType& list = m_OutputLevelSet->GetLayer( status );
    LevelSetLayerIterator layerIt = list.begin();

    while( layerIt != list.end() )
      {
      LevelSetInputType   currentIndex = layerIt->first;
      LevelSetOutputType  currentValue = layerIt->second;

      LevelSetOutputType M =
        NumericTraits<LevelSetOutputType>::NonpositiveMin();

      // flag = is there at least one neighbor q s.t. ( q.m_Status == status + 1 )
      bool IsThereNeighborEqualToStatusPlus1 = false;

      LevelSetInputType tempIndex = currentIndex;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk += 2 )
        {
//          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            char tempStatus = 0;
            LevelSetOutputType tempValue = 0.;

            tempIndex[dim] = currentIndex[dim] + kk;
            m_TempLevelSet->StatusAndValue( tempIndex, tempStatus, tempValue );

            if( tempStatus == status_plus_1 )
            {
              IsThereNeighborEqualToStatusPlus1 = true;
            }
            if( ( tempValue > M ) && ( tempStatus >= status_plus_1 ) )
            {
              M = tempValue;
            }
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
            list.erase( layerIt );
            m_TempLevelSet->GetLayer( status_minus_1 ).insert(
                  std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
          // else it is out of the layers (inside)
          else
          {
            // we only need to remove it from the layers.
            list.erase( layerIt );
          }
        }
      else
        {
        currentValue = M-1;

        if ( currentValue >= o1 )
          {
          list.erase( layerIt );
          m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
          }
        else
          {
          if ( currentValue < o2 )
            {
            if ( status_minus_1 > m_MinStatus )
              {
              list.erase( layerIt );
              m_TempLevelSet->GetLayer( status_minus_1 ).insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
              }
            // else it is out of the layers (inside)
            else
            {
              // we only need to remove it from the layers
              list.erase( layerIt );
            }
            }
          else
            {
            // we keep in the layer
            }
          }
        }
     ++layerIt;
     }

  }


  void UpdatePlusLevelSet( const char& status )
  {
//     LevelSetOutputRealType oldValue, newValue;
    const LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) - 0.5;
    const LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) + 0.5;

    const char status_minus_1 = status - 1;
    const char status_plus_1 = status + 1;

    LevelSetLayerType& list = m_OutputLevelSet->GetLayer( status );
    LevelSetLayerIterator layerIt = list.begin();

    while( layerIt != list.end() )
      {
      LevelSetInputType   currentIndex = layerIt->first;
      LevelSetOutputType  currentValue = layerIt->second;

      LevelSetOutputType M = NumericTraits<LevelSetOutputType>::max();

      bool IsThereNeighborEqualToStatusMinus1 = false;

      LevelSetInputType tempIndex = currentIndex;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk += 2 )
        {
//          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            char tempStatus = 0;
            LevelSetOutputType tempValue = 0.;

            tempIndex[dim] = currentIndex[dim] + kk;
            m_TempLevelSet->StatusAndValue( tempIndex, tempStatus, tempValue );

            if( tempStatus == status_minus_1 )
            {
              IsThereNeighborEqualToStatusMinus1 = true;
            }
            if ( ( tempValue < M ) &&
                ( tempStatus <= status_minus_1 ) )
            {
              M = tempValue;
            }
          }
        }
      }

      if ( !IsThereNeighborEqualToStatusMinus1 )
        {
        // let's make sure the layer at status_plus_1 is in the active layers
        // before pushing back the current node
        if ( status_plus_1 < m_MaxStatus )
        {
          list.erase( layerIt );
          m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, currentValue ) );
        }
        else
          {
          // then it is outside
          }
        }
      else
        {
        currentValue = M+1;
//        p.second.m_Value = M+1;

        if ( currentValue <= o1 )
        {
          list.erase( layerIt );
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
              list.erase( layerIt );
              m_TempLevelSet->GetLayer( status_plus_1 ).insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                        currentValue ) );
            }
            else
              {
              // then it is outside
              }
            }
          }
        }
      ++layerIt;
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

    m_OutputLevelSet->SetLabelObject( m_SparseLevelSet->GetLabelObject() );

    UpdateZeroLevelSet();
    UpdateMinusLevelSet( -1 );
    UpdatePlusLevelSet( 1 );
    UpdateMinusLevelSet( -2 );
    UpdatePlusLevelSet( 2 );

    UpdatePointsChangingStatus();
  }

  void UpdatePointsChangingStatus()
  {
    // Move points into the zero levelset
    LevelSetLayerType list0 = m_TempLevelSet->GetLayer( 0 );

    LevelSetLayerIterator it = list0.begin();

    while( it != list0.end() )
    {
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

      // add p to L_{iStatus}
      m_OutputLevelSet->GetLayer( iStatus ).insert(
            std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex,
                                                                currentValue ) );

      // remove p from S_{iStatus}
      listMinus1.erase( layerIt );

      if( m_MaxStatus - static_cast< char >( vnl_math_abs( iStatus ) ) > 1 )
      {
        LevelSetInputType tempIndex = currentIndex;

        for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
          for( int kk = -1; kk < 2; kk += 2 )
          {
//            if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
            {
              char tempStatus = 0;
              LevelSetOutputType tempValue = 0.;

              tempIndex[dim] = currentIndex[dim] + kk;

              m_OutputLevelSet->StatusAndValue( tempIndex, tempStatus, tempValue );

              if( tempValue == static_cast< LevelSetOutputType >( iStatus + 2 * iSign ) )
              {
                tempStatus += iSign;
                tempValue += static_cast< LevelSetOutputType >( iSign );
                m_TempLevelSet->GetLayer( tempStatus ).insert(
                      std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, tempValue ) );
              }
            }
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

  LevelSetPointer   m_TempLevelSet;
  char              m_MinStatus;
  char              m_MaxStatus;

private:
  UpdateWhitakerSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateWhitakerSparseLevelSet_h
