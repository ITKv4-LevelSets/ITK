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

#ifndef __itkUpdateSparseLevelSet_h
#define __itkUpdateSparseLevelSet_h

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
template< unsigned int VDimension, typename TLevelSetValueType >
class UpdateSparseLevelSet : public Object
{
public:
  typedef UpdateSparseLevelSet       Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( UpdateSparseLevelSet, Object );

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef TLevelSetValueType  LevelSetOutputType;

  typedef WhitakerSparseLevelSetBase< LevelSetOutputType, ImageDimension >
                                                       LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::InputType             LevelSetInputType;

  typedef std::list< LevelSetOutputType >              UpdateListType;


  typedef typename LevelSetType::SparseImageType       SparseImageType;
  typedef typename SparseImageType::Pointer            SparseImagePointer;
  typedef typename SparseImageType::IndexType          SparseImageIndexType;

  typedef typename LevelSetType::NodeAttributeType     LevelSetNodeAttributeType;
  typedef typename LevelSetType::NodeStatusType        LevelSetNodeStatusType;
  typedef typename LevelSetType::NodePairType          LevelSetNodePairType;
  typedef typename LevelSetType::NodeListType          LevelSetNodeListType;
  typedef typename LevelSetType::NodeListIterator      LevelSetNodeListIterator;
  typedef typename LevelSetType::NodeListConstIterator LevelSetNodeListConstIterator;

  typedef typename LevelSetType::SparseLayerMapType           SparseLayerMapType;
  typedef typename LevelSetType::SparseLayerMapIterator       SparseLayerMapIterator;
  typedef typename LevelSetType::SparseLayerMapConstIterator  SparseLayerMapConstIterator;

  typedef ImageRegionIteratorWithIndex< SparseImageType > SparseIteratorType;
  typedef ShapedNeighborhoodIterator< SparseImageType >   SparseNeighborhoodIteratorType;

  void DebugPrint()
    {
    std::cout << "Sn2.size() " << m_StatusLists->GetListNode( -2 )->size() <<std::endl;
    std::cout << "Sn1.size() " << m_StatusLists->GetListNode( -1 )->size() <<std::endl;
    std::cout << "Sz.size() " << m_StatusLists->GetListNode( 0 )->size() <<std::endl;
    std::cout << "Sp1.size() " << m_StatusLists->GetListNode( 1 )->size() <<std::endl;
    std::cout << "Sp2.size() " << m_StatusLists->GetListNode( 2 )->size() <<std::endl;

    std::cout << "Ln2.size() " << m_SparseLevelSet->GetListNode( -2 )->size() <<std::endl;
    std::cout << "Ln1.size() " << m_SparseLevelSet->GetListNode( -1 )->size() <<std::endl;
    std::cout << "Lz.size() " << m_SparseLevelSet->GetListNode( 0 )->size() <<std::endl;
    std::cout << "Lp1.size() " << m_SparseLevelSet->GetListNode( 1 )->size() <<std::endl;
    std::cout << "Lp2.size() " << m_SparseLevelSet->GetListNode( 2 )->size() <<std::endl;
    }

  // this is the same as Procedure 2
  // Input is a update image point m_UpdateImage
  // Input is also WhitakerSparseLevelSetBasePointer
  void UpdateZeroLevelSet()
  {
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();

    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_0 = m_SparseLevelSet->GetListNode( 0 );
    LevelSetNodePairType p;
    LevelSetOutputType update;
    LevelSetNodeAttributeType q;

    m_StatusLists->GetListNode( 1 )->clear();
    m_StatusLists->GetListNode( -1 )->clear();

    // for each point in Lz
    while( !m_Update->empty() )
      {
      p = list_0->front();

      // update the level set
      update = m_Update->front();
      p.second.m_Value += m_Dt * update;

      // if(phi(p)> .5), remove p from Lz, add p to Sp1
      if( p.second.m_Value > static_cast<LevelSetOutputType>( 0.5 ) )
        {
//        p.second.m_Status = 1;
        m_StatusLists->GetListNode( 1 )->push_back( p );
        m_SparseImage->SetPixel( p.first, p.second );
        }
      else
        {
        // if(phi(p)<-.5), remove p from Lz, add p to Sn1
        if( p.second.m_Value < static_cast<LevelSetOutputType>( -0.5 ) )
          {
//          p.second.m_Status = -1;
          m_StatusLists->GetListNode( -1 )->push_back( p );
          m_SparseImage->SetPixel( p.first, p.second );
          }
        // else keep it in Lz
        else
          {
          new_list_0.push_back( p );
          m_SparseImage->SetPixel( p.first, p.second );
          }
        }
      list_0->pop_front();
      m_Update->pop_front();
      }

    // update list_0
    while( !new_list_0.empty() )
      {
      list_0->push_back( new_list_0.front() );
      new_list_0.pop_front();
      }
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();
  }

  void UpdateMinusLevelSet( const LevelSetNodeStatusType& status )
  {
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();

    LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) + 0.5;
    LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) - 0.5;

    SparseImageIndexType idx;
    LevelSetNodePairType p;
    LevelSetNodeAttributeType q;
    ZeroFluxNeumannBoundaryCondition< SparseImageType > sp_nbc;

    typename SparseNeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    SparseNeighborhoodIteratorType
      sparseNeighborhoodIt( radius, m_SparseImage,
                            m_SparseImage->GetLargestPossibleRegion() );

    sparseNeighborhoodIt.OverrideBoundaryCondition( &sp_nbc );

    typename SparseNeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      sparseNeighborhoodIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      sparseNeighborhoodIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

    LevelSetNodeStatusType status_plus_1 = static_cast<LevelSetNodeStatusType>( status + 1 );
    LevelSetNodeStatusType status_minus_1 = static_cast<LevelSetNodeStatusType>( status - 1 );

    // for each point p in Ln1 -- status = -1
    LevelSetNodeListType* list = m_SparseLevelSet->GetListNode( status );
    LevelSetNodeListType temp_list;

    while( !list->empty() )
      {
      p = list->front();
      idx = p.first;
      sparseNeighborhoodIt.SetLocation( idx );

      LevelSetNodeStatusType M =
          NumericTraits<LevelSetNodeStatusType>::NonpositiveMin();

      // flag = is there at least one neighbor q s.t. ( q.m_Status == status + 1 )
      bool IsThereNeighborEqualToStatusPlus1 = false;

      for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        q = i.Get();
        if ( q.m_Status == status_plus_1 )
          {
          IsThereNeighborEqualToStatusPlus1 = true;
          }
        if ( ( q.m_Value > M ) && ( q.m_Status >= status_plus_1 ) )
          {
          M = q.m_Value;
          }
        }

      std::cout << "***" <<std::endl;
      std::cout << "M=" << static_cast< int >( M )<< std::endl;
      std::cout <<"IsThereNeighborEqualToStatusPlus1=" <<static_cast< int >( IsThereNeighborEqualToStatusPlus1 ) <<std::endl;
      if ( !IsThereNeighborEqualToStatusPlus1 )
        {
        // let's make sure the layer at status_minus_1 is in the active layers
        // before pushing back the current node
        if ( status_minus_1 > m_MinStatus )
          {
          std::cout << "status_minus_1 " <<static_cast< int >( status_minus_1 ) <<std::endl;
//          p.second.m_Status = status_minus_1;
//          m_SparseImage->SetPixel( p.first, p.second );
          m_StatusLists->GetListNode( status_minus_1 )->push_back( p );
          }
        else
          {
          p.second.m_Status = m_MinStatus;
          p.second.m_Value = static_cast< LevelSetOutputType >( m_MinStatus );
          m_SparseImage->SetPixel( idx, p.second );
          }
        }
      else
        {
        p.second.m_Value = static_cast< LevelSetOutputType >( M-1 );

        std::cout << p.second.m_Value << std::endl;
        if ( p.second.m_Value >= o1 )
          {
          std::cout << "entering" <<std::endl;
//          r.m_Status = status_plus_1;
          m_SparseImage->SetPixel( idx, p.second );
          m_StatusLists->GetListNode( status_plus_1 )->push_back(
                LevelSetNodePairType( idx, p.second ) );
          }
        else
          {
          if ( p.second.m_Value < o2 )
            {
            if ( status_minus_1 > m_MinStatus )
              {
//              r.m_Status = status_minus_1;
              m_SparseImage->SetPixel( idx, p.second );
              m_StatusLists->GetListNode( status_minus_1 )->push_back(
                  LevelSetNodePairType( idx, p.second ) );
              }
            else
              {
              p.second.m_Value = static_cast< LevelSetOutputType >( m_MinStatus );
              m_SparseImage->SetPixel( idx, p.second );
              }
            }
          else
            {
            temp_list.push_back( p );
            }
          }
        }
       list->pop_front();
       std::cout << "list->size() : " << list->size() <<std::endl;
      }

    while( !temp_list.empty() )
      {
      list->push_back( temp_list.front() );
      temp_list.pop_front();
      }
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();
    }


  void UpdatePlusLevelSet( const LevelSetNodeStatusType& status )
  {
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();

    LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) - 0.5;
    LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) + 0.5;

    SparseImageIndexType idx;
    LevelSetNodePairType p;
    LevelSetNodeAttributeType q;
    ZeroFluxNeumannBoundaryCondition< SparseImageType > sp_nbc;

    typename SparseNeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    SparseNeighborhoodIteratorType
    sparseNeighborhoodIt( radius, m_SparseImage,
                          m_SparseImage->GetLargestPossibleRegion() );
                          sparseNeighborhoodIt.OverrideBoundaryCondition( &sp_nbc );

    typename SparseNeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      sparseNeighborhoodIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      sparseNeighborhoodIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

    LevelSetNodeStatusType status_minus_1 = status - 1;
    LevelSetNodeStatusType status_plus_1 = status + 1;

    LevelSetNodeListType* list = m_SparseLevelSet->GetListNode( status );
    LevelSetNodeListType temp_list;

    while( !list->empty() )
      {
      p = list->front();
      idx = p.first;
      sparseNeighborhoodIt.SetLocation( idx );

      LevelSetNodeStatusType M = NumericTraits<LevelSetNodeStatusType>::max();

      bool IsThereNeighborEqualToStatusMinus1 = false;
      for( typename SparseNeighborhoodIteratorType::Iterator
            i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        q = i.Get();
        if ( q.m_Status == status_minus_1 )
          {
          IsThereNeighborEqualToStatusMinus1 = true;
          }
        if ( ( M > q.m_Value ) &&
            ( q.m_Status <= status_minus_1 ) )
          {
          M = q.m_Value;
          }
        }

      if ( !IsThereNeighborEqualToStatusMinus1 )
        {
        // let's make sure the layer at status_plus_1 is in the active layers
        // before pushing back the current node
        if ( status_plus_1 < m_MaxStatus )
          {
//          p.second.m_Status = status_plus_1;
          m_StatusLists->GetListNode( status_plus_1 )->push_back( p );
          m_SparseImage->SetPixel( p.first, p.second );
          }
        else
          {
          p.second.m_Status = m_MaxStatus;
          p.second.m_Value = static_cast< LevelSetOutputType >( m_MaxStatus );
          m_SparseImage->SetPixel( idx, p.second );
          }
        }
      else
        {
        p.second.m_Value = static_cast< LevelSetOutputType >( M+1 );

        if ( p.second.m_Value <= o1 )
          {
//          r.m_Status = status_minus_1;
          m_SparseImage->SetPixel( idx, p.second );
          m_StatusLists->GetListNode( status_minus_1 )->push_back(
                LevelSetNodePairType( idx, p.second ) );
          }
        else
          {
          if ( p.second.m_Value > o2 )
            {
            if ( status_plus_1 < m_MaxStatus )
              {
//              r.m_Status = status_plus_1;
              m_SparseImage->SetPixel( idx, p.second );
              m_StatusLists->GetListNode( status_plus_1 )->push_back(
                    LevelSetNodePairType( idx, p.second ) );
              }
            else
              {
              p.second.m_Status = m_MaxStatus;
              p.second.m_Value = static_cast< LevelSetOutputType >( m_MaxStatus );
              m_SparseImage->SetPixel( idx, p.second );
              }
            }
          else
            {
            temp_list.push_back( p );
            }
          }
        }
      list->pop_front();
      }

    while( !temp_list.empty() )
      {
      list->push_back( temp_list.front() );
      temp_list.pop_front();
      }
    std::cout <<__FILE__ <<" " <<__LINE__ <<std::endl;
    DebugPrint();
  }

  void Update()
  {
    if( m_SparseLevelSet.IsNull() )
      {
      itkGenericExceptionMacro( <<"m_SparseLevelSet is NULL" );
      }
    if( !m_Update )
      {
      itkGenericExceptionMacro( <<"m_Update is NULL" );
      }
    m_SparseImage = m_SparseLevelSet->GetImage();

    UpdateZeroLevelSet();
    UpdateMinusLevelSet( -1 );
    UpdatePlusLevelSet( 1 );
    UpdateMinusLevelSet( -2 );
    UpdatePlusLevelSet( 2 );

    std::cout << ":::::::::::" <<std::endl;
    DebugPrint();

    UpdatePointsChangingStatus();

    std::cout <<"$$$$$$$$$$$$$$$" << std::endl;
    DebugPrint();
  }

  void UpdatePointsChangingStatus()
  {
    // Move points into the zero levelset
    LevelSetNodeListType* list_0 = m_StatusLists->GetListNode( 0 );
    LevelSetNodePairType p;

    //for each point p in Sz
    while( !list_0->empty() )
      {
      p = list_0->front();

      //label(p) = 0
      p.second.m_Status = 0;
      //std::cout << p.second.m_Value << std::endl;
      // p.second.m_Value = 0;

      // add p to Lz
      m_SparseLevelSet->GetListNode( 0 )->push_back( p );

      // remove p from Sz
      list_0->pop_front();
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

  void UpdatePointsChangingStatus( const LevelSetNodeStatusType& iStatus )
  {
  int iSign = ( iStatus > 0 ) ? 1 : -1;

    // Move points into -1 and +1 level sets
  // and ensure -2, +2 neighbors
  LevelSetNodeListType* list_minus_1 = m_StatusLists->GetListNode( iStatus );
  LevelSetNodePairType p, q;

  SparseImageIndexType idx;

  typename SparseNeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  ZeroFluxNeumannBoundaryCondition< SparseImageType > sp_nbc;

  while( !list_minus_1->empty() )
    {
    p = list_minus_1->front();

    // label(p) = -1
    p.second.m_Status = iStatus;

    // add p to Ln1
    m_SparseLevelSet->GetListNode( iStatus )->push_back( p );

    // remove p from Sn1
    list_minus_1->pop_front();

    if( m_MaxStatus - static_cast< char >( vnl_math_abs( iStatus ) ) > 1 )
      {
      SparseNeighborhoodIteratorType
      sparseNeighborhoodIt( radius, m_SparseImage,
                            m_SparseImage->GetLargestPossibleRegion() );
      sparseNeighborhoodIt.OverrideBoundaryCondition( &sp_nbc );

      typename SparseNeighborhoodIteratorType::OffsetType sparse_offset;
      sparse_offset.Fill( 0 );

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
        sparse_offset[dim] = -1;
        sparseNeighborhoodIt.ActivateOffset( sparse_offset );
        sparse_offset[dim] = 1;
        sparseNeighborhoodIt.ActivateOffset( sparse_offset );
        sparse_offset[dim] = 0;
        }

      idx = p.first;
      sparseNeighborhoodIt.SetLocation( idx );

      // for each point q in N(p)
      for( typename SparseNeighborhoodIteratorType::Iterator
        i = sparseNeighborhoodIt.Begin(); !i.IsAtEnd(); ++i )
        {
        q.first = sparseNeighborhoodIt.GetIndex( i.GetNeighborhoodOffset() );
        q.second = i.Get();

        // if(phi(q)==-3)
        std::cout << q.second.m_Value <<std::endl;

        if( q.second.m_Value == static_cast< LevelSetOutputType >( iStatus + 2 * iSign ) )
          {
          // phi(q)=phi(p) + iSign
          // TODO: Check if p.second.m_Value is same as m_SparseImage->GetPixel( idx )
          q.second.m_Value = p.second.m_Value + iSign;

          // add q to Sn2
          m_SparseImage->SetPixel( q.first, q.second );
          m_StatusLists->GetListNode( iStatus + iSign )->push_back( q );
          }
        }
      }
    }
  std::cout <<"***" <<std::endl;
  DebugPrint();
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  void SetUpdate( UpdateListType* iUpdate )
    {
    m_Update = iUpdate;
    }

protected:
  UpdateSparseLevelSet() : m_Dt( NumericTraits< LevelSetOutputType >::One ),
    m_Update( NULL ), m_MinStatus( -3 ), m_MaxStatus( 3 )
    {
    m_StatusLists = LevelSetType::New();
    }
  ~UpdateSparseLevelSet() {}

  LevelSetOutputType m_Dt;

  UpdateListType*    m_Update;
  LevelSetPointer    m_SparseLevelSet;
  SparseImagePointer m_SparseImage;

  LevelSetPointer        m_StatusLists;
  LevelSetNodeStatusType m_MinStatus;
  LevelSetNodeStatusType m_MaxStatus;

private:
  UpdateSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateSparseLevelSet_h
