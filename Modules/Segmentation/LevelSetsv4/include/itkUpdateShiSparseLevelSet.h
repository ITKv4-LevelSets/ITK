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

#ifndef __itkUpdateShiSparseLevelSet_h
#define __itkUpdateShiSparseLevelSet_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkShiSparseLevelSetBase.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <list>
#include "itkObject.h"

namespace itk
{
template< unsigned int VDimension >
class UpdateShiSparseLevelSet : public Object
{
public:
  typedef UpdateShiSparseLevelSet       Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( UpdateShiSparseLevelSet, Object );

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef ShiSparseLevelSetBase< ImageDimension >      LevelSetType;
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

  // this is the same as Procedure 2
  // Input is a update image point m_UpdateImage
  // Input is also ShiSparseLevelSetBasePointer
  void UpdateL_out()
  {
    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_out = m_SparseLevelSet->GetListNode( 1 );
    LevelSetNodePairType p;
    LevelSetOutputType update;
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

    m_StatusLists->GetListNode( 1 )->clear();
    m_StatusLists->GetListNode( -1 )->clear();

    // for each point in Lz
    while( !m_Update[1]->empty() )
      {
      p = list_out->front();

      list_out->pop_front();

      // update the level set
      update = m_Update[1]->front();

      if( update > NumericTraits< LevelSetOutputType >::Zero )
        {
        if( Con( p.first, p.second.status, update ) )
          {
          // CheckIn
          p.second.m_Status = -1;
          p.second.m_Value = -1.;
          m_StatusLists->GetListNode( -1 )->push_back( p );

          sparseNeighborhoodIt.SetLocation( p.first );

          for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            q = i.Get();
            if ( q.m_Status == 3 )
              {
              q.second.m_Status = 1;
              q.second.m_Value = 1;
              m_StatusLists->GetListNode( 1 )->push_back( q );
              // compute the update here of q;
              }
            }
          }
        else
          {
          new_list_out.push_back( p );
          }
        }
      else
        {
        new_list_out.push_back( p );
        }
      }

    while( !new_list_out.empty() )
      {
      list_out->push_back( new_list_out.front() );
      new_list_out.pop_front();
      }
    }

  void UpdateL_in()
  {
    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_in = m_SparseLevelSet->GetListNode( -1 );
    LevelSetNodePairType p;
    LevelSetOutputType update;
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

    m_StatusLists->GetListNode( 1 )->clear();
    m_StatusLists->GetListNode( -1 )->clear();

    // for each point in Lz
    while( !m_Update[-1]->empty() )
      {
      p = list_in->front();

      list_in->pop_front();

      // update the level set
      update = m_Update[-1]->front();

      if( update < NumericTraits< LevelSetOutputType >::Zero )
        {
        if( Con( p.first, p.second.status, update ) )
          {
          // CheckOut
          p.second.m_Status = 1;
          p.second.m_Value = 1.;
          m_StatusLists->GetListNode( 1 )->push_back( p );

          sparseNeighborhoodIt.SetLocation( p.first );

          for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            q = i.Get();
            if ( q.m_Status == -3 )
              {
              q.second.m_Status = -1;
              q.second.m_Value = -1;
              m_StatusLists->GetListNode( -1 )->push_back( q );
              // compute the update here of q;
              }
            }
          }
        else
          {
          new_list_in.push_back( p );
          }
        }
      else
        {
        new_list_in.push_back( p );
        }
      }

    while( !new_list_in.empty() )
      {
      list_out->push_back( new_list_in.front() );
      new_list_in.pop_front();
      }
    }

  bool Con( const SparseImageIndexType& iIdx,
            const LevelSetNodeStatusType& iCurrentStatus,
            const LevelSetOutputType& iCurrentUpdate )
  {
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

    sparseNeighborhoodIt.SetLocation( iIdx );

    opposite_status = ( CurrentStatus == 1 ) ? -1 : 1;

    for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
      {
      q = i.Get();
      if ( q.m_Status == opposite_status )
        {
        if ( q.m_Value * CurrentUpdate < 0. )
          {
          return true;
          }
        }
      }
   return false;
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

    UpdatePointsChangingStatus();
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  void SetUpdate( UpdateListType* iUpdate )
    {
    m_Update = iUpdate;
    }

protected:
  UpdateShiSparseLevelSet() : m_Dt( NumericTraits< LevelSetOutputType >::One ),
    m_Update( NULL ), m_MinStatus( -3 ), m_MaxStatus( 3 )
    {
    m_StatusLists = LevelSetType::New();
    }
  ~UpdateShiSparseLevelSet() {}

  LevelSetOutputType m_Dt;

  UpdateListType*    m_Update;
  LevelSetPointer    m_SparseLevelSet;
  SparseImagePointer m_SparseImage;

  LevelSetPointer        m_StatusLists;
  LevelSetNodeStatusType m_MinStatus;
  LevelSetNodeStatusType m_MaxStatus;

private:
  UpdateShiSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateShiSparseLevelSet_h
