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
  typedef typename LevelSetType::OutputType            LevelSetOutputType;
  typedef typename LevelSetType::OutputRealType        LevelSetOutputRealType;

  typedef std::list< LevelSetOutputRealType >          UpdateListType;

  typedef typename LevelSetType::SparseImageType       SparseImageType;
  typedef typename SparseImageType::Pointer            SparseImagePointer;
  typedef typename SparseImageType::IndexType          SparseImageIndexType;

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
    LevelSetNodeListType* list_out = m_SparseLevelSet->GetListNode( 1 );
    LevelSetNodeListType* list_in = m_SparseLevelSet->GetListNode( -1 );
    LevelSetNodePairType p;
    LevelSetOutputRealType update;
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

    if( m_Update[1]->size() != list_out->size() )
      {
      itkGenericExceptionMacro( "m_Update[1]->size() != list_out->size()" );
      }

    // for each point in Lz
    while( !m_Update[1]->empty() )
      {
      p = list_out->front();
      list_out->pop_front();

      // update the level set
      update = m_Update[1]->front();
      m_Update[1]->pop_front();

      if( update > NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( Con( p.first, p.second , update ) )
          {
          // CheckIn
          p.second = -1;
          this->m_SparseImage->SetPixel( p.first, p.second );
          m_StatusLists->GetListNode( -1 )->push_back( p );

          sparseNeighborhoodIt.SetLocation( p.first );

          for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            LevelSetOutputType q = i.Get();
            if ( q == 3 )
              {
              LevelSetNodePairType temp;
              temp.first =
                sparseNeighborhoodIt.GetIndex( i.GetNeighborhoodOffset() );
              temp.second = 1;
              m_StatusLists->GetListNode( 1 )->push_back( temp);
              this->m_SparseImage->SetPixel( temp.first, temp.second );
              // compute the update here of q;
              }
            }
          }
        else
          {
          m_StatusLists->GetListNode( 1 )->push_back( p );
          }
        }
      else
        {
        m_StatusLists->GetListNode( 1 )->push_back( p );
        }
      }

    while( !m_StatusLists->GetListNode( -1 )->empty() )
      {
      list_in->push_back( m_StatusLists->GetListNode( -1 )->front() );
      m_StatusLists->GetListNode( -1 )->pop_front();
      }

    while( !m_StatusLists->GetListNode( 1 )->empty() )
      {
      list_out->push_back( m_StatusLists->GetListNode( 1 )->front() );
      m_StatusLists->GetListNode( 1 )->pop_front();
      }
    }

  void UpdateL_in()
  {
    LevelSetNodeListType* list_in = m_SparseLevelSet->GetListNode( -1 );
    LevelSetNodeListType* list_out = m_SparseLevelSet->GetListNode( 1 );
    LevelSetNodePairType p;
    LevelSetOutputRealType update;
    LevelSetOutputType q;
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
      m_Update[-1]->pop_front();

      if( update < NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( Con( p.first, p.second , update ) )
          {
          // CheckOut
          p.second = 1;
          this->m_SparseImage->SetPixel( p.first, p.second );
          m_StatusLists->GetListNode( 1 )->push_back( p );

          sparseNeighborhoodIt.SetLocation( p.first );

          for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            q = i.Get();
            if ( q == -3 )
              {
              LevelSetNodePairType temp;
              temp.first =
                sparseNeighborhoodIt.GetIndex( i.GetNeighborhoodOffset() );
              temp.second = -1;
              m_StatusLists->GetListNode( -1 )->push_back( temp);
              this->m_SparseImage->SetPixel( temp.first, temp.second );

              // compute the update here of q;
              }
            }
          }
        else
          {
          m_StatusLists->GetListNode( -1 )->push_back( p );
          }
        }
      else
        {
        m_StatusLists->GetListNode( -1 )->push_back( p );
        }
      }

    while( !m_StatusLists->GetListNode( -1 )->empty() )
      {
      list_in->push_back( m_StatusLists->GetListNode( -1 )->front() );
      m_StatusLists->GetListNode( -1 )->pop_front();
      }

    while( !m_StatusLists->GetListNode( 1 )->empty() )
      {
      list_out->push_back( m_StatusLists->GetListNode( 1 )->front() );
      m_StatusLists->GetListNode( 1 )->pop_front();
      }
    }

  bool Con( const SparseImageIndexType& iIdx,
            const LevelSetOutputType& iCurrentStatus,
            const LevelSetOutputRealType& iCurrentUpdate ) const
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

    LevelSetOutputType opposite_status = ( iCurrentStatus == 1 ) ? -1 : 1;

    LevelSetOutputType q;

    for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
      {
      q = i.Get();
      if ( q == opposite_status )
        {
        if ( q * iCurrentUpdate < NumericTraits< LevelSetOutputType >::Zero )
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
    if( m_Update.empty() )
      {
      itkGenericExceptionMacro( <<"m_Update is empty" );
      }
    m_SparseImage = m_SparseLevelSet->GetImage();

    UpdateL_out();

    // neighborhood iterator
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

    LevelSetNodePairType p;
    LevelSetOutputType q;

    // for each point x in L_in
    while( !m_SparseLevelSet->GetListNode( -1 )->empty() )
      {
      p = m_SparseLevelSet->GetListNode( -1 )->front();
      m_SparseLevelSet->GetListNode( -1 )->pop_front();

      sparseNeighborhoodIt.SetLocation( p.first );

      bool to_be_deleted = true;

      for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        q = i.Get();
        if ( q > NumericTraits< LevelSetOutputType >::Zero )
          {
          to_be_deleted = false;
          break;
          }
        }
      if( to_be_deleted )
        {
        p.second = -3;
        this->m_SparseImage->SetPixel( p.first, p.second );
        }
      else
        {
        m_StatusLists->GetListNode( -1 )->push_back( p );
        }
      }

    while( !m_StatusLists->GetListNode( -1 )->empty() )
      {
      m_SparseLevelSet->GetListNode( -1 )->push_back(
        m_StatusLists->GetListNode( -1 )->front() );
      m_StatusLists->GetListNode( -1 )->pop_front();
      }

    UpdateL_in();

    while( !m_SparseLevelSet->GetListNode( 1 )->empty() )
      {
      p = m_SparseLevelSet->GetListNode( 1 )->front();
      m_SparseLevelSet->GetListNode( 1 )->pop_front();

      sparseNeighborhoodIt.SetLocation( p.first );

      bool to_be_deleted = true;

      for( typename SparseNeighborhoodIteratorType::Iterator
              i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        q = i.Get();
        if ( q < NumericTraits< LevelSetOutputType >::Zero )
          {
          to_be_deleted = false;
          break;
          }
        }
      if( to_be_deleted )
        {
        p.second = 3;
        this->m_SparseImage->SetPixel( p.first, p.second );
        }
      else
        {
        m_StatusLists->GetListNode( 1 )->push_back( p );
        }
      }

    while( !m_StatusLists->GetListNode( 1 )->empty() )
      {
      m_SparseLevelSet->GetListNode( 1 )->push_back(
        m_StatusLists->GetListNode( 1 )->front() );
      m_StatusLists->GetListNode( 1 )->pop_front();
      }

  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  void SetUpdate( const std::map< char, UpdateListType* >& iUpdate )
    {
    if( ( iUpdate.find( -1 ) != iUpdate.end() ) &&
        ( iUpdate.find( 1 ) != iUpdate.end() ) &&
        ( iUpdate.size() == 2 ) )
      {
      m_Update = iUpdate;
      }
    }

protected:
  UpdateShiSparseLevelSet()
    {
    m_StatusLists = LevelSetType::New();
    }
  ~UpdateShiSparseLevelSet() {}

  std::map< LevelSetOutputType, UpdateListType* > m_Update;

  LevelSetPointer    m_SparseLevelSet;
  SparseImagePointer m_SparseImage;

  LevelSetPointer     m_StatusLists;

private:
  UpdateShiSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateShiSparseLevelSet_h
