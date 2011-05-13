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
template< typename TLevelSetValueType >
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

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

  typedef TLevelSetValueType  LevelSetOutputType;

  typedef std::pair< IndexType, LevelSetOutputType > UpdateType;
  typedef std::list< UpdateType >                    UpdateListType;

  typedef SparseLevelSetBase< LevelSetOutputType, ImageDimension >
                                                       LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::InputType             LevelSetInputType;

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
  // Input is also WhitakerSparseLevelSetBasePointer
  void UpdateZeroLevelSet()
  {
    /* Procedure 2 Update Level Set Lists
    // Update the zero level set
    1: for each point p in Lz
    2: add F(p) to phi(p)
    3: if(phi(p)> .5), remove p from Lz, add p to Sp1
    4: if(phi(p)<-.5), remove p from Lz, add p to Sn1
    */

    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_0 = m_SparseLevelSet->GetListNode( 0 );
    LevelSetNodePairType p;
    LevelSetOutputType update;
    LevelSetNodeAttributeType q;

    while( !list_0->empty() && !m_Update->empty() )
    {
      p = list_0->front();
      update = m_Update->front();

      q = p.second;
      q.m_Value += m_Dt * update;
      m_SparseImage->SetPixel( p.first, p );

      if( q.m_Value > static_cast<LevelSetOutputType>( 0.5 ) )
      {
        S_plus_1.push_back( p );
      }
      else
      {
        if( q.m_Value < static_cast<LevelSetOutputType>( -0.5 ) )
        {
          S_neg_1.push_back( p );
        }
        else
        {
          new_list_0.push_back( p );
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
  }

  void UpdateMinusLevelSet( LevelSetNodeStatusType& status )
  {
    LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) + 0.5;
    LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) - 0.5;

    bool flag;
    SparseImageIndexType idx;
    LevelSetNodePairType p;
    LevelSetNodeAttributeType q, r;
    ZeroFluxNeumannBoundaryCondition< SparseImageType > sp_nbc;

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

    LevelSetNodeStatusType M;
    LevelSetNodeListType* list = m_SparseLevelSet->GetListNode( status );
    while( !list->empty() )
    {
      p = list->front();
      idx = p.first;
      sparseNeighborhoodIt.SetLocation( idx );
      M = NumericTraits<LevelSetNodeStatusType>::Min();

      flag = true;
      for( typename SparseNeighborhoodIteratorType::Iterator i = sparseNeighborhoodIt.Begin();
        !i.IsAtEnd(); ++i )
      {
        q = i.Get();
        if ( q.m_Status == static_cast<LevelSetNodeStatusType>( status + 1 ) )
        {
          flag = false;
        }
        if ( ( M < q.m_Status ) && ( q.m_Status >= static_cast<LevelSetNodeStatusType>( status + 1 ) ) )
        {
          M = q.m_Status;
        }
      }

      if (!flag)
      {
        if ( status < m_MaxStatus )
        {
          m_StatusLists->GetListNode( status-1 )->push_back( p );
        }
        else
        {
          r.m_Status = m_MinStatus;
          r.m_Value = -3.0;
          m_SparseImage->SetPixel( idx, r );
        }
      }
      else
      {
        r.m_Status = p.m_Status;
        r.m_Value = static_cast< LevelSetNodeStatusType >( M-1 );
        m_SparseImage->SetPixel( idx, r );

        if ( r.m_Value >= static_cast<LevelSetOutputType>( o1 ) )
        {
          m_StatusLists->GetListNode( status+1 )->push_back( r );
        }
        if ( r.m_Value < static_cast<LevelSetOutputType>( o2 ) )
        {
          if ( status > m_MinStatus )
          {
            m_StatusLists->GetListNode( status-1 )->push_back( r );
          }
          else
          {
            r.m_Status = m_MinStatus;
            r.m_Value = -3.0;
            m_SparseImage->SetPixel( idx, r );
          }
        }
      }
      list->pop_front();
    }
  }


  void UpdatePlusLevelSet( LevelSetNodeStatusType& status )
  {
    LevelSetOutputType o1 = static_cast<LevelSetOutputType>(status) - 0.5;
    LevelSetOutputType o2 = static_cast<LevelSetOutputType>(status) + 0.5;

    bool flag;
    SparseImageIndexType idx;
    LevelSetNodePairType p;
    LevelSetNodeAttributeType q, r;
    ZeroFluxNeumannBoundaryCondition< SparseImageType > sp_nbc;

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

    LevelSetNodeStatusType M;
    LevelSetNodeListType* list_plus = m_SparseLevelSet->GetListNode( status );
    while( !list_plus->empty() )
    {
      p = list_plus->front();
      idx = p.first;
      sparseNeighborhoodIt.SetLocation( idx );
      M = NumericTraits<LevelSetNodeStatusType>::Max();

      flag = true;
      for( typename SparseNeighborhoodIteratorType::Iterator i = sparseNeighborhoodIt.Begin();
        !i.IsAtEnd(); ++i )
      {
        q = i.Get();
        if ( q.m_Status == static_cast<LevelSetNodeStatusType>( status - 1 ) )
        {
          flag = false;
        }
        if ( ( M > q.m_Status ) && ( q.m_Status <= static_cast<LevelSetNodeStatusType>( status - 1 ) ) )
        {
          M = q.m_Status;
        }
      }

      if (!flag)
      {
        if ( status < m_MaxStatus )
        {
          m_StatusLists->GetListNode( status+1 )->push_back( p );
        }
        else
        {
          r.m_Status = m_MaxStatus;
          r.m_Value = 3.0;
          m_SparseImage->SetPixel( idx, r );
        }
      }
      else
      {
        r.m_Status = p.m_Status;
        r.m_Value = static_cast< LevelSetNodeStatusType >( M+1 );
        m_SparseImage->SetPixel( idx, r );

        if ( r.m_Value <= static_cast<LevelSetOutputType>( o1 ) )
        {
          m_StatusLists->GetListNode( status-1 )->push_back( r );
        }
        if ( r.m_Value > static_cast<LevelSetOutputType>( o2 ) )
        {
          if ( status < m_MaxStatus )
          {
            m_StatusLists->GetListNode( status+1 )->push_back( r );
          }
          else
          {
            r.m_Status = m_MaxStatus;
            r.m_Value = 3.0;
            m_SparseImage->SetPixel( idx, r );
          }
        }
      }
      list_plus->pop_front();
    }
  }

  void GenerateData()
  {
    m_MinStatus = -3;
    m_MaxStatus = 3;

    UpdateZeroLevelSet();
    UpdateMinusLevelSet( -1 );
    UpdatePlusLevelSet( 1 );
    UpdateMinusLevelSet( -2 );
    UpdatePlusLevelSet( 2 );
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set get the input image
  itkSetObjectMacro( Update, UpdateListType );
  itkGetObjectMacro( Update, UpdateListType );

protected:
  UpdateSparseLevelSet() {}
  ~UpdateSparseLevelSet() {}

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
