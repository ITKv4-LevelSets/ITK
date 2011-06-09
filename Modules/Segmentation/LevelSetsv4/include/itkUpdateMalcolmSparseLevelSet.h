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

#ifndef __itkUpdateMalcolmSparseLevelSet_h
#define __itkUpdateMalcolmSparseLevelSet_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkMalcolmSparseLevelSetBase.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <list>
#include "itkObject.h"

namespace itk
{
template< unsigned int VDimension >
class UpdateMalcolmSparseLevelSet : public Object
{
public:
  typedef UpdateMalcolmSparseLevelSet   Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( UpdateMalcolmSparseLevelSet, Object );

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef MalcolmSparseLevelSetBase< ImageDimension >  LevelSetType;
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

  typedef ImageRegionIteratorWithIndex< SparseImageType > SparseIteratorType;
  typedef ShapedNeighborhoodIterator< SparseImageType >   SparseNeighborhoodIteratorType;

  void UnPhasedPropagation()
    {
    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_0 = m_SparseLevelSet->GetListNode();

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
    LevelSetOutputRealType update;

    while( !m_Update->empty() )
      {
      update = m_Update->front();
      m_Update->pop_front();

      p = list_0->front();
      list_0->pop_front();

      if( update != NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( update > NumericTraits< LevelSetOutputRealType >::Zero )
          {
          p.second = -1;
          }
        if( update < NumericTraits< LevelSetOutputRealType >::Zero )
          {
          p.second = 1;
          }
        m_SparseImage->SetPixel( p.first, p.second );
        sparseNeighborhoodIt.SetLocation( p.first );

        for( typename SparseNeighborhoodIteratorType::Iterator
                i = sparseNeighborhoodIt.Begin();
            !i.IsAtEnd(); ++i )
          {
          q = i.Get();

          if( q * p.second == -1 )
            {
            LevelSetNodePairType temp;
            temp.first =
              sparseNeighborhoodIt.GetIndex( i.GetNeighborhoodOffset() );
            temp.second = 0;
            new_list_0.push_back( temp );
            m_SparseImage->SetPixel( temp.first, temp.second );
            }
          }
        }
      else
        {
        new_list_0.push_back( p );
        }
      }

    while( !new_list_0.empty() )
      {
      list_0->push_back( new_list_0.front() );
      new_list_0.pop_front();
      }
    }

  void PhasedPropagation( LevelSetNodeListType& ioList,
                         UpdateListType& ioUpdate,
                         const bool& iContraction )
    {
    if( ioList.size() != ioUpdate.size() )
      {
      itkGenericExceptionMacro( "ioList and ioUpdate have different size!" );
      }

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
    LevelSetOutputRealType update;

    while( !ioUpdate.empty() )
      {
      update = ioUpdate.front();
      ioUpdate.pop_front();

      p = ioList.front();
      ioList.pop_front();

      bool to_be_updated = false;

      if( update != NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( iContraction ) // contraction
          {
          // only allow positive forces
          if( update > NumericTraits< LevelSetOutputRealType >::Zero )
            {
            p.second = -1;
            to_be_updated = true;
            }
          }
        else // Dilation
          {
          // only allow negative forces
          if( update < NumericTraits< LevelSetOutputRealType >::Zero )
            {
            p.second = 1;
            to_be_updated = true;
            }
          }
        if( to_be_updated )
          {
          m_SparseImage->SetPixel( p.first, p.second );
          sparseNeighborhoodIt.SetLocation( p.first );

          for( typename SparseNeighborhoodIteratorType::Iterator
                  i = sparseNeighborhoodIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            q = i.Get();

            if( q * p.second == -1 )
              {
              LevelSetNodePairType temp;
              temp.first =
                sparseNeighborhoodIt.GetIndex( i.GetNeighborhoodOffset() );
              temp.second = 0;
              m_StatusLists->GetListNode()->push_back( temp );
              m_SparseImage->SetPixel( temp.first, temp.second );
              }
            }
          }
        else
          {
          m_StatusLists->GetListNode()->push_back( p );
          }
        }
      else
        {
        m_StatusLists->GetListNode()->push_back( p );
        }
      }
    }

  void MinimalInterface()
    {
    LevelSetNodeListType new_list_0;
    LevelSetNodeListType* list_0 = m_SparseLevelSet->GetListNode();

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

    while( !list_0->empty() )
      {
      p = list_0->front();
      list_0->pop_front();

      sparseNeighborhoodIt.SetLocation( p.first );

      bool positive = false;
      bool negative = false;

      for( typename SparseNeighborhoodIteratorType::Iterator
          i = sparseNeighborhoodIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        q = i.Get();

        if( q != NumericTraits< LevelSetOutputType >::Zero )
          {
          if( q == -1 )
            {
            negative = true;
            if( positive )
              {
              break;
              }
            }
          else // ( q == 1 )
            {
            positive = true;
            if( negative )
              {
              break;
              }
            }
          }
        }
      if( negative && !positive )
        {
        p.second = -1;
        m_SparseImage->SetPixel( p.first, p.second );
        }
      else
        {
        if( positive && !negative )
          {
          p.second = 1;
          m_SparseImage->SetPixel( p.first, p.second );
          }
        else
          {
          new_list_0.push_back( p );
          }
        }
      }

    while( !new_list_0.empty() )
      {
      list_0->push_back( new_list_0.front() );
      new_list_0.pop_front();
      }
    }

  void Update()
  {
    if( m_SparseLevelSet.IsNull() )
      {
      itkGenericExceptionMacro( <<"m_SparseLevelSet is NULL" );
      }
    if( m_Update->empty() )
      {
      itkGenericExceptionMacro( <<"m_Update is empty" );
      }
    m_SparseImage = m_SparseLevelSet->GetImage();

    if( m_UnPhased )
      {
      UnPhasedPropagation();
      MinimalInterface();
      }
    else
      {
      LevelSetNodeListType* list_0 = m_SparseLevelSet->GetListNode();

      LevelSetNodeListType list_pos;
      UpdateListType update_pos;

      LevelSetNodeListType list_neg;
      UpdateListType update_neg;

      LevelSetNodePairType p;
      LevelSetOutputRealType update;

      while( !m_Update->empty() )
        {
        p = list_0->front();
        list_0->pop_front();

        update = m_Update->front();
        m_Update->pop_front();

        if( update > 0 )
          {
          list_pos.push_back( p );
          update_pos.push_back( update );
          }
        else
          {
          list_neg.push_back( p );
          update_neg.push_back( update );
          }
        }

      // contraction
      PhasedPropagation( list_pos, update_pos, true );
      MinimalInterface();

      // dilation
      PhasedPropagation( list_neg, update_neg, false );
      MinimalInterface();

      while( !m_StatusLists->GetListNode()->empty() )
        {
        m_SparseLevelSet->GetListNode()->push_back(
              m_StatusLists->GetListNode()->front() );
        m_StatusLists->GetListNode()->pop_front();
        }

      MinimalInterface();
      }
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  void SetUpdate( UpdateListType* & iUpdate )
    {
    if( !iUpdate->empty() )
      {
      m_Update = iUpdate;
      }
    }

protected:
  UpdateMalcolmSparseLevelSet() : m_UnPhased( false )
    {
    m_StatusLists = LevelSetType::New();
    }

  ~UpdateMalcolmSparseLevelSet() {}

  UpdateListType* m_Update;

  LevelSetPointer    m_StatusLists;
  LevelSetPointer    m_SparseLevelSet;
  SparseImagePointer m_SparseImage;

  bool m_UnPhased;

private:
  UpdateMalcolmSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateMalcolmSparseLevelSet_h
