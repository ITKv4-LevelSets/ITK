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


#ifndef __itkBinaryImageToSparseLevelSetAdaptor_h
#define __itkBinaryImageToSparseLevelSetAdaptor_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkSparseLevelSetBase.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include <list>
#include "itkObject.h"

namespace itk
{
template< class TInputImage, class TLevelSetType >
class BinaryImageToSparseLevelSetAdaptor : public Object
{
public:
  typedef BinaryImageToSparseLevelSetAdaptor  Self;
  typedef SmartPointer< Self >                Pointer;
  typedef SmartPointer< const Self >          ConstPointer;
  typedef Object                              Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToSparseLevelSetAdaptor, Object );

  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  typedef typename InputImageType::IndexType    InputImageIndexType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

  typedef TLevelSetType                                LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::ImageType             LevelSetImageType;
  typedef typename LevelSetType::InputType             LevelSetInputType;
  typedef typename LevelSetType::OutputType            LevelSetOutputType;
  typedef typename LevelSetImageType::Pointer          LevelSetImagePointer;

  typedef typename LevelSetType::SparseImageType       SparseImageType;
  typedef typename SparseImageType::Pointer            SparseImagePointer;

  typedef typename LevelSetType::NodeStatusType        LevelSetNodeStatusType;
  typedef typename LevelSetType::NodePairType          LevelSetNodePairType;
  typedef typename LevelSetType::NodeListType          LevelSetNodeListType;
  typedef typename LevelSetType::NodeListIterator      LevelSetNodeListIterator;
  typedef typename LevelSetType::NodeListConstIterator LevelSetNodeListConstIterator;

  typedef typename LevelSetType::SparseLayerMapType           SparseLayerMapType;
  typedef typename LevelSetType::SparseLayerMapIterator       SparseLayerMapIterator;
  typedef typename LevelSetType::SparseLayerMapConstIterator  SparseLayerMapConstIterator;

  typedef ImageRegionIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef ImageRegionIteratorWithIndex< SparseImageType > SparseIteratorType;
  typedef NeighborhoodIterator< InputImageType >          InputNeighborhoodIteratorType;
  typedef NeighborhoodIterator< SparseImageType >         SparseNeighborhoodIteratorType;

  // this is the same as Procedure 1
  // Input is a binary image m_InputImage
  // Output is a WhitakerSparseLevelSetBasePointer
  void Initialization()
  {
    m_SparseLevelSet = LevelSetType::New();

    SparseImagePointer sparseImage = m_SparseLevelSet->GetImage();
    sparseImage->CopyInformation( m_InputImage );
    sparseImage->SetLargestPossibleRegion( m_InputImage->GetLargestPossibleRegion() );
    sparseImage->SetBufferedRegion( m_InputImage->GetBufferedRegion() );
    sparseImage->SetRequestedRegion( m_InputImage->GetRequestedRegion() );
    sparseImage->Allocate();

    // Precondition labelmap and phi
    LevelSetNodePairType p, q;
    SparseIteratorType sIt( sparseImage, sparseImage->GetRequestedRegion() );
    sIt.GoToBegin();
    InputIteratorType iIt( m_InputImage, m_InputImage->GetRequestedRegion() );
    iIt.GoToBegin();
    while( !iIt.IsAtEnd() )
    {
      if ( iIt.Get() == 0 )
      {
        p.first( 3 );
        p.second( 3.0 );
        sIt.Set( p );
      }

      if ( iIt.Get() == 1 )
      {
        p.first( -3 );
        p.second( -3.0 );
        sIt.Set( p );
      }
      ++iIt;
      ++sIt;
    }

    // Find the zero levelset
    bool flag = false;
    InputNeighborhoodIteratorType inputNeighborhoodIt ( 1, m_InputImage, m_InputImage->GetRequestedRegion() );
    inputNeighborhoodIt.GoToBegin();
    sIt.GoToBegin();
    while( !inputNeighborhoodIt.IsAtEnd() )
    {
      flag = false;

      // Iterate through all the pixels in the neighborhood
      for( unsigned int i = 0; i < inputNeighborhoodIt.Size(); i++ )
      {
        if( inputNeighborhoodIt.GetPixel( i ) == 0 )
        {
          flag = true;
          break;
        }
      }

      if ( ( inputNeighborhoodIt.GetCenterPixel() == 1 ) && ( flag ) )
      {
        p.first( 0 );
        p.second( 0.0 );
        m_SparseLevelSet->GetListNode( 0 )->push_back( p );
        sIt.Set( p );
      }
      ++inputNeighborhoodIt;
      ++sIt;
    }

    // Find the +1 and -1 levelset
    InputImageIndexType idx;
    LevelSetNodeListType* list_of_nodes = m_SparseLevelSet->GetListNode( 0 );
    LevelSetNodeListIterator node_it = list_of_nodes->begin();
    LevelSetNodeListIterator node_end = list_of_nodes->end();
    SparseNeighborhoodIteratorType sparseNeighborhoodIt ( 1, sparseImage, sparseImage->GetRequestedRegion() );
    while( node_it != node_end )
    {
      idx = (*node_it)->first;
      sparseNeighborhoodIt.SetLocation( idx );

      // Iterate through all the pixels in the neighborhood
      for( unsigned int i = 0; i < sparseNeighborhoodIt.Size(); i++ )
      {
        q = sparseNeighborhoodIt.GetPixel( i );
        if ( q.first == -3 )
        {
          q.first = -1;
          q.second = -1.0;
          m_SparseLevelSet->GetListNode( -1 )->push_back( q );
        }

        if ( q.first == 3 )
        {
          q.first = 1;
          q.second = 1.0;
          m_SparseLevelSet->GetListNode( 1 )->push_back( q );
        }
      }
      ++node_it;
    }

    // Find the +2 and -2 levelset
    list_of_nodes = m_SparseLevelSet->GetListNode( -1 );
    node_it = list_of_nodes->begin();
    node_end = list_of_nodes->end();
    while( node_it != node_end )
    {
      idx = (*node_it)->first;
      sparseNeighborhoodIt.SetLocation( idx );

      // Iterate through all the pixels in the neighborhood
      for( unsigned int i = 0; i < sparseNeighborhoodIt.Size(); i++ )
      {
        q = sparseNeighborhoodIt.GetPixel( i );
        if ( q.first == -3 )
        {
          q.first = -2;
          q.second = -2.0;
          m_SparseLevelSet->GetListNode( -2 )->push_back( q );
        }
      }
      ++node_it;
    }

    list_of_nodes = m_SparseLevelSet->GetListNode( 1 );
    node_it = list_of_nodes->begin();
    node_end = list_of_nodes->end();
    while( node_it != node_end )
    {
      idx = (*node_it)->first;
      sparseNeighborhoodIt.SetLocation( idx );

      // Iterate through all the pixels in the neighborhood
      for( unsigned int i = 0; i < sparseNeighborhoodIt.Size(); i++ )
      {
        q = sparseNeighborhoodIt.GetPixel( i );
        if ( q.first == 3 )
        {
          q.first = 2;
          q.second = 2.0;
          m_SparseLevelSet->GetListNode( 2 )->push_back( q );
        }
      }
      ++node_it;
    }
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set get the input image
  itkSetObjectMacro( InputImage, InputImageType );
  itkGetObjectMacro( InputImage, InputImageType );

protected:
  BinaryImageToSparseLevelSetAdaptor() {}
  ~BinaryImageToSparseLevelSetAdaptor() {}

  InputImagePointer m_InputImage;
  LevelSetPointer   m_SparseLevelSet;

private:
  BinaryImageToSparseLevelSetAdaptor( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkBinaryImageToSparseLevelSetAdaptor_h
