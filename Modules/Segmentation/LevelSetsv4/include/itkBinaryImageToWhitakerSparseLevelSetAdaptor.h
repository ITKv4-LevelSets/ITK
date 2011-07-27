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

#ifndef __itkBinaryImageToWhitakerSparseLevelSetAdaptor_h
#define __itkBinaryImageToWhitakerSparseLevelSetAdaptor_h

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
template< class TInputImage, typename TLevelSetValueType >
class BinaryImageToWhitakerSparseLevelSetAdaptor : public Object
{
public:
  typedef BinaryImageToWhitakerSparseLevelSetAdaptor  Self;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;
  typedef Object                                      Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToWhitakerSparseLevelSetAdaptor, Object );

  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  typedef typename InputImageType::IndexType    InputImageIndexType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

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

  typedef typename LevelSetType::LabelMapType          LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer       LevelSetLabelMapPointer;

  typedef typename LevelSetType::LayerType             LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator         LevelSetLayerIterator;
  typedef typename LevelSetType::LayerConstIterator    LevelSetLayerConstIterator;

  typedef ImageRegionIteratorWithIndex< InputImageType >  InputIteratorType;
//  typedef ImageRegionIteratorWithIndex< SparseImageType > SparseIteratorType;
  typedef ShapedNeighborhoodIterator< InputImageType >    InputNeighborhoodIteratorType;
//  typedef ShapedNeighborhoodIterator< SparseImageType >   SparseNeighborhoodIteratorType;

  // this is the same as Procedure 1
  // Input is a binary image m_InputImage
  // Output is a WhitakerSparseLevelSetBasePointer
  void Initialize()
  {
    m_SparseLevelSet = LevelSetType::New();

    LevelSetLabelMapPointer labelMap = LevelSetLabelMapType::New();
    labelMap->SetBackgroundValue( 3 );
    labelMap->CopyInformation( m_InputImage );

    // Precondition labelmap and phi
    InputIteratorType iIt( m_InputImage, m_InputImage->GetRequestedRegion() );
    iIt.GoToBegin();

    while( !iIt.IsAtEnd() )
      {
      if ( iIt.Get() != NumericTraits< InputImagePixelType >::Zero )
        {
        labelMap->SetPixel( iIt.GetIndex(), -3 );
        }
      ++iIt;
      }

    labelMap->Optimize();

    m_SparseLevelSet->SetLabelMap( labelMap );

    FindActiveLayer();

    const LevelSetLayerType layer0 = m_SparseLevelSet->GetLayer( 0 );

    LevelSetLayerConstIterator nodeIt = layer0.begin();
    LevelSetLayerConstIterator nodeEnd = layer0.end();

    while( nodeIt != nodeEnd )
      {
      labelMap->SetPixel( nodeIt->first, 0 );
      ++nodeIt;
      }

    FindPlusOneMinusOneLayer();

    const LevelSetLayerType layerMinus1 = m_SparseLevelSet->GetLayer( -1 );

    nodeIt = layerMinus1.begin();
    nodeEnd = layerMinus1.end();

    while( nodeIt != nodeEnd )
      {
      labelMap->SetPixel( nodeIt->first, -1 );
      ++nodeIt;
      }

    const LevelSetLayerType layerPlus1 = m_SparseLevelSet->GetLayer( 1 );

    nodeIt = layerPlus1.begin();
    nodeEnd = layerPlus1.end();

    while( nodeIt != nodeEnd )
      {
      labelMap->SetPixel( nodeIt->first, 1 );
      ++nodeIt;
      }

    FindMinusTwoLayer();

    const LevelSetLayerType layerMinus2 = m_SparseLevelSet->GetLayer( -2 );

    nodeIt = layerMinus2.begin();
    nodeEnd = layerMinus2.end();

    while( nodeIt != nodeEnd )
      {
      labelMap->SetPixel( nodeIt->first, -2 );
      ++nodeIt;
      }

    FindPlusTwoLayer();

    const LevelSetLayerType layerPlus2 = m_SparseLevelSet->GetLayer( 2 );

    nodeIt = layerPlus2.begin();
    nodeEnd = layerPlus2.end();

    while( nodeIt != nodeEnd )
      {
      labelMap->SetPixel( nodeIt->first, 2 );
      ++nodeIt;
      }
  }

  void FindMinusTwoLayer()
  {
    LevelSetLabelMapPointer labelMap = m_SparseLevelSet->GetLabelMap();

    const LevelSetLayerType layer0 = m_SparseLevelSet->GetLayer( 0 );
    const LevelSetLayerType layerMinus1 = m_SparseLevelSet->GetLayer( -1 );

    LevelSetLayerType& layerMinus2 = m_SparseLevelSet->GetLayer( -2 );
    const LevelSetOutputType minus2 =
        static_cast< LevelSetOutputType >( -2 ) * NumericTraits< LevelSetOutputType >::One;

    LevelSetLayerConstIterator nodeIt = layerMinus1.begin();
    LevelSetLayerConstIterator nodeEnd = layerMinus1.end();

    typename InputImageType::SizeType imageSize = m_InputImage->GetRequestedRegion().GetSize();

    LevelSetInputType idx, tempIdx;

    while( nodeIt != nodeEnd )
    {
      idx = nodeIt->first;
      tempIdx = idx;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk += 2 )
        {
          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            tempIdx[dim] += kk;
            if( layerMinus1.find( tempIdx ) == layerMinus1.end() )
            {
              if( layer0.find( tempIdx ) == layer0.end() )
              {
//                if( labelObject->HasIndex( tempIdx ) )
                {
                  layerMinus2.insert(
                        std::pair< LevelSetInputType, LevelSetOutputType >( tempIdx, minus2 ) );
                }
              }
            }
          }
          tempIdx[dim] = idx[dim];
        }
      }
      ++nodeIt;
    }
  }

  void FindPlusTwoLayer()
  {
    LevelSetLabelMapPointer labelMap = m_SparseLevelSet->GetLabelMap();

    const LevelSetLayerType layer0 = m_SparseLevelSet->GetLayer( 0 );
    const LevelSetLayerType layerPlus1 = m_SparseLevelSet->GetLayer( 1 );

    LevelSetLayerType& layerPlus2 = m_SparseLevelSet->GetLayer( 2 );
    const LevelSetOutputType plus2 =
        static_cast< LevelSetOutputType >( 2 ) * NumericTraits< LevelSetOutputType >::One;

    LevelSetLayerConstIterator nodeIt = layerPlus1.begin();
    LevelSetLayerConstIterator nodeEnd = layerPlus1.end();

    typename InputImageType::SizeType imageSize = m_InputImage->GetRequestedRegion().GetSize();

    LevelSetInputType idx, tempIdx;

    while( nodeIt != nodeEnd )
    {
      idx = nodeIt->first;
      tempIdx = idx;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk+=2 )
        {
          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            tempIdx[dim] += kk;
            if( layerPlus1.find( tempIdx ) == layerPlus1.end() )
            {
              if( layer0.find( tempIdx ) == layer0.end() )
              {
//                if( !labelObject->HasIndex( tempIdx ) )
                {
                  layerPlus2.insert(
                        std::pair< LevelSetInputType, LevelSetOutputType >( tempIdx, plus2 ) );
                }
              }
            }
          }
          tempIdx[dim] = idx[dim];
        }
      }
      ++nodeIt;
    }
  }

  void FindActiveLayer()
  {
    LevelSetLabelMapPointer labelMap = m_SparseLevelSet->GetLabelMap();
    LevelSetLabelObjectPointer labelObject = labelMap->GetLabelObject( -3 );
    LevelSetLabelObjectPointer labelObject0 = LevelSetLabelObjectType::New();
    labelObject0->SetLabel( 0 );

    LevelSetLayerType& layer0 = m_SparseLevelSet->GetLayer( 0 );

    LevelSetLabelObjectLineContainerType
            lineContainer = labelObject->GetLineContainer();

    typename LevelSetLabelObjectLineContainerType::const_iterator
      lineIt = lineContainer.begin();

    const LevelSetOutputType zero = NumericTraits< LevelSetOutputType >::Zero;

    typename InputImageType::SizeType imageSize = m_InputImage->GetRequestedRegion().GetSize();

    while( lineIt != lineContainer.end() )
      {
      LevelSetLabelObjectLineType line = *lineIt;
      typename LevelSetLabelObjectLineType::IndexType
        index = line.GetIndex();

      typename LevelSetLabelObjectType::LengthType
        length = line.GetLength();

      typename LevelSetLabelObjectType::LengthType counter = 0;

      while( counter < length )
        {
        bool isOnBorder = false;
        typename LevelSetLabelObjectLineType::IndexType tempIdx = index;

        if( ( counter == 0 ) && ( tempIdx[0] != 0 ) )
          {
          --tempIdx[0];
          if( !labelObject->HasIndex( tempIdx ) )
            {
            isOnBorder = true;
            }
          tempIdx[0] = index[0];
          }
        if( isOnBorder )
        {
          layer0.insert(
              std::pair< LevelSetInputType, LevelSetOutputType >( index, zero ) );
          ++index[0];
          ++counter;
          continue;
        }

        if( ( counter == length-1 ) && ( tempIdx[0] != imageSize[0] ) )
          {
          ++tempIdx[0];
          if( !labelObject->HasIndex( tempIdx ) )
            {
            isOnBorder = true;
            }
          tempIdx[0] = index[0];
          }
        if( isOnBorder )
        {
          layer0.insert(
              std::pair< LevelSetInputType, LevelSetOutputType >( index, zero ) );
          ++index[0];
          ++counter;
          continue;
        }

        for( unsigned int dim = 1; ( dim < ImageDimension ) && !isOnBorder; dim++ )
        {
          for( int kk = -1; ( kk < 2 ) && !isOnBorder; kk += 2)
          {
            if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
            {
              tempIdx[dim] += kk;
              if( !labelObject->HasIndex( tempIdx ) )
                {
                isOnBorder = true;
                }
              tempIdx[dim] = index[dim];
            }
          }
        }
        if( isOnBorder )
          {
          layer0.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( index, zero ) );
          }

        ++index[0];
        ++counter;
        }
      ++lineIt;
      }
  }

  void FindPlusOneMinusOneLayer()
  {
    LevelSetLabelMapPointer labelMap = m_SparseLevelSet->GetLabelMap();
    LevelSetLabelObjectPointer labelObjectMinus3 = labelMap->GetLabelObject( -3 );

    const LevelSetOutputType minus1 = - NumericTraits< LevelSetOutputType >::One;
    const LevelSetOutputType plus1 = NumericTraits< LevelSetOutputType >::One;

    const LevelSetLayerType layer0 = m_SparseLevelSet->GetLayer( 0 );
    LevelSetLayerType& layerMinus1 = m_SparseLevelSet->GetLayer( -1 );
    LevelSetLayerType& layerPlus1 = m_SparseLevelSet->GetLayer( 1 );

    LevelSetLayerConstIterator nodeIt = layer0.begin();
    LevelSetLayerConstIterator nodeEnd = layer0.end();

    typename InputImageType::SizeType imageSize = m_InputImage->GetRequestedRegion().GetSize();

    LevelSetInputType idx, tempIdx;

    while( nodeIt != nodeEnd )
    {
      idx = nodeIt->first;
      tempIdx = idx;

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
        for( int kk = -1; kk < 2; kk+=2 )
        {
          if( ( tempIdx[dim] > 0 ) && ( tempIdx[dim] < imageSize[dim] - 1 ) )
          {
            tempIdx[dim] += kk;
            if( layer0.find( tempIdx ) == layer0.end() )
            {
              if( labelObjectMinus3->HasIndex( tempIdx ) )
              {
                layerMinus1.insert(
                      std::pair< LevelSetInputType, LevelSetOutputType >( tempIdx, minus1 ) );
              }
              else
              {
                layerPlus1.insert(
                      std::pair< LevelSetInputType, LevelSetOutputType >( tempIdx, plus1 ) );
              }
            }
          }
          tempIdx[dim] = idx[dim];
        }
      }
      ++nodeIt;
    }
  }

  // Set/Get the sparse levet set image
//  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set get the input image
  itkSetObjectMacro( InputImage, InputImageType );
  itkGetObjectMacro( InputImage, InputImageType );

protected:
  BinaryImageToWhitakerSparseLevelSetAdaptor() {}
  ~BinaryImageToWhitakerSparseLevelSetAdaptor() {}

  InputImagePointer  m_InputImage;
  LevelSetPointer    m_SparseLevelSet;
//  SparseImagePointer m_SparseImage;

private:
  BinaryImageToWhitakerSparseLevelSetAdaptor( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkBinaryImageToWhitakerSparseLevelSetAdaptor_h
