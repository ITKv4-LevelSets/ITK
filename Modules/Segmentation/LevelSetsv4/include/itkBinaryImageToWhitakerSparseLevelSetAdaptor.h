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
  typedef typename LevelSetType::LabelObjectLengthType LevelSetLabelObjectLengthType;
  typedef typename LevelSetType::LabelObjectLineType   LevelSetLabelObjectLineType;
  typedef typename LevelSetType::LabelObjectLineContainerType
                                                       LevelSetLabelObjectLineContainerType;

  typedef typename LevelSetType::LabelMapType          LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer       LevelSetLabelMapPointer;

  typedef typename LevelSetType::LayerType             LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator         LevelSetLayerIterator;
  typedef typename LevelSetType::LayerConstIterator    LevelSetLayerConstIterator;

  // this is the same as Procedure 1
  // Input is a binary image m_InputImage
  // Output is a WhitakerSparseLevelSetBasePointer
  void Initialize()
  {
    m_LabelMap = LevelSetLabelMapType::New();
    m_LabelMap->SetBackgroundValue( static_cast< char >( 3 ) );
    m_LabelMap->CopyInformation( m_InputImage );

    m_InternalImage = InternalImageType::New();
    m_InternalImage->CopyInformation( m_InputImage );
    m_InternalImage->SetRegions( m_InputImage->GetBufferedRegion() );
    m_InternalImage->Allocate();
    m_InternalImage->FillBuffer( static_cast< char >( 3 ) );

    LevelSetLabelObjectPointer InnerPart = LevelSetLabelObjectType::New();
    InnerPart->SetLabel( static_cast< char >( -3 ) );

    // Precondition labelmap and phi
    InputIteratorType iIt( m_InputImage, m_InputImage->GetLargestPossibleRegion() );
    iIt.GoToBegin();

    InternalIteratorType labelIt( m_InternalImage,
                                  m_InternalImage->GetLargestPossibleRegion() );
    labelIt.GoToBegin();

    while( !iIt.IsAtEnd() )
      {
      if ( iIt.Get() != NumericTraits< InputImagePixelType >::Zero )
        {
        InnerPart->AddIndex( iIt.GetIndex() );
        labelIt.Set( static_cast< char >( -3 ) );
        }
      ++labelIt;
      ++iIt;
      }

    InnerPart->Optimize();
    m_LabelMap->AddLabelObject( InnerPart );

    FindActiveLayer();

    FindPlusOneMinusOneLayer();

    FindMinusTwoLayer();

    FindPlusTwoLayer();

    m_InternalImage = 0;

    m_SparseLevelSet->SetLabelMap( m_LabelMap );
  }

  // Set/Get the sparse levet set image
//  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set get the input image
  itkSetObjectMacro( InputImage, InputImageType );
  itkGetObjectMacro( InputImage, InputImageType );

protected:
  BinaryImageToWhitakerSparseLevelSetAdaptor()
    {
    m_SparseLevelSet = LevelSetType::New();
    }
  ~BinaryImageToWhitakerSparseLevelSetAdaptor() {}

  InputImagePointer       m_InputImage;
  LevelSetPointer         m_SparseLevelSet;
  LevelSetLabelMapPointer m_LabelMap;

  typedef Image< char, ImageDimension >         InternalImageType;
  typedef typename InternalImageType::Pointer   InternalImagePointer;

  InternalImagePointer m_InternalImage;

  typedef ImageRegionIteratorWithIndex< InputImageType >      InputIteratorType;
  typedef ImageRegionIteratorWithIndex< InternalImageType >   InternalIteratorType;

  typedef ShapedNeighborhoodIterator< InternalImageType > NeighborhoodIteratorType;

  void FindMinusTwoLayer()
  {
    const LevelSetLayerType layerMinus1 = m_SparseLevelSet->GetLayer( static_cast< char >( -1 ) );

    LevelSetLayerType& layerMinus2 = m_SparseLevelSet->GetLayer( static_cast< char >( -2 ) );
    const LevelSetOutputType minus2 = static_cast< LevelSetOutputType >( -2 );

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    ZeroFluxNeumannBoundaryCondition< InternalImageType > im_nbc;

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &im_nbc );

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


    LevelSetLayerConstIterator nodeIt = layerMinus1.begin();
    LevelSetLayerConstIterator nodeEnd = layerMinus1.end();

    while( nodeIt != nodeEnd )
    {
      LevelSetInputType idx = nodeIt->first;
      neighIt.SetLocation( idx );

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        if( it.Get() == static_cast< char >( -3 ) )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          layerMinus2.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, minus2 ) );
          }
        }
      ++nodeIt;
    }

    LevelSetLabelObjectPointer ObjectMinus2 = LevelSetLabelObjectType::New();
    ObjectMinus2->SetLabel(static_cast< char >( -2 ) );

    nodeIt = layerMinus2.begin();
    nodeEnd = layerMinus2.end();

    while( nodeIt != nodeEnd )
      {
      m_LabelMap->GetLabelObject( static_cast< char >( -3 ) )->RemoveIndex( nodeIt->first );
      ObjectMinus2->AddIndex( nodeIt->first );
      m_InternalImage->SetPixel( nodeIt->first, -2 );
      ++nodeIt;
      }

    ObjectMinus2->Optimize();
    m_LabelMap->AddLabelObject( ObjectMinus2 );
  }

  void FindPlusTwoLayer()
  {
    const LevelSetLayerType layerPlus1 = m_SparseLevelSet->GetLayer( static_cast< char >( 1 ) );

    LevelSetLayerType& layerPlus2 = m_SparseLevelSet->GetLayer( static_cast< char >( 2 ) );
    const LevelSetOutputType plus2 = static_cast< LevelSetOutputType >( 2 );

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    ZeroFluxNeumannBoundaryCondition< InternalImageType > im_nbc;

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &im_nbc );

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


    LevelSetLayerConstIterator nodeIt = layerPlus1.begin();
    LevelSetLayerConstIterator nodeEnd = layerPlus1.end();

    while( nodeIt != nodeEnd )
    {
      LevelSetInputType idx = nodeIt->first;
      neighIt.SetLocation( idx );

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        if( it.Get() == static_cast< char >( 3 ) )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          m_InternalImage->SetPixel( nodeIt->first, static_cast< char >( 2 ) );

          layerPlus2.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, plus2 ) );
          }
        }
      ++nodeIt;
    }

    LevelSetLabelObjectPointer ObjectPlus2 = LevelSetLabelObjectType::New();
    ObjectPlus2->SetLabel( static_cast< char >( 2 ) );

    nodeIt = layerPlus2.begin();
    nodeEnd = layerPlus2.end();

    while( nodeIt != nodeEnd )
      {
      ObjectPlus2->AddIndex( nodeIt->first );
      m_InternalImage->SetPixel( nodeIt->first, static_cast< char >( 2 ) );
      ++nodeIt;
      }

    ObjectPlus2->Optimize();
    m_LabelMap->AddLabelObject( ObjectPlus2 );
  }

  void FindActiveLayer()
  {
    LevelSetLabelObjectPointer labelObject = m_LabelMap->GetLabelObject( -3 );

    LevelSetLabelObjectPointer labelObject0 = LevelSetLabelObjectType::New();
    labelObject0->SetLabel( static_cast< char >( 0 ) );

    const LevelSetOutputType zero = NumericTraits< LevelSetOutputType >::Zero;

    LevelSetLayerType& layer0 = m_SparseLevelSet->GetLayer( static_cast< char >( 0 ) );

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    ZeroFluxNeumannBoundaryCondition< InternalImageType > im_nbc;

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &im_nbc );

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

    LevelSetLabelObjectLineContainerType
            lineContainer = labelObject->GetLineContainer();

    typename LevelSetLabelObjectLineContainerType::const_iterator
      lineIt = lineContainer.begin();
    typename LevelSetLabelObjectLineContainerType::const_iterator
      lineEnd = lineContainer.end();

    while( lineIt != lineEnd )
      {
      const LevelSetInputType &     firstIdx = lineIt->GetIndex();
      LevelSetLabelObjectLengthType length = lineIt->GetLength();

      typedef typename LevelSetInputType::IndexValueType IndexValueType;
      IndexValueType endIdx0 = firstIdx[0] + static_cast< IndexValueType >( length );

      for( LevelSetInputType idx = firstIdx; idx[0] < endIdx0; ++idx[0] )
        {
        neighIt.SetLocation( idx );

        for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
             !it.IsAtEnd();
             ++it )
          {
          if( it.Get() == static_cast< char >( 3 ) )
            {
            layer0.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( idx, zero ) );
            break;
            }
          }
        }

      ++lineIt;
      }

    if( !layer0.empty() )
      {
      LevelSetLabelObjectPointer ZeroSet = LevelSetLabelObjectType::New();
      ZeroSet->SetLabel( static_cast< char >( 0 ) );

      LevelSetLayerConstIterator nodeIt = layer0.begin();
      LevelSetLayerConstIterator nodeEnd = layer0.end();

      while( nodeIt != nodeEnd )
        {
        m_LabelMap->GetLabelObject( static_cast< char >( -3 ) )->RemoveIndex( nodeIt->first );
        ZeroSet->AddIndex( nodeIt->first );
        m_InternalImage->SetPixel( nodeIt->first, 0 );
        ++nodeIt;
        }

      ZeroSet->Optimize();
      m_LabelMap->AddLabelObject( ZeroSet );
      }
  }

  void FindPlusOneMinusOneLayer()
  {
    const LevelSetOutputType minus1 = - NumericTraits< LevelSetOutputType >::One;
    const LevelSetOutputType plus1 = NumericTraits< LevelSetOutputType >::One;

    const LevelSetLayerType layer0 = m_SparseLevelSet->GetLayer( static_cast< char >( 0 ) );
    LevelSetLayerType& layerMinus1 = m_SparseLevelSet->GetLayer( static_cast< char >( -1 ) );
    LevelSetLayerType& layerPlus1 = m_SparseLevelSet->GetLayer( static_cast< char >( 1 ) );

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    ZeroFluxNeumannBoundaryCondition< InternalImageType > im_nbc;

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &im_nbc );

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

    LevelSetLayerConstIterator nodeIt   = layer0.begin();
    LevelSetLayerConstIterator nodeEnd  = layer0.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType idx = nodeIt->first;
      neighIt.SetLocation( idx );

      for( typename NeighborhoodIteratorType::Iterator it = neighIt.Begin();
           !it.IsAtEnd();
           ++it )
        {
        if( it.Get() == static_cast< char >( -3 ) )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          layerMinus1.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, minus1 ) );
          }
        if( it.Get() == static_cast< char >( 3 ) )
          {
          LevelSetInputType tempIndex =
              neighIt.GetIndex( it.GetNeighborhoodOffset() );

          layerPlus1.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, plus1 ) );
          }
        }
      ++nodeIt;
      }

    LevelSetLabelObjectPointer ObjectMinus1 = LevelSetLabelObjectType::New();
    ObjectMinus1->SetLabel( static_cast< char >( -1 ) );

    nodeIt = layerMinus1.begin();
    nodeEnd = layerMinus1.end();

    while( nodeIt != nodeEnd )
      {
      m_LabelMap->GetLabelObject( static_cast< char >( -3 ) )->RemoveIndex( nodeIt->first );
      ObjectMinus1->AddIndex( nodeIt->first );
      m_InternalImage->SetPixel( nodeIt->first, static_cast< char >( -1 ) );
      ++nodeIt;
      }

    ObjectMinus1->Optimize();
    m_LabelMap->AddLabelObject( ObjectMinus1 );

    LevelSetLabelObjectPointer ObjectPlus1 = LevelSetLabelObjectType::New();
    ObjectPlus1->SetLabel( static_cast< char >( 1 ) );

    nodeIt = layerPlus1.begin();
    nodeEnd = layerPlus1.end();

    while( nodeIt != nodeEnd )
      {
      ObjectPlus1->AddIndex( nodeIt->first );
      m_InternalImage->SetPixel( nodeIt->first, static_cast< char >( -1 ) );
      ++nodeIt;
      }

    ObjectPlus1->Optimize();
    m_LabelMap->AddLabelObject( ObjectPlus1 );
  }

private:
  BinaryImageToWhitakerSparseLevelSetAdaptor( const Self& );
  void operator = ( const Self& );


};
}
#endif // __itkBinaryImageToWhitakerSparseLevelSetAdaptor_h
