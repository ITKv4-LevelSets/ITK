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

#ifndef __itkBinaryImageToMalcolmSparseLevelSetAdaptor_h
#define __itkBinaryImageToMalcolmSparseLevelSetAdaptor_h

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
template< class TInputImage >
class BinaryImageToMalcolmSparseLevelSetAdaptor : public Object
{
public:
  typedef BinaryImageToMalcolmSparseLevelSetAdaptor   Self;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;
  typedef Object                                      Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( BinaryImageToMalcolmSparseLevelSetAdaptor, Object );

  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  typedef typename InputImageType::IndexType    InputImageIndexType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename NumericTraits< InputImagePixelType >::RealType
                                                InputPixelRealType;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       InputImageType::ImageDimension );

  typedef MalcolmSparseLevelSetBase< ImageDimension >   LevelSetType;
  typedef typename LevelSetType::Pointer                LevelSetPointer;
  typedef typename LevelSetType::InputType              LevelSetInputType;
  typedef typename LevelSetType::OutputType             LevelSetOutputType;

  typedef typename LevelSetType::SparseImageType        SparseImageType;
  typedef typename SparseImageType::Pointer             SparseImagePointer;

  typedef typename LevelSetType::NodeAttributeType      LevelSetNodeAttributeType;
  typedef typename LevelSetType::NodeStatusType         LevelSetNodeStatusType;
  typedef typename LevelSetType::NodePairType           LevelSetNodePairType;
  typedef typename LevelSetType::NodeListType           LevelSetNodeListType;
  typedef typename LevelSetType::NodeListIterator       LevelSetNodeListIterator;
  typedef typename LevelSetType::NodeListConstIterator  LevelSetNodeListConstIterator;

  typedef typename LevelSetType::SparseLayerMapType           SparseLayerMapType;
  typedef typename LevelSetType::SparseLayerMapIterator       SparseLayerMapIterator;
  typedef typename LevelSetType::SparseLayerMapConstIterator  SparseLayerMapConstIterator;

  typedef ImageRegionIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef ImageRegionIteratorWithIndex< SparseImageType > SparseIteratorType;
  typedef ShapedNeighborhoodIterator< InputImageType >    InputNeighborhoodIteratorType;
  typedef ShapedNeighborhoodIterator< SparseImageType >   SparseNeighborhoodIteratorType;

  // this is the same as Procedure 1
  // Input is a binary image m_InputImage
  // Output is a WhitakerSparseLevelSetBasePointer
  void Initialize()
  {
    std::cout << "Initialize()" << std::endl;
    m_SparseLevelSet = LevelSetType::New();

    m_SparseImage = SparseImageType::New();
    m_SparseImage->CopyInformation( m_InputImage );
    m_SparseImage->SetLargestPossibleRegion( m_InputImage->GetLargestPossibleRegion() );
    m_SparseImage->SetBufferedRegion( m_InputImage->GetBufferedRegion() );
    m_SparseImage->SetRequestedRegion( m_InputImage->GetRequestedRegion() );
    m_SparseImage->Allocate();
    m_SparseLevelSet->SetImage( m_SparseImage );
    std::cout << "Allocated sparse image" << std::endl;

    // Precondition labelmap and phi
    SparseIteratorType sIt( m_SparseImage, m_SparseImage->GetRequestedRegion() );
    sIt.GoToBegin();
    InputIteratorType iIt( m_InputImage, m_InputImage->GetRequestedRegion() );
    iIt.GoToBegin();

    LevelSetNodeAttributeType p_plus1;
    p_plus1.m_Status = 1;
    p_plus1.m_Value = 1;

    LevelSetNodeAttributeType p_minus1;
    p_minus1.m_Status = -1;
    p_minus1.m_Value = -1;

    while( !iIt.IsAtEnd() )
      {
      if ( iIt.Get() == NumericTraits< InputImagePixelType >::Zero )
        {
        sIt.Set( p_plus1 );
        }
      else
        {
        sIt.Set( p_minus1 );
        }
      ++iIt;
      ++sIt;
      }

    FindActiveLayer();
  }

  void FindActiveLayer()
  {
    LevelSetNodePairType nodePair;

    typename InputNeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    ZeroFluxNeumannBoundaryCondition< InputImageType > im_nbc;

    // Find the zero levelset
    typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType > InputFaceCalculatorType;
    InputFaceCalculatorType faceCalculator;
    typename InputFaceCalculatorType::FaceListType faceList;
    faceList = faceCalculator(m_InputImage, m_InputImage->GetLargestPossibleRegion(), radius);
    typename InputFaceCalculatorType::FaceListType::iterator faceListIterator;

    typename InputNeighborhoodIteratorType::OffsetType temp_offset;
    temp_offset.Fill( 0 );

    for ( faceListIterator=faceList.begin(); faceListIterator != faceList.end(); ++faceListIterator)
    {
      InputNeighborhoodIteratorType
      inputNeighborhoodIt( radius, m_InputImage,
                           *faceListIterator );
                           inputNeighborhoodIt.OverrideBoundaryCondition( &im_nbc );

      inputNeighborhoodIt.GoToBegin();
      SparseIteratorType sIt( m_SparseImage, *faceListIterator );
      sIt.GoToBegin();
      InputIteratorType iIt( m_InputImage, *faceListIterator );
      iIt.GoToBegin();

      for( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
        temp_offset[dim] = -1;
        inputNeighborhoodIt.ActivateOffset( temp_offset );
        temp_offset[dim] = 1;
        inputNeighborhoodIt.ActivateOffset( temp_offset );
        temp_offset[dim] = 0;
        }

      LevelSetNodeAttributeType p_zero;
      p_zero.m_Status = 0;
      p_zero.m_Value = static_cast< LevelSetOutputType >( 0.0 );

      while( !iIt.IsAtEnd() )
        {
        bool flag = false;

        // Iterate through all the pixels in the neighborhood
        for( typename InputNeighborhoodIteratorType::Iterator i = inputNeighborhoodIt.Begin();
        !i.IsAtEnd(); ++i )
        {
        if( i.Get() == NumericTraits< InputImagePixelType >::Zero )
          {
          flag = true;
          break;
          }
        }

        if ( ( iIt.Get() != NumericTraits< InputImagePixelType >::Zero ) && ( flag ) )
          {
          nodePair.first = inputNeighborhoodIt.GetIndex();
          nodePair.second = p_zero;
          m_SparseLevelSet->GetListNode( 0 )->push_back( nodePair );
          sIt.Set( p_zero );
          }
        ++inputNeighborhoodIt;
        ++sIt;
        ++iIt;
      }
    }
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  // Set get the input image
  itkSetObjectMacro( InputImage, InputImageType );
  itkGetObjectMacro( InputImage, InputImageType );

protected:
  BinaryImageToMalcolmSparseLevelSetAdaptor() {}
  ~BinaryImageToMalcolmSparseLevelSetAdaptor() {}

  InputImagePointer  m_InputImage;
  LevelSetPointer    m_SparseLevelSet;
  SparseImagePointer m_SparseImage;

private:
  BinaryImageToMalcolmSparseLevelSetAdaptor( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkBinaryImageToMalcolmSparseLevelSetAdaptor_h
