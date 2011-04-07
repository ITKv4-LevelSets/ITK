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
#ifndef __itkLevelSetDomainPartitionBase_h
#define __itkLevelSetDomainPartitionBase_h

#include "itkLightObject.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTreeGenerator.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
/** \class LevelSetDomainPartitionBase
 *
 * \brief Helper class used to partition domain and efficiently compute overlap.
 *
 */
template< class TInputImage, class TFeatureImage >
class LevelSetDomainPartitionBase : public LightObject
{
public:

  typedef LevelSetDomainPartitionBase             Self;
  typedef LightObject                           Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, TFeatureImage::ImageDimension);

  itkTypeMacro(LevelSetDomainPartitionBase, LightObject);

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      InputPixelType;
  typedef typename InputImageType::RegionType     InputRegionType;
  typedef typename InputImageType::SizeType       InputSizeType;
  typedef typename InputSizeType::SizeValueType   InputSizeValueType;
  typedef typename InputImageType::SpacingType    InputSpacingType;
  typedef typename InputImageType::IndexType      InputIndexType;
  typedef typename InputIndexType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::PointType      InputPointType;

  typedef TFeatureImage                           FeatureImageType;
  typedef typename FeatureImageType::Pointer      FeatureImagePointer;
  typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;
  typedef typename FeatureImageType::PixelType    FeaturePixelType;
  typedef typename FeatureImageType::RegionType   FeatureRegionType;
  typedef typename FeatureImageType::SizeType     FeatureSizeType;
  typedef typename FeatureSizeType::SizeValueType FeatureSizeValueType;
  typedef typename FeatureImageType::SpacingType  FeatureSpacingType;
  typedef typename FeatureImageType::IndexType    FeatureIndexType;
  typedef typename FeatureImageType::PointType    FeaturePointType;

  typedef std::list< unsigned int > ListPixelType;
  typedef Image< ListPixelType, itkGetStaticConstMacro(ImageDimension) >
  ListImageType;
  typedef typename ListImageType::Pointer               ListImagePointer;
  typedef typename ListImageType::ConstPointer          ListImageConstPointer;
  typedef typename ListImageType::RegionType            ListRegionType;
  typedef typename ListImageType::SizeType              ListSizeType;
  typedef typename ListSizeType::SizeValueType          ListSizeValueType;
  typedef typename ListImageType::SpacingType           ListSpacingType;
  typedef typename ListImageType::IndexType             ListIndexType;
  typedef typename ListIndexType::IndexValueType        ListIndexValueType;
  typedef typename ListImageType::PointType             ListPointType;
  typedef ImageRegionIteratorWithIndex< ListImageType > ListIteratorType;

  typedef Vector< float, itkGetStaticConstMacro(ImageDimension) >
    CentroidVectorType;
  typedef itk::Statistics::ListSample< CentroidVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType >    TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer               TreePointer;
  typedef typename TreeGeneratorType::KdTreeType            TreeType;
  typedef typename TreeType::Pointer                        KdTreePointer;


  void SetFunctionCount(const unsigned int & n)
  {
    this->m_FunctionCount = n;
  }

  void SetNumberOfNeighbors(const unsigned int & n)
  {
    this->m_NumberOfNeighbors = n;
  }

  void SetKdTree(KdTreePointer kdtree)
  {
    this->m_KdTree = kdtree;
  }

  void AllocateListImage(const FeatureImageType *featureImage)
  {
    this->m_NearestNeighborListImage = ListImageType::New();
    this->m_NearestNeighborListImage->CopyInformation(featureImage);
    this->m_NearestNeighborListImage->SetRegions( featureImage->GetLargestPossibleRegion() );
    this->m_NearestNeighborListImage->Allocate();
  }

  virtual void PopulateListImage() = 0;

  unsigned int     m_FunctionCount;
  unsigned int     m_NumberOfNeighbors;
  ListImagePointer m_NearestNeighborListImage;
  KdTreePointer    m_KdTree;
  bool             m_UsePartitionedDomain;

protected:
  LevelSetDomainPartitionBase()
  {
    m_NumberOfNeighbors = 6;
    m_NearestNeighborListImage = 0;
    m_KdTree = 0;
    m_UsePartitionedDomain = true;
  }

  ~LevelSetDomainPartitionBase(){}

private:
  LevelSetDomainPartitionBase(const Self &); //purposely not
                                                       // implemented
  void operator=(const Self &);                        //purposely not
                                                       // implemented
};
} //end namespace itk

#endif
