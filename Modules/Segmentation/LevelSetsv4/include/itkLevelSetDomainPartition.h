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
#ifndef __itkLevelSetDomainPartition_h
#define __itkLevelSetDomainPartition_h

#include "itkLevelSetDomainPartitionBase.h"

namespace itk
{
/** \class LevelSetDomainPartition
 *
 * \brief Helper class used to share data in the ScalarChanAndVeseLevelSetFunction.
 *
 */
template< class TInputImage, class TFeatureImage >
class LevelSetDomainPartition:
  public LevelSetDomainPartitionBase< TInputImage, TFeatureImage >
{
public:

  typedef LevelSetDomainPartition Self;
  typedef LevelSetDomainPartitionBase< TInputImage, TFeatureImage >
  Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, TFeatureImage::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkTypeMacro(LevelSetDomainPartition, LevelSetDomainPartitionBase );

  typedef TInputImage                                 InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename Superclass::InputPixelType         InputPixelType;
  typedef typename Superclass::InputRegionType        InputRegionType;
  typedef typename Superclass::InputSizeType          InputSizeType;
  typedef typename Superclass::InputSizeValueType     InputSizeValueType;
  typedef typename Superclass::InputSpacingType       InputSpacingType;
  typedef typename Superclass::InputIndexType         InputIndexType;
  typedef typename Superclass::InputIndexValueType    InputIndexValueType;
  typedef typename Superclass::InputPointType         InputPointType;

  typedef TFeatureImage                                 FeatureImageType;
  typedef typename Superclass::FeatureImagePointer      FeatureImagePointer;
  typedef typename Superclass::FeatureImageConstPointer FeatureImageConstPointer;
  typedef typename Superclass::FeaturePixelType         FeaturePixelType;
  typedef typename Superclass::FeatureRegionType        FeatureRegionType;
  typedef typename Superclass::FeatureSizeType          FeatureSizeType;
  typedef typename Superclass::FeatureSizeValueType     FeatureSizeValueType;
  typedef typename Superclass::FeatureSpacingType       FeatureSpacingType;
  typedef typename Superclass::FeatureIndexType         FeatureIndexType;
  typedef typename Superclass::FeaturePointType         FeaturePointType;

  typedef typename Superclass::ListPixelType         ListPixelType;
  typedef typename Superclass::ListImageType         ListImageType;
  typedef typename Superclass::ListImagePointer      ListImagePointer;
  typedef typename Superclass::ListImageConstPointer ListImageConstPointer;
  typedef typename Superclass::ListRegionType        ListRegionType;
  typedef typename Superclass::ListSizeType          ListSizeType;
  typedef typename Superclass::ListSizeValueType     ListSizeValueType;
  typedef typename Superclass::ListSpacingType       ListSpacingType;
  typedef typename Superclass::ListIndexType         ListIndexType;
  typedef typename Superclass::ListIndexValueType    ListIndexValueType;
  typedef typename Superclass::ListPointType         ListPointType;
  typedef typename Superclass::ListIteratorType      ListIteratorType;

  typedef typename Superclass::CentroidVectorType CentroidVectorType;
  typedef typename Superclass::SampleType         SampleType;
  typedef typename Superclass::TreeGeneratorType  TreeGeneratorType;
  typedef typename Superclass::TreePointer        TreePointer;
  typedef typename Superclass::TreeType           TreeType;
  typedef typename Superclass::KdTreePointer      KdTreePointer;

  void PopulateListImage()
  {
    ListPixelType L;

    for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
      {
      L.push_back(i);
      }
    this->m_NearestNeighborListImage->FillBuffer(L);
  }

protected:
  LevelSetDomainPartition():Superclass(){}
  ~LevelSetDomainPartition(){}
private:
  LevelSetDomainPartition(const Self &); //purposely
                                                                    // not
                                                                    // implemented
  void operator=(const Self &);                                     //purposely
                                                                    // not
                                                                    // implemented
};
} //end namespace itk

#endif
