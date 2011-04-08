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
#ifndef __itkLevelSetDomainPartitionImageWithKdTree_h
#define __itkLevelSetDomainPartitionImageWithKdTree_h

#include "itkLevelSetDomainPartitionBase.h"

namespace itk
{
/** \class LevelSetDomainPartitionImageWithKdTree
 *
 * \brief Helper class used to share data in the ScalarChanAndVeseLevelSetFunction.
 *
 */
template< class TImage >
class LevelSetDomainPartitionImageWithKdTree:
  public LevelSetDomainPartitionImageBase< TImage >
{
public:

  typedef LevelSetDomainPartitionImageWithKdTree Self;
  typedef LevelSetDomainPartitionBase< TImage >
    Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkTypeMacro( LevelSetDomainPartitionImageWithKdTree, 
                LevelSetDomainPartitionImageBase);

  typedef TImage                                   InputImageType;
  typedef typename Superclass::InputImagePointer        InputImagePointer;
  typedef typename Superclass::InputImageConstPointer   InputImageConstPointer;
  typedef typename Superclass::InputPixelType           InputPixelType;
  typedef typename Superclass::InputRegionType          InputRegionType;
  typedef typename Superclass::InputSizeType            InputSizeType;
  typedef typename Superclass::InputSizeValueType       InputSizeValueType;
  typedef typename Superclass::InputSpacingType         InputSpacingType;
  typedef typename Superclass::InputIndexType           InputIndexType;
  typedef typename Superclass::InputIndexValueType      InputIndexValueType;
  typedef typename Superclass::InputPointType           InputPointType;

  typedef typename Superclass::ListPixelType            ListPixelType;
  typedef typename Superclass::ListImageType            ListImageType;
  typedef typename Superclass::ListImagePointer         ListImagePointer;
  typedef typename Superclass::ListImageConstPointer    ListImageConstPointer;
  typedef typename Superclass::ListRegionType           ListRegionType;
  typedef typename Superclass::ListSizeType             ListSizeType;
  typedef typename Superclass::ListSizeValueType        ListSizeValueType;
  typedef typename Superclass::ListSpacingType          ListSpacingType;
  typedef typename Superclass::ListIndexType            ListIndexType;
  typedef typename Superclass::ListIndexValueType       ListIndexValueType;
  typedef typename Superclass::ListPointType            ListPointType;
  typedef typename Superclass::ListIteratorType         ListIteratorType;

  typedef typename Superclass::CentroidVectorType CentroidVectorType;
  typedef typename Superclass::SampleType         SampleType;
  typedef typename Superclass::TreeGeneratorType  TreeGeneratorType;
  typedef typename Superclass::TreePointer        TreePointer;
  typedef typename Superclass::TreeType           TreeType;
  typedef typename Superclass::KdTreePointer      KdTreePointer;


  void PopulateListImage()
  {
    ListSpacingType spacing = this->m_NearestNeighborListImage->GetSpacing();

    ListRegionType region = this->m_NearestNeighborListImage->GetLargestPossibleRegion();

    ListIteratorType lIt(this->m_NearestNeighborListImage, region);

    if ( this->m_KdTree.IsNotNull() )
      {
      for ( lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
        {
        ListIndexType ind = lIt.GetIndex();

        float queryPoint[ImageDimension];
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          queryPoint[i] = ind[i] * spacing[i];
          }

        typename TreeType::InstanceIdentifierVectorType neighbors;
        this->m_KdTree->Search(queryPoint, this->m_NumberOfNeighbors, neighbors);

        ListPixelType L;
        for ( unsigned int i = 0; i < this->m_NumberOfNeighbors; i++ )
          {
          if ( this->m_LevelSetDataPointerVector[i]->VerifyInsideRegion(ind) )
            {
            L.push_back(neighbors[i]);
            }
          }
        lIt.Set(L);
        }
      }
    else
      {
      for ( lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
        {
        ListIndexType ind = lIt.GetIndex();
        ListPixelType L;
        for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
          {
          if ( this->m_LevelSetDataPointerVector[i]->VerifyInsideRegion(ind) )
            {
            L.push_back(i);
            }
          }
        lIt.Set(L);
        }
      }
  }

protected:
  LevelSetDomainPartitionImageWithKdTree():Superclass(){}
  ~LevelSetDomainPartitionImageWithKdTree(){}

private:
  //purposely not implemented
  LevelSetDomainPartitionImageWithKdTree(const Self &);
  //purposely not implemented
  void operator=(const Self &);
};
} //end namespace itk

#endif
