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
#ifndef __itkLevelSetDomainPartitionImageBase_h
#define __itkLevelSetDomainPartitionImageBase_h

#include "itkLevelSetDomainPartitionBase.h"

#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
/** \class LevelSetDomainPartitionImageBase
 *
 * \brief Helper class used to partition domain and efficiently compute overlap.
 *
 */
template< class TImage >
class LevelSetDomainPartitionImageBase : 
  public LevelSetDomainPartitionBase 
{
public:

  typedef LevelSetDomainPartitionImageBase      Self;
  typedef LevelSetDomainPartitionBase           Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  itkTypeMacro( LevelSetDomainPartitionImageBase, 
                LevelSetDomainPartitionBase );

  typedef TImage                             ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::SpacingType    SpacingType;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename IndexType::IndexValueType IndexValueType;
  typedef typename ImageType::PointType      PointType;

  typedef typename Superclass::ListPixelType ListPixelType;

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

  void SetImage( ImagePointer iImage )
  {
    m_Image = iImage;
  }
/*
  typedef Vector< float, itkGetStaticConstMacro(ImageDimension) >
    CentroidVectorType;
  typedef itk::Statistics::ListSample< CentroidVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType >    TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer               TreePointer;
  typedef typename TreeGeneratorType::KdTreeType            TreeType;
  typedef typename TreeType::Pointer                        KdTreePointer;

  void SetKdTree(KdTreePointer kdtree)
  {
    this->m_KdTree = kdtree;
  }

  unsigned int     m_NumberOfNeighbors;
  KdTreePointer    m_KdTree;
  bool             m_UsePartitionedDomain;
*/

protected:
  LevelSetDomainPartitionImageBase()
  {
//    m_NearestNeighborListImage = 0;
//    m_KdTree = 0;
//    m_UsePartitionedDomain = true;
  }

  virtual ~LevelSetDomainPartitionImageBase(){}

  ImagePointer     m_Image;
  ListImagePointer m_NearestNeighborListImage;


  virtual void PopulateListDomain()
  {
    ListSpacingType spacing = this->m_NearestNeighborListImage->GetSpacing();

    ListRegionType region = this->m_NearestNeighborListImage->GetLargestPossibleRegion();

    ListIteratorType lIt(this->m_NearestNeighborListImage, region);

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

  void AllocateListDomain()
  {
    this->m_NearestNeighborListImage = ListImageType::New();
    this->m_NearestNeighborListImage->CopyInformation( m_Image );
    this->m_NearestNeighborListImage->SetRegions( m_Image->GetLargestPossibleRegion() );
    this->m_NearestNeighborListImage->Allocate();
  }

private:
  LevelSetDomainPartitionImageBase(const Self &); //purposely not
                                                       // implemented
  void operator=(const Self &);                        //purposely not
                                                       // implemented
};
} //end namespace itk

#endif
