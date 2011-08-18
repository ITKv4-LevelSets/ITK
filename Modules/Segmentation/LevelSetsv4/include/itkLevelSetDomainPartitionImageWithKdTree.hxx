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
#ifndef __itkLevelSetDomainPartitionImageWithKdTree_hxx
#define __itkLevelSetDomainPartitionImageWithKdTree_hxx

#include "itkLevelSetDomainPartitionImageWithKdTree.h"

namespace itk
{
template< class TImage >
void LevelSetDomainPartitionImageWithKdTree< TImage >
::PopulateListDomain()
  {
    if( this->m_KdTree.IsNotNull() )
      {
      this->PopulateDomainWithKdTree();
      }
    else
      {
      Superclass::PopulateListDomain();
      }
  }

template< class TImage >
void LevelSetDomainPartitionImageWithKdTree< TImage >
::PopulateDomainWithKdTree()
  {
    ListSpacingType spacing = this->m_ListDomain->GetSpacing();

    ListRegionType region = this->m_ListDomain->GetLargestPossibleRegion();

    ListIteratorType lIt(this->m_ListDomain, region);

    for ( lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
      {
      ListIndexType ind = lIt.GetIndex();
      ListPointType pt;

      this->m_ListDomain->TransformIndexToPhysicalPoint( ind, pt );

      CentroidVectorType queryPoint = pt.GetVectorFromOrigin();

      typename TreeType::InstanceIdentifierVectorType neighbors;
      this->m_KdTree->Search(queryPoint, this->m_NumberOfNeighbors, neighbors);

      IdentifierListType L;
      for ( unsigned int i = 0; i < this->m_NumberOfNeighbors; ++i )
        {
        // this is not yet defined, but it will have to be !!!
        if ( this->m_LevelSetDataPointerVector[i]->VerifyInsideRegion(ind) )
          {
          L.push_back(neighbors[i]);
          }
        }
      lIt.Set(L);
      }
  }

} //end namespace itk
#endif
