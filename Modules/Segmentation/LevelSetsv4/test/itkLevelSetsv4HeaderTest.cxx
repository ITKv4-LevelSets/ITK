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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>

#include "itkSparseCityBlockNeighborList.h"
#include "itkLevelSetDomainPartitionBase.h"
#include "itkLevelSetDomainPartition.h"
#include "itkLevelSetDomainPartitionImageWithKdTree.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkLevelSetImageBase.h"
#include "itkLevelSetQuadEdgeMeshBase.h"
#include "itkLevelSetBase.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationTermBase.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"
#include "itkLevelSetEquationChanAndVeseExternalTerm.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEvolutionBase.h"
#include "itkWhitakerSparseLevelSetBase.h"
#include "itkShiSparseLevelSetBase.h"
#include "itkMalcolmSparseLevelSetBase.h"
// #include "itkLevelSetWhitakerSparseEvolutionBase.h"
#include "itkBinaryImageToWhitakerSparseLevelSetAdaptor.h"
#include "itkBinaryImageToMalcolmSparseLevelSetAdaptor.h"
#include "itkBinaryImageToShiSparseLevelSetAdaptor.h"
#include "itkUpdateWhitakerSparseLevelSet.h"
#include "itkUpdateShiSparseLevelSet.h"

int itkLevelSetsv4HeaderTest ( int , char * [] )
{
  return EXIT_SUCCESS;
}
