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

#include "itkMalcolmSparseLevelSetBase.h"

int itkMalcolmSparseLevelSetBaseTest( int , char* [] )
{
  const unsigned int Dimension = 2;
  typedef itk::MalcolmSparseLevelSetBase< Dimension > SparseLevelSetType;

  SparseLevelSetType::Pointer phi = SparseLevelSetType::New();

  typedef SparseLevelSetType::SparseImageType SparseImageType;

  SparseImageType::IndexType index;
  index[0] = 0;
  index[1] = 0;

  SparseImageType::SizeType size;
  size[0] = 10;
  size[1] = 10;

  SparseImageType::RegionType region;
  region.SetIndex( index );
  region.SetSize( size );

  SparseLevelSetType::OutputType value = -1;

  SparseImageType::Pointer image = SparseImageType::New();
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer( value );

  phi->SetImage( image );

  return EXIT_SUCCESS;
}
