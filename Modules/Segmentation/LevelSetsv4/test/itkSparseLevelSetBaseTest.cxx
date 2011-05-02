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

#include "itkSparseLevelSetBase.h"

namespace itk
{
template< typename TOutput, unsigned int VDimension >
class SparseLevelSetBaseTestHelper :
    public SparseLevelSetBase< TOutput, VDimension >
{
public:
  typedef SparseLevelSetBaseTestHelper  Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef SparseLevelSetBase< TOutput, VDimension >
                                        Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SparseLevelSetBaseTestHelper, SparseLevelSetBase);

protected:
  SparseLevelSetBaseTestHelper() : Superclass() {}
  ~SparseLevelSetBaseTestHelper() {}

  void InitializeLayers() {}

private:
  SparseLevelSetBaseTestHelper( const Self& );
  void operator = ( const Self& );
};
  }

int itkSparseLevelSetBaseTest( int , char* [] )
{
  typedef double OutputType;
  const unsigned int Dimension = 3;
  typedef itk::SparseLevelSetBaseTestHelper< OutputType, Dimension >
      SparseLevelSetType;

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

  SparseLevelSetType::NodeAttributeType value;
  value.m_Status = -3;
  value.m_Value = -3.;

  SparseImageType::Pointer image = SparseImageType::New();
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer( value );

  phi->SetImage( image );

  return EXIT_SUCCESS;
}
