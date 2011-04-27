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

#include "itkLevelSetDomainPartitionBase.h"
#include "itkImage.h"
#include "itkTestingMacros.h"

namespace itk
{

template < class TDomain >
class LevelSetDomainPartitionBaseHelper
  : public LevelSetDomainPartitionBase< TDomain >
{
public:
  /** Standard class typedefs. */
  typedef LevelSetDomainPartitionBaseHelper      Self;
  typedef LevelSetDomainPartitionBase< TDomain > Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( LevelSetDomainPartitionBaseHelper, LevelSetDomainPartitionBase );

  itkNewMacro( Self );

protected:
  void AllocateListDomain() {}
  void PopulateListDomain() {}
};

}


int itkLevelSetDomainPartitionBaseTest( int, char* [] )
{
  const unsigned int Dimension = 3;

  typedef itk::Image< double, Dimension >         ImageType;

  typedef itk::LevelSetDomainPartitionBaseHelper< ImageType >
    DomainPartitionBaseHelperType;

  unsigned int count = 2;

  DomainPartitionBaseHelperType::Pointer function = DomainPartitionBaseHelperType::New();
  function->SetFunctionCount( count );

  if( function->GetFunctionCount() != count )
    {
    return EXIT_FAILURE;
    }

//   std::cout << "GetNameOfClass() = " << function->GetNameOfClass() << std::endl;
//   function->Print( std::cout );

  return EXIT_SUCCESS;
}