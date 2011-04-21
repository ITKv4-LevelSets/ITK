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

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEquationChanAndVeseExternalTerm.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"
#include "itkLevelSetEvolutionBase.h"

int itkMultiLevelSetEvolutionTest( int , char* [] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char                                       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >             InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > InputIteratorType;

  typedef float                                          PixelType;
  typedef itk::Image< PixelType, Dimension >             ImageType;
  typedef itk::LevelSetImageBase< ImageType >            LevelSetType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::IdentifierType                            IdentifierType;
  typedef std::list< IdentifierType >                    IdListType;
  typedef itk::Image< IdListType, Dimension >            IdListImageType;
  typedef itk::Image< short, Dimension >                 CacheImageType;
  typedef itk::LevelSetDomainMapImageFilter< IdListImageType, CacheImageType >
                                                         DomainMapImageFilterType;

  typedef itk::LevelSetContainerBase< IdentifierType, LevelSetType >  LevelSetContainerType;
  typedef itk::LevelSetEquationChanAndVeseInternalTerm< InputImageType, LevelSetContainerType >
                                                                      ChanAndVeseInternalTermType;
  typedef itk::LevelSetEquationChanAndVeseExternalTerm< InputImageType, LevelSetContainerType >
                                                                      ChanAndVeseExternalTermType;
  typedef itk::LevelSetEquationTermContainerBase< InputImageType, LevelSetContainerType >
                                                                      TermContainerType;

  typedef itk::LevelSetEquationContainerBase< TermContainerType >     EquationContainerType;

  typedef itk::LevelSetEvolutionBase< EquationContainerType >             LevelSetEvolutionType;
  typedef itk::AtanRegularizedHeavisideStepFunction< PixelType, PixelType >
                                                                      HeavisideFunctionBaseType;

  ImageType::IndexType index;
  index[0] = 0;
  index[1] = 0;

  ImageType::SizeType size;
  size[0] = 10;
  size[1] = 10;

  ImageType::RegionType region;
  region.SetIndex( index );
  region.SetSize( size );

  PixelType value = 0.;

  InputImageType::Pointer input = InputImageType::New();
  input->SetRegions( region );
  input->Allocate();
  input->FillBuffer( 1 );

//  InputIteratorType it( input, input->GetLargestPossibleRegion() );
//  it.GoToBegin();
//  while( !it.IsAtEnd() )
//  {
//    it.Set( 1 );
//    ++it;
//  }

  ImageType::Pointer input1 = ImageType::New();
  input1->SetRegions( region );
  input1->Allocate();
  input1->FillBuffer( value );

  ImageType::Pointer input2 = ImageType::New();
  input2->SetRegions( region );
  input2->Allocate();
  input2->FillBuffer( value );

  ImageType::IndexType idx;
  IdListType list_ids;

  IdListImageType::Pointer id_image = IdListImageType::New();
  id_image->SetRegions( region );
  id_image->Allocate();
  id_image->FillBuffer( list_ids );

  IteratorType it1( input1, input1->GetLargestPossibleRegion() );
  IteratorType it2( input2, input2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();

  while( !it1.IsAtEnd() )
    {
    idx = it1.GetIndex();
    list_ids.clear();

    if( ( idx[0] < 5 ) && ( idx[1] < 5 ) )
      {
      list_ids.push_back( 1 );
      }

    if( ( idx[0] > 1 ) && ( idx[1] > 1 ) &&
        ( idx[0] < 8 ) && ( idx[1] < 8 ) )
      {
      list_ids.push_back( 2 );
      }

    id_image->SetPixel( idx, list_ids );

    it1.Set( vcl_sqrt(
             static_cast< float> ( ( idx[0] - 2 ) * ( idx[0] - 2 ) +
                                   ( idx[1] - 2 ) * ( idx[1] - 2 ) ) ) );

    it2.Set( vcl_sqrt(
             static_cast< float> ( ( idx[0] - 5 ) * ( idx[0] - 5 ) +
                                   ( idx[1] - 5 ) * ( idx[1] - 5 ) ) ) );
    ++it1;
    ++it2;
    }

  // Map of levelset bases
  std::map< itk::IdentifierType, LevelSetType::Pointer > level_set;
  level_set[1] = LevelSetType::New();
  level_set[1]->SetImage( input1 );

  level_set[2] = LevelSetType::New();
  level_set[2]->SetImage( input2 );

  // Insert the levelsets in a levelset container
  LevelSetContainerType::Pointer lscontainer = LevelSetContainerType::New();
  bool LevelSetNotYetAdded = lscontainer->AddLevelSet( 0, level_set[1], false );

  if ( !LevelSetNotYetAdded )
    {
    return EXIT_FAILURE;
    }

  LevelSetNotYetAdded = lscontainer->AddLevelSet( 1, level_set[2], false );
  if ( !LevelSetNotYetAdded )
    {
    return EXIT_FAILURE;
    }
  std::cout << "Level set container created" << std::endl;

  // Define the Heaviside function
  HeavisideFunctionBaseType::Pointer heaviside = HeavisideFunctionBaseType::New();
  heaviside->SetEpsilon( 1.0 );

  // **************** CREATE ALL TERMS ****************

  // -----------------------------
  // *** 1st Level Set phi_{1} ***

  // Create ChanAndVese internal term for phi_{1}
  ChanAndVeseInternalTermType::Pointer cvInternalTerm0 = ChanAndVeseInternalTermType::New();
  cvInternalTerm0->SetHeaviside( heaviside );
  cvInternalTerm0->SetInput( input );
  cvInternalTerm0->SetCoefficient( 1.0 );
  cvInternalTerm0->SetCurrentLevelSet( 0 );
  cvInternalTerm0->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 1: CV internal term created" << std::endl;

  // Create ChanAndVese external term for phi_{1}
  ChanAndVeseExternalTermType::Pointer cvExternalTerm0 = ChanAndVeseExternalTermType::New();
  cvExternalTerm0->SetHeaviside( heaviside );
  cvExternalTerm0->SetInput( input );
  cvExternalTerm0->SetCoefficient( 1.0 );
  cvExternalTerm0->SetCurrentLevelSet( 0 );
  cvExternalTerm0->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 1: CV external term created" << std::endl;

  // -----------------------------
  // *** 2nd Level Set phi_{2} ***

  // Create ChanAndVese internal term for phi_{1}
  ChanAndVeseInternalTermType::Pointer cvInternalTerm1 = ChanAndVeseInternalTermType::New();
  cvInternalTerm1->SetHeaviside( heaviside );
  cvInternalTerm1->SetInput( input );
  cvInternalTerm1->SetCoefficient( 1.0 );
  cvInternalTerm1->SetCurrentLevelSet( 1 );
  cvInternalTerm1->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 2: CV internal term created" << std::endl;

  // Create ChanAndVese external term for phi_{2}
  ChanAndVeseExternalTermType::Pointer cvExternalTerm1 = ChanAndVeseExternalTermType::New();
  cvExternalTerm1->SetHeaviside( heaviside );
  cvExternalTerm1->SetInput( input );
  cvExternalTerm1->SetCoefficient( 1.0 );
  cvExternalTerm1->SetCurrentLevelSet( 1 );
  cvExternalTerm1->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 2: CV external term created" << std::endl;

  // **************** CREATE ALL EQUATIONS ****************

  // Create Term Container
  TermContainerType::Pointer termContainer0 = TermContainerType::New();
  termContainer0->SetInput( input );

  TermContainerType::TermPointer temp;
  temp = dynamic_cast< TermContainerType::TermType* >( cvInternalTerm0.GetPointer() );
  termContainer0->AddTerm( 0, temp );

  temp = dynamic_cast< TermContainerType::TermType* >( cvExternalTerm0.GetPointer() );
  termContainer0->AddTerm( 1, temp );
  std::cout << "Term container 0 created" << std::endl;

  EquationContainerType::Pointer equationContainer = EquationContainerType::New();
  equationContainer->AddEquation( 0, termContainer0 );

  TermContainerType::Pointer termContainer1 = TermContainerType::New();
  termContainer1->SetInput( input );

  temp = dynamic_cast< TermContainerType::TermType* >( cvInternalTerm1.GetPointer() );
  termContainer1->AddTerm( 0, temp );

  temp = dynamic_cast< TermContainerType::TermType* >( cvExternalTerm1.GetPointer() );
  termContainer1->AddTerm( 1, temp );
  std::cout << "Term container 1 created" << std::endl;

  equationContainer->AddEquation( 1, termContainer1 );

  DomainMapImageFilterType::Pointer domainMapFilter = DomainMapImageFilterType::New();
  domainMapFilter->SetInput( id_image );
  domainMapFilter->Update();

  LevelSetEvolutionType::Pointer evolution = LevelSetEvolutionType::New();
  evolution->SetEquationContainer( equationContainer );
  evolution->SetNumberOfIterations( 2 );
  evolution->SetLevelSetContainer( lscontainer );
  evolution->SetDomainMapFilter( domainMapFilter );
  evolution->Update();

  return EXIT_SUCCESS;
}
