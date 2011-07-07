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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEquationChanAndVeseExternalTerm.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"
#include "itkSinRegularizedHeavisideStepFunction.h"
#include "itkLevelSetSparseEvolutionBase.h"
#include "itkBinaryImageToWhitakerSparseLevelSetAdaptor.h"
#include "itkNumericTraits.h"

int itkTwoLevelSetSparse2DTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned short                                    InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >           InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType >
                                                            InputIteratorType;
  typedef itk::ImageFileReader< InputImageType >            ReaderType;

  typedef float                                             PixelType;

  typedef itk::BinaryImageToWhitakerSparseLevelSetAdaptor< InputImageType, PixelType >
                                                            BinaryToSparseAdaptorType;

  typedef itk::IdentifierType                               IdentifierType;
  typedef BinaryToSparseAdaptorType::LevelSetType           SparseLevelSetType;
  typedef SparseLevelSetType::ImageType                     SparseImageType;
  typedef SparseLevelSetType::NodeAttributeType             NodeAttributeType;

  typedef itk::LevelSetContainerBase< IdentifierType, SparseLevelSetType >
                                                            LevelSetContainerType;

  typedef std::list< IdentifierType >                       IdListType;
  typedef itk::Image< IdListType, Dimension >               IdListImageType;
  typedef itk::Image< short, Dimension >                    CacheImageType;
  typedef itk::LevelSetDomainMapImageFilter< IdListImageType, CacheImageType >
                                                            DomainMapImageFilterType;

  typedef itk::LevelSetEquationChanAndVeseInternalTerm< InputImageType, LevelSetContainerType >
                                                            ChanAndVeseInternalTermType;
  typedef itk::LevelSetEquationChanAndVeseExternalTerm< InputImageType, LevelSetContainerType >
                                                            ChanAndVeseExternalTermType;
  typedef itk::LevelSetEquationTermContainerBase< InputImageType, LevelSetContainerType >
                                                            TermContainerType;

  typedef itk::LevelSetEquationContainerBase< TermContainerType >
                                                            EquationContainerType;

  typedef itk::LevelSetSparseEvolutionBase< EquationContainerType >
                                                            LevelSetEvolutionType;

  typedef SparseLevelSetType::OutputRealType                      LevelSetOutputRealType;
  typedef itk::SinRegularizedHeavisideStepFunction< LevelSetOutputRealType, LevelSetOutputRealType >
                                                            HeavisideFunctionBaseType;
  typedef itk::ImageRegionIteratorWithIndex< SparseImageType >    IteratorType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType >     InputIteratorType;

  // load binary input for segmentation
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  InputImageType::Pointer input = reader->GetOutput();

  // Create a binary initialization
  InputImageType::Pointer binary = InputImageType::New();
  binary->SetRegions( input->GetLargestPossibleRegion() );
  binary->CopyInformation( input );
  binary->Allocate();
  binary->FillBuffer( itk::NumericTraits<InputPixelType>::Zero );

  InputImageType::RegionType region;
  InputImageType::IndexType index;
  InputImageType::SizeType size;

  index.Fill( 10 );
  size.Fill( 30 );

  region.SetIndex( index );
  region.SetSize( size );

  InputIteratorType iIt( binary, region );
  iIt.GoToBegin();
  while( !iIt.IsAtEnd() )
  {
    iIt.Set( itk::NumericTraits<InputPixelType>::One );
    ++iIt;
  }

  // Convert binary mask to sparse level set
  BinaryToSparseAdaptorType::Pointer adaptor0 = BinaryToSparseAdaptorType::New();
  adaptor0->SetInputImage( binary );
  adaptor0->Initialize();
  std::cout << "Finished converting to sparse format" << std::endl;

  SparseLevelSetType::Pointer level_set0 = adaptor0->GetSparseLevelSet();
  SparseImageType::Pointer sparseImage0 = level_set0->GetImage();

  BinaryToSparseAdaptorType::Pointer adaptor1 = BinaryToSparseAdaptorType::New();
  adaptor1->SetInputImage( binary );
  adaptor1->Initialize();
  std::cout << "Finished converting to sparse format" << std::endl;

  SparseLevelSetType::Pointer level_set1 = adaptor1->GetSparseLevelSet();
  SparseImageType::Pointer sparseImage1 = level_set1->GetImage();

  // Create a list image specifying both level set ids
  IdListType list_ids;
  list_ids.push_back( 1 );
  list_ids.push_back( 2 );

  IdListImageType::Pointer id_image = IdListImageType::New();
  id_image->SetRegions( input->GetLargestPossibleRegion() );
  id_image->Allocate();
  id_image->FillBuffer( list_ids );

  DomainMapImageFilterType::Pointer domainMapFilter = DomainMapImageFilterType::New();
  domainMapFilter->SetInput( id_image );
  domainMapFilter->Update();
  std::cout << "Domain map computed" << std::endl;

  // Define the Heaviside function
  HeavisideFunctionBaseType::Pointer heaviside = HeavisideFunctionBaseType::New();
  heaviside->SetEpsilon( 1.0 );

  // Insert the levelsets in a levelset container
  LevelSetContainerType::Pointer lscontainer = LevelSetContainerType::New();
  lscontainer->SetHeaviside( heaviside );
  lscontainer->SetDomainMapFilter( domainMapFilter );

  bool LevelSetNotYetAdded = lscontainer->AddLevelSet( 0, level_set0, false );
  if ( !LevelSetNotYetAdded )
    {
    return EXIT_FAILURE;
    }

  LevelSetNotYetAdded = lscontainer->AddLevelSet( 1, level_set1, false );
  if ( !LevelSetNotYetAdded )
  {
    return EXIT_FAILURE;
  }
  std::cout << "Level set container created" << std::endl;

  // **************** CREATE ALL TERMS ****************

  // -----------------------------
  // *** 1st Level Set phi ***

  // Create ChanAndVese internal term for phi_{1}
  ChanAndVeseInternalTermType::Pointer cvInternalTerm0 = ChanAndVeseInternalTermType::New();
  cvInternalTerm0->SetInput( input );
  cvInternalTerm0->SetCoefficient( 1.0 );
  cvInternalTerm0->SetCurrentLevelSet( 0 );
  cvInternalTerm0->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 1: CV internal term created" << std::endl;

  // Create ChanAndVese external term for phi_{1}
  ChanAndVeseExternalTermType::Pointer cvExternalTerm0 = ChanAndVeseExternalTermType::New();
  cvExternalTerm0->SetInput( input );
  cvExternalTerm0->SetCoefficient( 1.0 );
  cvExternalTerm0->SetCurrentLevelSet( 0 );
  cvExternalTerm0->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 1: CV external term created" << std::endl;

  // -----------------------------
  // *** 2nd Level Set phi ***
  ChanAndVeseInternalTermType::Pointer cvInternalTerm1 = ChanAndVeseInternalTermType::New();
  cvInternalTerm1->SetInput( input );
  cvInternalTerm1->SetCoefficient( 1.0 );
  cvInternalTerm1->SetCurrentLevelSet( 1 );
  cvInternalTerm1->SetLevelSetContainer( lscontainer );
  std::cout << "LevelSet 2: CV internal term created" << std::endl;

  // Create ChanAndVese external term for phi_{1}
  ChanAndVeseExternalTermType::Pointer cvExternalTerm1 = ChanAndVeseExternalTermType::New();
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

  // Create Term Container
  TermContainerType::Pointer termContainer1 = TermContainerType::New();
  termContainer1->SetInput( input );

  temp = dynamic_cast< TermContainerType::TermType* >( cvInternalTerm1.GetPointer() );
  termContainer1->AddTerm( 0, temp );

  temp = dynamic_cast< TermContainerType::TermType* >( cvExternalTerm1.GetPointer() );
  termContainer1->AddTerm( 1, temp );
  std::cout << "Term container 1 created" << std::endl;

  // Create equation container
  EquationContainerType::Pointer equationContainer = EquationContainerType::New();
  equationContainer->AddEquation( 0, termContainer0 );
  equationContainer->AddEquation( 1, termContainer1 );

  LevelSetEvolutionType::Pointer evolution = LevelSetEvolutionType::New();
  evolution->SetEquationContainer( equationContainer );
  evolution->SetNumberOfIterations( 40 );
  evolution->SetLevelSetContainer( lscontainer );
  evolution->SetDomainMapFilter( domainMapFilter );

  try
    {
    evolution->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }

  PixelType internalmean1 = cvInternalTerm0->GetMean();
  PixelType internalmean2 = cvInternalTerm1->GetMean();
  if ( ( internalmean1 < 42100 ) || ( internalmean1 > 42150 ) )
  {
    std::cout << "( ( mean1 < 42100 ) || ( mean1 > 42150 ) )" <<std::endl;
    std::cout << "internalmean1 = " << internalmean1 <<std::endl;
    return EXIT_FAILURE;
  }

  PixelType externalmean1 = cvExternalTerm0->GetMean();
  PixelType externalmean2 = cvExternalTerm1->GetMean();
  if ( ( externalmean1 < 1500 ) || ( externalmean1 > 1550 ) )
  {
    std::cout << "( ( externalmean1 < 1500 ) || ( externalmean1 > 1550 ) )" <<std::endl;
    std::cout << "externalmean1 = " << externalmean1 <<std::endl;
    return EXIT_FAILURE;
  }


  if ( ( internalmean1 != internalmean2  ) || ( externalmean1 != externalmean2 ) )
  {
    std::cout << "internalmean = " << internalmean1 <<std::endl;
    std::cout << "externalmean = " << externalmean1 <<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
