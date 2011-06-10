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

int itkSingleLevelSetSparse2DTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char                                     InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >           InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType >
                                                            InputIteratorType;
  typedef itk::ImageFileReader< InputImageType >            ReaderType;

  typedef float                                             OutputPixelType;

  typedef itk::BinaryImageToWhitakerSparseLevelSetAdaptor< InputImageType, OutputPixelType >
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
  typedef itk::SinRegularizedHeavisideStepFunction< OutputPixelType, OutputPixelType >
                                                            HeavisideFunctionBaseType;
  typedef itk::ImageRegionIteratorWithIndex< SparseImageType >    IteratorType;

  // load binary mask
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  InputImageType::Pointer input = reader->GetOutput();

  // Convert binary mask to sparse level set
  BinaryToSparseAdaptorType::Pointer adaptor = BinaryToSparseAdaptorType::New();
  adaptor->SetInputImage( input );
  adaptor->Initialize();
  std::cout << "Finished converting to sparse format" << std::endl;

  SparseLevelSetType::Pointer level_set = adaptor->GetSparseLevelSet();
  SparseImageType::Pointer sparseImage = level_set->GetImage();

  IdListType list_ids;
  list_ids.push_back( 1 );

  IdListImageType::Pointer id_image = IdListImageType::New();
  id_image->SetRegions( input->GetLargestPossibleRegion() );
  id_image->Allocate();
  id_image->FillBuffer( list_ids );

  // Insert the levelsets in a levelset container
  LevelSetContainerType::Pointer lscontainer = LevelSetContainerType::New();
  bool LevelSetNotYetAdded = lscontainer->AddLevelSet( 0, level_set, false );

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
  // *** 1st Level Set phi ***

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

  DomainMapImageFilterType::Pointer domainMapFilter = DomainMapImageFilterType::New();
  domainMapFilter->SetInput( id_image );
  domainMapFilter->Update();

  LevelSetEvolutionType::Pointer evolution = LevelSetEvolutionType::New();
  evolution->SetEquationContainer( equationContainer );
  evolution->SetNumberOfIterations( 2 );
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

  PixelType mean = cvInternalTerm0->GetMean();
  if ( ( mean < 24900 ) || ( mean > 24910 ) )
  {
    return EXIT_FAILURE;
  }

  mean = cvExternalTerm0->GetMean();
  if ( ( mean < 1350 ) || ( mean > 1360 ) )
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
