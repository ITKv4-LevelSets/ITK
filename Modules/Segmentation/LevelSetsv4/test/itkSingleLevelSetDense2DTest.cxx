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
#include "itkFastMarchingImageFilter.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkDenseLevelSetContainer.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEquationChanAndVeseExternalTerm.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"
#include "itkLevelSetEvolutionBase.h"

int itkSingleLevelSetDense2DTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned short                                      InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >             InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > InputIteratorType;
  typedef itk::ImageFileReader< InputImageType >              ReaderType;

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

  typedef itk::DenseLevelSetContainer< IdentifierType, LevelSetType >  LevelSetContainerType;
  typedef itk::LevelSetEquationChanAndVeseInternalTerm< InputImageType, LevelSetContainerType >
                                                                      ChanAndVeseInternalTermType;
  typedef itk::LevelSetEquationChanAndVeseExternalTerm< InputImageType, LevelSetContainerType >
                                                                      ChanAndVeseExternalTermType;
  typedef itk::LevelSetEquationTermContainerBase< InputImageType, LevelSetContainerType >
                                                                      TermContainerType;

  typedef itk::LevelSetEquationContainerBase< TermContainerType >     EquationContainerType;

  typedef itk::LevelSetEvolutionBase< EquationContainerType >         LevelSetEvolutionType;
  typedef itk::AtanRegularizedHeavisideStepFunction< PixelType, PixelType >
                                                                      HeavisideFunctionBaseType;

  typedef  itk::FastMarchingImageFilter< ImageType, ImageType > FastMarchingFilterType;
  typedef FastMarchingFilterType::NodeContainer                 NodeContainer;
  typedef FastMarchingFilterType::NodeType                      NodeType;

  // Read the image to be segmented
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  InputImageType::Pointer input = reader->GetOutput();

  FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();

  NodeContainer::Pointer seeds = NodeContainer::New();

  ImageType::IndexType  seedPosition;
  seedPosition[0] = atoi( argv[2] );
  seedPosition[1] = atoi( argv[3] );

  const double initialDistance = atof( argv[4] );
  const double seedValue = - initialDistance;

  NodeType node;
  node.SetValue( seedValue );
  node.SetIndex( seedPosition );

  //  The list of nodes is initialized and then every node is inserted using
  //  the \code{InsertElement()}.
  //
  seeds->Initialize();
  seeds->InsertElement( 0, node );

  //  The set of seed nodes is passed now to the
  //  FastMarchingImageFilter with the method
  //  \code{SetTrialPoints()}.
  //
  fastMarching->SetTrialPoints(  seeds  );

  //  Since the FastMarchingImageFilter is used here just as a
  //  Distance Map generator. It does not require a speed image as input.
  //  Instead the constant value $1.0$ is passed using the
  //  \code{SetSpeedConstant()} method.
  //
  fastMarching->SetSpeedConstant( 1.0 );

  //  The FastMarchingImageFilter requires the user to specify the
  //  size of the image to be produced as output. This is done using the
  //  \code{SetOutputSize()}. Note that the size is obtained here from the
  //  output image of the smoothing filter. The size of this image is valid
  //  only after the \code{Update()} methods of this filter has been called
  //  directly or indirectly.
  //
  fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
  fastMarching->Update();

  IdListType list_ids;
  list_ids.push_back( 1 );

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

  // Map of levelset bases
  LevelSetType::Pointer  level_set = LevelSetType::New();
  level_set->SetImage( fastMarching->GetOutput() );

  // Insert the levelsets in a levelset container
  LevelSetContainerType::Pointer lscontainer = LevelSetContainerType::New();
  lscontainer->SetHeaviside( heaviside );
  lscontainer->SetDomainMapFilter( domainMapFilter );

  bool LevelSetNotYetAdded = lscontainer->AddLevelSet( 0, level_set, false );
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

  LevelSetEvolutionType::Pointer evolution = LevelSetEvolutionType::New();
  evolution->SetEquationContainer( equationContainer );
  evolution->SetNumberOfIterations( 5 );
  evolution->SetLevelSetContainer( lscontainer );

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
