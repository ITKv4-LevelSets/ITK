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

#include <string>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryImageToWhitakerSparseLevelSetAdaptor.h"
#include "itkUpdateWhitakerSparseLevelSet.h"
#include "itkImageRegionIterator.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"

int itkUpdateWhitakerSparseLevelSetTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char InputPixelType;
  typedef double        OutputPixelType;

  typedef itk::Image< InputPixelType, Dimension >   InputImageType;

  typedef itk::ImageFileReader< InputImageType >  InputReaderType;

  typedef itk::IdentifierType                         IdentifierType;

  InputImageType::Pointer input = InputImageType::New();
  InputImageType::RegionType region;

  InputImageType::SizeType size;
  size.Fill( 25 );

  InputImageType::IndexType start;
  start.Fill( 0 );

  region.SetSize( size );
  region.SetIndex( start );
  input->SetRegions( region );
  input->Allocate();
  input->FillBuffer( 0 );

  size[0] = 10;
  region.SetSize( size );

  typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;
  InputIteratorType i_it( input, region );
  i_it.GoToBegin();

  while( !i_it.IsAtEnd() )
    {
    i_it.Set( 255 );
    ++i_it;
    }

  std::cout << "Input image computed" << std::endl;

  typedef itk::BinaryImageToWhitakerSparseLevelSetAdaptor< InputImageType,
      OutputPixelType > BinaryToSparseAdaptorType;

  BinaryToSparseAdaptorType::Pointer adaptor = BinaryToSparseAdaptorType::New();
  adaptor->SetInputImage( input );
  adaptor->Initialize();

  typedef BinaryToSparseAdaptorType::LevelSetType     SparseLevelSetType;
  SparseLevelSetType::Pointer sparseLevelSet = adaptor->GetSparseLevelSet();

  for( char lyr = -2; lyr < 3; lyr++ )
  {
    SparseLevelSetType::LayerType layer = sparseLevelSet->GetLayer( lyr );
    SparseLevelSetType::LayerIterator lIt = layer.begin();

    std::cout << "*** " << static_cast< int >( lyr ) << " ***" <<std::endl;

    while( lIt != layer.end() )
    {
      std::cout << lIt->first << std::endl;
      ++lIt;
    }
    std::cout << std::endl;
  }

  typedef itk::LevelSetContainerBase< IdentifierType, SparseLevelSetType >
    LevelSetContainerType;
  typedef itk::LevelSetEquationTermContainerBase< InputImageType, LevelSetContainerType >
    TermContainerType;
  typedef itk::LevelSetEquationContainerBase< TermContainerType >
    EquationContainerType;

  typedef itk::UpdateWhitakerSparseLevelSet< Dimension, OutputPixelType, EquationContainerType >
    UpdateLevelSetType;
  UpdateLevelSetType::Pointer update_levelset = UpdateLevelSetType::New();
  update_levelset->SetInputLevelSet( sparseLevelSet );

  UpdateLevelSetType::LevelSetLayerType update_list;

  SparseLevelSetType::LayerIterator list_it = sparseLevelSet->GetLayer( 0 ).begin();
  SparseLevelSetType::LayerIterator list_end = sparseLevelSet->GetLayer( 0 ).end();

  typedef SparseLevelSetType::InputType LevelSetInputType;

  size_t k = 0;

  while( list_it != list_end )
  {
    LevelSetInputType tempIdx = list_it->first;
    if( atoi( argv[1]) == 2 )
    {
      update_list.insert(
            std::pair< LevelSetInputType, OutputPixelType >( tempIdx, -1.0  ) );
    }
    else
    {
      if( atoi( argv[1] ) == 0 )
      {
        update_list.insert(
              std::pair< LevelSetInputType, OutputPixelType >( tempIdx, 1.0  ) );
      }
      else
      {
        if( ( ( list_it->first )[1] % 20 ) < 10 )
        {
          update_list.insert(
                std::pair< LevelSetInputType, OutputPixelType >( tempIdx, -1.0  ) );
        }
        else
        {
          update_list.insert(
                std::pair< LevelSetInputType, OutputPixelType >( tempIdx, 1.0  ) );
        }
      }
    }
  ++k;
  ++list_it;
  }

  update_levelset->SetUpdate( update_list );
  update_levelset->Update();

  typedef itk::Image< OutputPixelType, Dimension >    OutputImageType;
  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( input->GetLargestPossibleRegion() );
  output->CopyInformation( input );
  output->Allocate();
  output->FillBuffer( 0.0 );

  typedef itk::Image< char, Dimension >               StatusImageType;
  StatusImageType::Pointer status = StatusImageType::New();
  status->SetRegions( input->GetLargestPossibleRegion() );
  status->CopyInformation( input );
  status->Allocate();
  status->FillBuffer( 0 );

  typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;
  OutputIteratorType oIt( output, output->GetLargestPossibleRegion() );
  oIt.GoToBegin();

  typedef itk::ImageRegionIteratorWithIndex< StatusImageType > StatusIteratorType;
  StatusIteratorType sIt( status, status->GetLargestPossibleRegion() );
  sIt.GoToBegin();

  StatusImageType::IndexType idx;

  while( !oIt.IsAtEnd() )
    {
    idx = oIt.GetIndex();
    oIt.Set( update_levelset->GetOutputLevelSet()->Evaluate( idx ) );
    sIt.Set( update_levelset->GetOutputLevelSet()->Status( idx ) );
    ++oIt;
    ++sIt;
    }

  typedef itk::ImageFileWriter< OutputImageType >     OutputWriterType;
  OutputWriterType::Pointer outputWriter = OutputWriterType::New();
  outputWriter->SetFileName( argv[2] );
  outputWriter->SetInput( output );

  try
    {
    outputWriter->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }

  typedef itk::ImageFileWriter< StatusImageType >     StatusWriterType;
  StatusWriterType::Pointer statusWriter = StatusWriterType::New();
  statusWriter->SetFileName( argv[3] );
  statusWriter->SetInput( status );

  try
    {
    statusWriter->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }

  return EXIT_SUCCESS;
}
