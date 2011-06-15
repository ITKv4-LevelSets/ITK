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
#include <iostream>
#include <fstream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryImageToWhitakerSparseLevelSetAdaptor.h"
#include "itkUpdateWhitakerSparseLevelSet.h"
#include "itkImageRegionIterator.h"

int itkUpdateWhitakerSparseLevelSetTest2( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char InputPixelType;
  typedef double        OutputPixelType;

  typedef itk::Image< InputPixelType, Dimension >   InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

  typedef itk::ImageFileReader< InputImageType >  InputReaderType;
  typedef itk::ImageFileWriter< OutputImageType > OutputWriterType;

  typedef itk::BinaryImageToWhitakerSparseLevelSetAdaptor< InputImageType,
      OutputPixelType > BinaryToSparseAdaptorType;

  typedef BinaryToSparseAdaptorType::LevelSetType            SparseLevelSetType;
  typedef SparseLevelSetType::ImageType                      SparseImageType;
  typedef SparseLevelSetType::NodeAttributeType              NodeAttributeType;
  typedef BinaryToSparseAdaptorType::LevelSetNodeStatusType  StatusPixelType;

  typedef itk::Image< StatusPixelType, Dimension >  StatusImageType;
  typedef itk::ImageFileWriter< StatusImageType >   StatusWriterType;

  typedef itk::ImageRegionIterator< SparseImageType > SparseIteratorType;
  typedef itk::ImageRegionIterator< OutputImageType > OutputIteratorType;
  typedef itk::ImageRegionIterator< StatusImageType > StatusIteratorType;

  InputImageType::RegionType region;

  InputImageType::SizeType size;
  size.Fill( 51 );

  InputImageType::IndexType start;
  start.Fill( 0 );

  region.SetSize( size );
  region.SetIndex( start );

  // TODO: Spacing not taken into account in levelset
  InputImageType::SpacingType spacing;
  spacing.Fill( 0.5 );

  // Create a binary image
  InputImageType::Pointer input = InputImageType::New();
  input->SetRegions( region );
  input->SetSpacing( spacing );
  input->Allocate();
  input->FillBuffer( 0 );

  size.Fill( 30 );
  region.SetSize( size );

  start.Fill( 10 );
  region.SetIndex( start );

  typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;
  InputIteratorType i_it( input, region );
  i_it.GoToBegin();

  while( !i_it.IsAtEnd() )
    {
    i_it.Set( 255 );
    ++i_it;
    }
  std::cout << "Input image computed" << std::endl;

  BinaryToSparseAdaptorType::Pointer adaptor = BinaryToSparseAdaptorType::New();
  adaptor->SetInputImage( input );
  adaptor->Initialize();
  SparseLevelSetType::Pointer sparseLevelSet = adaptor->GetSparseLevelSet();
  std::cout << "Converted to sparse levelset" << std::endl;

  typedef itk::UpdateWhitakerSparseLevelSet< Dimension, OutputPixelType > UpdateLevelSetType;
  UpdateLevelSetType::UpdateListType* update_list =
      new UpdateLevelSetType::UpdateListType;

  SparseLevelSetType::NodeListIterator list_it = sparseLevelSet->GetListNode( 0 )->begin();
  SparseLevelSetType::NodeListIterator list_end = sparseLevelSet->GetListNode( 0 )->end();

  std::ifstream file;
  file.open( argv[1] );

  if( !file.is_open() )
    {
    std::cout << argv[1] <<" can't be opened" << std::endl;
    return EXIT_FAILURE;
    }

  OutputPixelType t;
  while( ( list_it != list_end ) && ( file.good() ) )
    {
    file >> t;
    update_list->push_back( t );
    ++list_it;
    }
  file.close();

  UpdateLevelSetType::Pointer update_levelset = UpdateLevelSetType::New();
  update_levelset->SetSparseLevelSet( sparseLevelSet );
  update_levelset->SetUpdate( update_list );
  update_levelset->SetDt( 1.0 );
  update_levelset->Update();

  delete update_list;

  OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( sparseLevelSet->GetOutputImage() );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }

  StatusWriterType::Pointer status_writer = StatusWriterType::New();
  status_writer->SetFileName( argv[3] );
  status_writer->SetInput( sparseLevelSet->GetStatusImage() );

  try
    {
    status_writer->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }

  return EXIT_SUCCESS;
}
