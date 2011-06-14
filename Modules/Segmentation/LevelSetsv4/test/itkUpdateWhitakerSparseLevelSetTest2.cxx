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

  InputImageType::Pointer input = InputImageType::New();
  InputImageType::RegionType region;

//   InputImageType::SizeType size;
//   size.Fill( 25 );
//
//   InputImageType::IndexType start;
//   start.Fill( 0 );
//
//   region.SetSize( size );
//   region.SetIndex( start );
//   input->SetRegions( region );
//   input->Allocate();
//   input->FillBuffer( 0 );
//
//   size[0] = 10;
//   region.SetSize( size );
//
//   typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;
//   InputIteratorType i_it( input, region );
//   i_it.GoToBegin();
//
//   while( !i_it.IsAtEnd() )
//     {
//     i_it.Set( 255 );
//     ++i_it;
//     }

  InputImageType::SizeType size;
  size.Fill( 51 );

  InputImageType::IndexType start;
  start.Fill( 0 );

  region.SetSize( size );
  region.SetIndex( start );
  input->SetRegions( region );
  input->Allocate();
  input->FillBuffer( 0 );

  size[0] = 30;
  size[1] = 30;
  region.SetSize( size );

  start[0] = 10;
  start[1] = 10;
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
  std::cout << "Finished converting to sparse format" << std::endl;

  SparseLevelSetType::Pointer sparseLevelSet = adaptor->GetSparseLevelSet();

  typedef itk::UpdateWhitakerSparseLevelSet< Dimension, OutputPixelType > UpdateLevelSetType;
  UpdateLevelSetType::Pointer update_levelset = UpdateLevelSetType::New();
  update_levelset->SetSparseLevelSet( sparseLevelSet );

  UpdateLevelSetType::UpdateListType* update_list =
      new UpdateLevelSetType::UpdateListType;

  SparseLevelSetType::NodeListIterator list_it = sparseLevelSet->GetListNode( 0 )->begin();
  SparseLevelSetType::NodeListIterator list_end = sparseLevelSet->GetListNode( 0 )->end();

  size_t k = 0;

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

//     if( atoi( argv[1]) == 2 )
//       {
//       update_list->push_back( -1. );
//       }
//     else
//       {
//       if( atoi( argv[1] ) == 0 )
//         {
//         update_list->push_back( 1. );
//         }
//       else
//         {
//         if( ( ( list_it->first )[1] % 20 ) < 10 )
//           {
//           update_list->push_back( -1. );
//           }
//         else
//           {
//           update_list->push_back( 1. );
//           }
//         }
//       }
//     ++k;

  update_levelset->SetUpdate( update_list );
  update_levelset->Update();

  delete update_list;

  NodeAttributeType p;
  SparseIteratorType ls_It( sparseLevelSet->GetImage(),
                         sparseLevelSet->GetImage()->GetLargestPossibleRegion() );
  ls_It.GoToBegin();

  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( sparseLevelSet->GetImage()->GetLargestPossibleRegion() );
  output->CopyInformation( sparseLevelSet->GetImage() );
  output->Allocate();
  output->FillBuffer( 0.0 );

  StatusImageType::Pointer status = StatusImageType::New();
  status->SetRegions( sparseLevelSet->GetImage()->GetLargestPossibleRegion() );
  status->CopyInformation( sparseLevelSet->GetImage() );
  status->Allocate();
  status->FillBuffer( 0 );

  OutputIteratorType oIt( output,
                          output->GetLargestPossibleRegion() );
  oIt.GoToBegin();

  StatusIteratorType sIt( status, status->GetLargestPossibleRegion() );
  sIt.GoToBegin();

  while( !oIt.IsAtEnd() )
    {
    p = ls_It.Get();
    oIt.Set( p.m_Value );
    sIt.Set( p.m_Status );
    ++ls_It;
    ++oIt;
    ++sIt;
    }

  OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );

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
  status_writer->SetInput( status );

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
