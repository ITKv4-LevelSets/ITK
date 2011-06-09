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
#include "itkBinaryImageToMalcolmSparseLevelSetAdaptor.h"
#include "itkImageRegionIterator.h"

int itkBinaryImageToMalcolmSparseLevelSetAdaptorTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char InputPixelType;
  typedef double        OutputPixelType;

  typedef itk::Image< InputPixelType, Dimension >                  InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >                 OutputImageType;

  typedef itk::ImageFileReader< InputImageType >                   InputReaderType;
  typedef itk::ImageFileWriter< OutputImageType >                  OutputWriterType;
  typedef itk::BinaryImageToMalcolmSparseLevelSetAdaptor< InputImageType >
    BinaryToSparseAdaptorType;

  typedef BinaryToSparseAdaptorType::SparseImageType  SparseImageType;
  typedef BinaryToSparseAdaptorType::LevelSetType     SparseLevelSetType;
  typedef itk::ImageRegionIterator< SparseImageType > SparseIteratorType;
  typedef itk::ImageRegionIterator< OutputImageType > OutputIteratorType;

  InputReaderType::Pointer reader = InputReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject& err )
    {
    std::cout << err << std::endl;
    }
  InputImageType::Pointer input = reader->GetOutput();
  std::cout << "Input image read" << std::endl;

  BinaryToSparseAdaptorType::Pointer adaptor = BinaryToSparseAdaptorType::New();
  adaptor->SetInputImage( input );
  adaptor->Initialize();
  std::cout << "Finished converting to sparse format" << std::endl;

  SparseLevelSetType::Pointer sparseLevelSet = adaptor->GetSparseLevelSet();
  SparseImageType::Pointer sparseImage = sparseLevelSet->GetImage();

  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( sparseImage->GetLargestPossibleRegion() );
  output->CopyInformation( sparseImage );
  output->Allocate();
  output->FillBuffer( 0.0 );

  typedef BinaryToSparseAdaptorType::LevelSetOutputType LevelSetOutputType;

  LevelSetOutputType p;
  SparseIteratorType sIt( sparseImage, sparseImage->GetLargestPossibleRegion() );
  sIt.GoToBegin();

  OutputIteratorType oIt( output, output->GetLargestPossibleRegion() );
  oIt.GoToBegin();
  while( !oIt.IsAtEnd() )
    {
    p = sIt.Get();
    oIt.Set( p );
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

  return EXIT_SUCCESS;
}
