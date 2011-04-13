/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-11

 Copyright (c) 2009-11, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"

int itkMultiLevelSetImageTest( int , char* [] )
{
  const unsigned int Dimension = 2;

  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::LevelSetImageBase< ImageType > LevelSetType;

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

  ImageType::Pointer input1 = ImageType::New();
  input1->SetRegions( region );
  input1->Allocate();
  input1->FillBuffer( value );

  ImageType::Pointer input2 = ImageType::New();
  input2->SetRegions( region );
  input2->Allocate();
  input2->FillBuffer( value );

  ImageType::Pointer input3 = ImageType::New();
  input3->SetRegions( region );
  input3->Allocate();
  input3->FillBuffer( value );

  itk::ImageRegionIteratorWithIndex< ImageType > it1( input1,
                                              input1->GetLargestPossibleRegion() );

  itk::ImageRegionIteratorWithIndex< ImageType > it2( input2,
                                              input2->GetLargestPossibleRegion() );

  itk::ImageRegionIteratorWithIndex< ImageType > it3( input3,
                                              input3->GetLargestPossibleRegion() );

  typedef std::list< itk::IdentifierType > IdListType;

  typedef itk::Image< IdListType, Dimension > IdListImageType;

  IdListType list_ids;

  IdListImageType::Pointer id_image = IdListImageType::New();
  id_image->SetRegions( region );
  id_image->Allocate();
  id_image->FillBuffer( list_ids );

  it1.GoToBegin();
  it2.GoToBegin();
  it3.GoToBegin();

  ImageType::IndexType idx;

  while( !it1.IsAtEnd() )
    {
    idx = it1.GetIndex();
    list_ids.clear();

    if( ( idx[0] < 6 ) && ( idx[1] < 6 ) )
      {
      list_ids.push_back( 1 );
      }

    if( ( idx[0] > 1 ) && ( idx[1] > 1 ) &&
        ( idx[0] < 9 ) && ( idx[1] < 9 ) )
      {
      list_ids.push_back( 2 );
      }

    if( ( idx[0] > 4 ) && ( idx[1] > 4 ) )
      {
      list_ids.push_back( 3 );
      }

    id_image->SetPixel( idx, list_ids );

    it1.Set( vcl_sqrt(
             static_cast< float> ( ( idx[0] - 2 ) * ( idx[0] - 2 ) +
                                   ( idx[1] - 2 ) * ( idx[1] - 2 ) ) ) );

    it2.Set( vcl_sqrt(
             static_cast< float> ( ( idx[0] - 5 ) * ( idx[0] - 5 ) +
                                   ( idx[1] - 5 ) * ( idx[1] - 5 ) ) ) );

    it3.Set( vcl_sqrt(
             static_cast< float> ( ( idx[0] - 7 ) * ( idx[0] - 7 ) +
                                   ( idx[1] - 7 ) * ( idx[1] - 7 ) ) ) );

    ++it1;
    ++it2;
    ++it3;
    }

  std::map< itk::IdentifierType, LevelSetType::Pointer > level_set;
  level_set[1] = LevelSetType::New();
  level_set[1]->SetImage( input1 );

  level_set[2] = LevelSetType::New();
  level_set[2]->SetImage( input2 );

  level_set[3] = LevelSetType::New();
  level_set[3]->SetImage( input3 );

  typedef itk::Image< itk::IdentifierType, Dimension > CacheImageType;

  typedef itk::LevelSetDomainMapImageFilter< IdListImageType, CacheImageType >
    DomainMapImageFilterType;

  DomainMapImageFilterType::Pointer filter = DomainMapImageFilterType::New();
  filter->SetInput( id_image );
  filter->Update();

  CacheImageType::Pointer output = filter->GetOutput();

  itk::ImageRegionConstIteratorWithIndex<CacheImageType >
      it( output, output->GetLargestPossibleRegion() );

  it.GoToBegin();

  CacheImageType::IndexType out_index;
  CacheImageType::PixelType out_id;

  while( !it.IsAtEnd() )
    {
    out_index = it.GetIndex();
    out_id = it.Get();

    std::cout <<"***" << std::endl;
    std::cout << out_index <<std::endl;

    if( out_id != 0 )
      {
      IdListType lout = filter->m_LevelSetList[out_id];
      if( lout.empty() )
        {
        return EXIT_FAILURE;
        }
      else
        {
        for( IdListType::iterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
          {
          std::cout << *lIt <<" " << level_set[*lIt]->Evaluate( out_index )
                    << std::endl;
          }
        std::cout << std::endl;
        }
      }
    ++it;
    }


  return EXIT_SUCCESS;
}
