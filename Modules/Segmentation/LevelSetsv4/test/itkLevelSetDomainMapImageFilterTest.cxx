/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 817 $  // Revision of last commit
  Date: $Date: 2009-11-07 17:21:12 -0500 (Sat, 07 Nov 2009) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
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

#include "itkLevelSetDomainMapImageFilter.h"
#include <list>
#include <vector>

int itkLevelSetDomainMapImageFilterTest( int, char* [] )
{
  const unsigned int Dimension = 2;

  typedef std::list<int> ListPixelType;
  typedef std::list<int>::iterator ListIteratorType;
  typedef itk::Image< ListPixelType, Dimension >  InputImageType;
  typedef itk::Image< unsigned short, Dimension >  OutputImageType;

  typedef itk::LevelSetDomainMapImageFilter< InputImageType, OutputImageType >
    DomainMapImageFilterType;

  InputImageType::IndexType index;
  index[0] = 0;
  index[1] = 0;

  InputImageType::SizeType size;
  size[0] = 10;
  size[1] = 10;

  InputImageType::RegionType region;
  region.SetIndex( index );
  region.SetSize( size );

  ListPixelType l;

  InputImageType::Pointer input = InputImageType::New();
  input->SetRegions( region );
  input->Allocate();
  input->FillBuffer( l );

  for( unsigned int i = 0; i < 10; i++ )
  {
    ListPixelType ll;
    ll.push_back(i);
    ll.push_back(i+1);

    index[0] = index[1] = i;
    input->SetPixel( index, ll );
  }

  DomainMapImageFilterType::Pointer filter = DomainMapImageFilterType::New();
  filter->SetInput( input );
  filter->Update();

  OutputImageType::Pointer output = filter->GetOutput();

  itk::ImageRegionConstIteratorWithIndex<OutputImageType >
      it( output, output->GetLargestPossibleRegion() );

  it.GoToBegin();

  OutputImageType::IndexType out_index;
  OutputImageType::PixelType out_id;

  while( !it.IsAtEnd() )
    {
    out_index = it.GetIndex();
    out_id = it.Get();

    if( out_id > 0 )
      {
      std::cout << "*** " <<std::endl;
      std::cout << out_index << " # " << out_id <<std::endl;
      std::cout << filter->m_SetOfRegions[out_id];

      ListPixelType lout = filter->m_LevelSetList[out_id];
      if( lout.empty() )
        {
        return EXIT_FAILURE;
        }
      else
        {
        for( ListIteratorType lIt = lout.begin(); lIt != lout.end(); ++lIt )
          {
          std::cout << *lIt << " ";
          }
        std::cout << std::endl;
        }
      }
    ++it;
    }

  /*
  OutputImageType::Pointer output = filter->GetOutput();
  for( unsigned int i = 0; i < 10; i++ )
  {
    // Output the image
    index[0] = index[1] = i;
    std::cout << i << ' ' << output->GetPixel( index ) << std::endl;

    // Output the list of function ids
    ListPixelType lout = filter->m_LevelSetList[i];
    ListIteratorType lIt;
    for( lIt = lout.begin(); lIt != lout.end(); ++lIt )
    {
      std::cout <<  *lIt << std::endl;
    }

    // Output the regions
    std::cout << filter->m_SetOfRegions[i] << std::endl;
  }*/

  return EXIT_SUCCESS;
}
