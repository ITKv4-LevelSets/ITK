/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1658 $  // Revision of last commit
  Date: $Date: 2010-06-14 15:49:25 -0400 (Mon, 14 Jun 2010) $  // Date of last commit
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

#ifndef __itkLevelSetDomainMapImageFilter_txx
#define __itkLevelSetDomainMapImageFilter_txx

#include "itkLevelSetDomainMapImageFilter.h"

namespace itk
{
template < class TInputImage, class TOutputImage >
LevelSetDomainMapImageFilter< TInputImage, TOutputImage >
::LevelSetDomainMapImageFilter()
{
  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0,TOutputImage::New() );
}

template < class TInputImage, class TOutputImage >
void
LevelSetDomainMapImageFilter< TInputImage, TOutputImage >::
GenerateData()
{
  InputImageConstPointer input = this->GetInput();
  InputImageSizeType size = input->GetLargestPossibleRegion().GetSize();

  OutputImagePointer output = OutputImageType::New();
  output->CopyInformation( input );
  output->SetRegions( input->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  InputImagePixelType inputPixel, nextPixel;
  OutputImagePixelType outputPixel;
  InputImageIndexType startIdx, stopIdx;
  InputImageIndexType end;
  InputImageSizeType sizeOfRegion;
  InputImageRegionType subRegion;
  unsigned int segmentId = 1;
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    end[i] = size[i] - 1;
  }

  InputConstIteratorType iIt( input, input->GetLargestPossibleRegion() );
  OutputIndexIteratorType oIt( output, output->GetLargestPossibleRegion() );

  iIt.GoToBegin();
  oIt.GoToBegin();
  while( !iIt.IsAtEnd() )
  {
    startIdx = iIt.GetIndex();
    stopIdx = startIdx;
    inputPixel = iIt.Get();
    outputPixel = oIt.Get();

    if ( ( inputPixel.size() != 0 ) && ( outputPixel == 0 ) )
    {
      for( unsigned int i = 0; i < ImageDimension; i++ )
      {
        bool flag = true;
        stopIdx = startIdx;
        while ( ( flag ) && ( stopIdx[i] <= end[i] ) )
        {
          nextPixel = input->GetPixel( stopIdx );
          if ( nextPixel != inputPixel )
          {
            flag = false;
          }
          else
          {
            stopIdx[i]++;
          }
        }
        sizeOfRegion[i] = stopIdx[i] - startIdx[i];
      }

      subRegion.SetSize( sizeOfRegion );
      subRegion.SetIndex( startIdx );
      m_SetOfRegions.push_back( subRegion );
      m_LevelSetList.push_back( inputPixel );

      OutputIndexIteratorType ooIt( output, subRegion );
      ooIt.GoToBegin();
      while( !ooIt.IsAtEnd() )
      {
        ooIt.Set( segmentId );
        ++ooIt;
      }
      segmentId++;
    }
    ++iIt;
    ++oIt;
  }

  this->GraftOutput ( output );
}

template < class TInputImage, class TOutputImage >
void
LevelSetDomainMapImageFilter< TInputImage, TOutputImage >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:        " << GetNameOfClass() << std::endl;
}

} /* end namespace itk */

#endif
