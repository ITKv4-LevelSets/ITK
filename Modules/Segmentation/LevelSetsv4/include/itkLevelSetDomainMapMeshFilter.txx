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

#ifndef __itkLevelSetDomainMapMeshFilter_txx
#define __itkLevelSetDomainMapMeshFilter_txx

#include "itkLevelSetDomainMapMeshFilter.h"

namespace itk
{
template < class TInputMesh, class TOutputMesh >
LevelSetDomainMapMeshFilter< TInputMesh, TOutputMesh >
::LevelSetDomainMapMeshFilter()
{
  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, OutputMeshType::New() );
}


template < class TInputMesh, class TOutputMesh >
void
LevelSetDomainMapMeshFilter< TInputMesh, TOutputMesh >::
GenerateData()
{
  InputMeshConstPointer input = this->GetInput();

  OutputMeshPointer output = this->GetOutput();

  OutputPointDataContainerPointer output_data = output->GetPointData();
  OutputPointDataContainerIterator o_it = output_data->Begin();

  while( o_it != output_data->End() )
    {
    o_it->Value() = NumericTraits< OutputMeshPixelType >::Zero;
    ++o_it;
    }

  InputPointDataContainerPointer input_data = input->GetPointData();
  InputPointDataContainerIterator i_it = input_data->Begin();

  o_it = output_data->Begin();

  InputPointDataContainerIterator i_end = input_data->End();
  --i_end;

  InputPointIdentifier end_idx = i_end->Index();
  i_end = input_data->End();

  IdentifierType segmentId = 1;

  while( i_it != i_end )
    {
    InputPointIdentifier startIdx = i_it->Index();
    InputPointIdentifier stopIdx = startIdx;

    InputMeshPixelType inputPixel = i_it->Value();
    OutputMeshPixelType outputPixel = o_it->Value();

    // outputPixel is null when it has not been processed yet,
    // or there is nothing to be processed
    if ( ( !inputPixel.empty() ) && ( outputPixel == 0 ) )
      {
      bool flag = true;
      stopIdx = startIdx;

      while ( ( flag ) && ( stopIdx <= end_idx ) )
        {
        InputMeshPixelType nextPixel;
        OutputMeshPixelType currentOutputPixel;

        input->GetPointData( stopIdx, &nextPixel );

        // Check if the input list pixels are different, or
        // the output image already has been assigned to another region
        if ( ( nextPixel != inputPixel ) ||
             ( currentOutputPixel != NumericTraits< OutputMeshPixelType >::Zero ) )
          {
          flag = false;
          }
        else
          {
          output->SetPointData( stopIdx, segmentId );
          ++stopIdx;
          }
        }
      }

    m_LevelSetMap[segmentId] = NounToBeDefined( segmentId, startIdx,
                                                stopIdx, inputPixel );

    ++i_it;
    ++o_it;
    }

  this->GraftOutput ( output );
}

template < class TInputMesh, class TOutputMesh >
void
LevelSetDomainMapMeshFilter< TInputMesh, TOutputMesh >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:        " << GetNameOfClass() << std::endl;
}

} /* end namespace itk */

#endif
