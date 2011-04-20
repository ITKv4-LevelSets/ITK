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
#include "itkNumericTraits.h"

namespace itk
{
template < class TInputMesh >
LevelSetDomainMapMeshFilter< TInputMesh >
::LevelSetDomainMapMeshFilter()
{}


template < class TInputMesh >
void
LevelSetDomainMapMeshFilter< TInputMesh >::
GenerateData()
{
  InputPointDataContainerConstPointer input_data = m_Input->GetPointData();
  PointDataContainerConstIterator i_it = input_data->Begin();

  PointDataContainerConstIterator i_end = input_data->End();

  while( i_it != i_end )
    {
    m_Output[ i_it->Index() ] = 0;
    ++i_it;
    }

  --i_end;

  InputPointIdentifier end_idx = i_end->Index();
  i_end = input_data->End();

  IdentifierType segmentId = 1;

  while( i_it != i_end )
    {
    InputPointIdentifier startIdx = i_it->Index();
    InputPointIdentifier stopIdx = startIdx;

    InputMeshPixelType inputPixel = i_it->Value();
    IdentifierType outputPixel = m_Output[startIdx];

    // outputPixel is null when it has not been processed yet,
    // or there is nothing to be processed
    if ( ( !inputPixel.empty() ) && ( outputPixel == 0 ) )
      {
      bool flag = true;
      stopIdx = startIdx;

      while ( ( flag ) && ( stopIdx <= end_idx ) )
        {
        InputMeshPixelType nextPixel;
        IdentifierType currentOutputPixel = m_Output[stopIdx];

        m_Input->GetPointData( stopIdx, &nextPixel );

        // Check if the input list pixels are different, or
        // the output image already has been assigned to another region
        if ( ( nextPixel != inputPixel ) ||
             ( currentOutputPixel != 0 ) )
          {
          flag = false;
          }
        else
          {
          m_Output[ stopIdx ] = segmentId;
          ++stopIdx;
          }
        }
      }

    m_LevelSetMap[segmentId] = NounToBeDefined( segmentId, startIdx,
                                                stopIdx, inputPixel );

    ++i_it;
    }
}

template < class TInputMesh >
void
LevelSetDomainMapMeshFilter< TInputMesh >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:        " << GetNameOfClass() << std::endl;
}

} /* end namespace itk */

#endif
