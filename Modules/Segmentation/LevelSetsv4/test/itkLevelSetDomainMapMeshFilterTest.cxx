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

#include "itkLevelSetDomainMapMeshFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include <list>
#include <vector>

int itkLevelSetDomainMapMeshFilterTest( int, char* [] )
{
  const unsigned int Dimension = 3;

  typedef std::list<int>            ListPixelType;
  typedef std::list<int>::iterator  ListIteratorType;

  typedef itk::QuadEdgeMesh< ListPixelType, Dimension >   InputMeshType;

  typedef itk::LevelSetDomainMapMeshFilter< InputMeshType >
    DomainMapMeshFilterType;

  InputMeshType::PointType center;
  center.Fill( 0. );

  typedef itk::RegularSphereMeshSource< InputMeshType > SphereSourceType;
  SphereSourceType::Pointer sphere_filter = SphereSourceType::New();
  sphere_filter->SetCenter( center );
  sphere_filter->SetResolution( 3 );
  sphere_filter->Update();

  InputMeshType::Pointer input = sphere_filter->GetOutput();

  ListPixelType l;

  InputMeshType::PointDataContainerPointer
      point_data = input->GetPointData();

  if( point_data.IsNull() )
    {
    point_data = InputMeshType::PointDataContainer::New();
    point_data->Reserve( input->GetNumberOfPoints() );
    input->SetPointData( point_data );
    }

  InputMeshType::PointDataContainerIterator
      i_it = point_data->Begin();

  while( i_it != point_data->End() )
    {
    l.clear();
    if( i_it->Index() < 7 )
      {
      l.push_back( 1 );
      }
    if( i_it->Index() > 3 )
      {
      l.push_back( 2 );
      }
    point_data->SetElement( i_it->Index(), l );
    ++i_it;
    }

  DomainMapMeshFilterType::Pointer filter = DomainMapMeshFilterType::New();
  filter->SetInput( input );
  filter->Update();

  std::map< InputMeshType::PointIdentifier, itk::IdentifierType >
      output = filter->m_Output;

  std::map< InputMeshType::PointIdentifier, itk::IdentifierType >::const_iterator
      o_it = output.begin();

  while( o_it != output.end() )
    {
    InputMeshType::PointIdentifier out_index = o_it->first;
    unsigned short out_id = o_it->second;

    if( out_id > 0 )
      {
      std::cout << "*** " <<std::endl;
      std::cout << out_index << " # " << out_id <<std::endl;
      std::cout << filter->m_LevelSetMap[out_id].m_Start << std::endl;
      std::cout << filter->m_LevelSetMap[out_id].m_End << std::endl;

      ListPixelType lout = filter->m_LevelSetMap[out_id].m_List;
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
    ++o_it;
    }

  return EXIT_SUCCESS;
}
