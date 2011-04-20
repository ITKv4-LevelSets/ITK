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

#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkLevelSetQuadEdgeMeshBase.h"
#include "itkLevelSetDomainMapMeshFilter.h"

int itkMultiLevelSetMeshTest( int , char* [] )
{
  const unsigned int Dimension = 3;

  typedef float                                           PixelType;
  typedef itk::QuadEdgeMesh< PixelType, Dimension >       MeshType;
  typedef itk::LevelSetQuadEdgeMeshBase< MeshType >       LevelSetType;
  typedef std::list< itk::IdentifierType >                IdListType;
  typedef itk::QuadEdgeMesh< IdListType, Dimension >      IdListMeshType;
  typedef itk::LevelSetDomainMapMeshFilter< IdListMeshType >
                                                          DomainMapMeshFilterType;

  IdListMeshType::PointType center;
  center.Fill( 0. );

  typedef itk::RegularSphereMeshSource< IdListMeshType > SphereSourceType;
  SphereSourceType::Pointer sphere_filter = SphereSourceType::New();
  sphere_filter->SetCenter( center );
  sphere_filter->SetResolution( 3 );
  sphere_filter->Update();

  IdListMeshType::Pointer id_Mesh = sphere_filter->GetOutput();

  IdListType l;

  IdListMeshType::PointDataContainerPointer
      point_data = id_Mesh->GetPointData();

  if( point_data.IsNull() )
    {
    point_data = IdListMeshType::PointDataContainer::New();
    point_data->Reserve( id_Mesh->GetNumberOfPoints() );
    id_Mesh->SetPointData( point_data );
    }

  IdListMeshType::PointDataContainerIterator
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

  typedef itk::RegularSphereMeshSource< MeshType > SphereSourceType2;
  SphereSourceType2::Pointer sphere_filter2 = SphereSourceType2::New();
  sphere_filter2->SetCenter( center );
  sphere_filter2->SetResolution( 3 );
  sphere_filter2->Update();

  MeshType::Pointer input1 = sphere_filter2->GetOutput();

{
  MeshType::PointDataContainerPointer data = input1->GetPointData();
  MeshType::PointDataContainerIterator d_it = data->Begin();

  while( d_it != data->End() )
    {
    float value = static_cast< float >( i_it->Index() );
    data->SetElement( d_it->Index(), value );
    ++d_it;
    }

  input1->SetPointData( data );
}

  SphereSourceType2::Pointer sphere_filter3 = SphereSourceType2::New();
  sphere_filter3->SetCenter( center );
  sphere_filter3->SetResolution( 3 );
  sphere_filter3->Update();

  MeshType::Pointer input2 = sphere_filter3->GetOutput();

  {
    MeshType::PointDataContainerPointer data = input2->GetPointData();
    MeshType::PointDataContainerIterator d_it = data->Begin();

    while( d_it != data->End() )
      {
      float value = static_cast< float >( i_it->Index() * i_it->Index() );
      data->SetElement( d_it->Index(), value );
      ++d_it;
      }

  input2->SetPointData( data );
  }


  std::map< itk::IdentifierType, LevelSetType::Pointer > level_set;
  level_set[1] = LevelSetType::New();
  level_set[1]->SetMesh( input1 );

  level_set[2] = LevelSetType::New();
  level_set[2]->SetMesh( input2 );

  DomainMapMeshFilterType::Pointer filter = DomainMapMeshFilterType::New();
  filter->SetInput( id_Mesh );
  filter->Update();

  std::map< IdListMeshType::PointIdentifier, itk::IdentifierType >
      output = filter->m_Output;

  std::map< IdListMeshType::PointIdentifier, itk::IdentifierType >::iterator
      it = output.begin();

  IdListMeshType::PointIdentifier out_index;
  itk::IdentifierType out_id;

  while( it != output.end() )
    {
    out_index = it->first;
    out_id = it->second;

    IdListType solution;
    if( out_index < 7 )
      {
      solution.push_back( 1 );
      }
    if( out_index > 3 )
      {
      solution.push_back( 2 );
      }
    solution.sort();

    std::cout <<"***" << std::endl;
    std::cout << out_index <<std::endl;

    if( out_id != 0 )
      {
      IdListType lout = filter->m_LevelSetMap[out_id].m_List;
      std::cout << filter->m_LevelSetMap[out_id].m_Start << std::endl;
      std::cout << filter->m_LevelSetMap[out_id].m_End << std::endl;

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

        lout.sort();
        if( lout != solution )
          {
          std::cout <<"FAILURE!!!" <<std::endl;
          return EXIT_FAILURE;
          }
        }
      }
    ++it;
    }

  std::map< itk::IdentifierType, DomainMapMeshFilterType::NounToBeDefined >::iterator map_it = filter->m_LevelSetMap.begin();
  std::map< itk::IdentifierType, DomainMapMeshFilterType::NounToBeDefined >::iterator map_end = filter->m_LevelSetMap.end();

  while( map_it != map_end )
    {
    IdListMeshType::PointIdentifier temp_start = map_it->second.m_Start;
    IdListMeshType::PointIdentifier temp_end = map_it->second.m_End;

    IdListMeshType::PointIdentifier temp_id = temp_start;

    while( temp_id != temp_end )
      {
      std::cout << temp_id << std::endl;
      IdListType lout = map_it->second.m_List;

      if( lout.empty() )
        {
        return EXIT_FAILURE;
        }

      for( IdListType::iterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
        {
        std::cout << *lIt <<" " << level_set[*lIt]->Evaluate( temp_id )
                  << std::endl;
        }
      std::cout << std::endl;
      ++temp_id;
      }
    ++map_it;
    }

  return EXIT_SUCCESS;
}
