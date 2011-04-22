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

#include "itkImage.h"
#include "itkListPixel.h"
#include "itkLevelSetImageBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"

int itkMultiLevelSetChanAndVeseInternalTermTest( int , char* [] )
{
  const unsigned int Dimension = 2;

  typedef unsigned char                                       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >             InputImageType;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > InputIteratorType;

  typedef float                                          PixelType;
  typedef itk::Image< PixelType, Dimension >             ImageType;
  typedef itk::LevelSetImageBase< ImageType >            LevelSetType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::IdentifierType                            IdentifierType;
  typedef itk::ListPixel< IdentifierType >               IdListType;
  typedef itk::Image< IdListType, Dimension >            IdListImageType;
  typedef itk::LevelSetDomainMapImageFilter< IdListImageType >
                                                         DomainMapImageFilterType;

  typedef itk::LevelSetContainerBase< IdentifierType, LevelSetType >  LevelSetContainerType;
  typedef itk::LevelSetEquationChanAndVeseInternalTerm< InputImageType, LevelSetContainerType > ChanAndVeseTermType;
  typedef itk::LevelSetEquationTermContainerBase< InputImageType, LevelSetContainerType > TermContainerType;

  typedef itk::AtanRegularizedHeavisideStepFunction< PixelType, PixelType >  HeavisideFunctionBaseType;

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

  InputImageType::Pointer input = InputImageType::New();
  input->SetRegions( region );
  input->Allocate();
  input->FillBuffer( 1 );

  ImageType::Pointer input1 = ImageType::New();
  input1->SetRegions( region );
  input1->Allocate();
  input1->FillBuffer( value );

  ImageType::Pointer input2 = ImageType::New();
  input2->SetRegions( region );
  input2->Allocate();
  input2->FillBuffer( value );

  ImageType::IndexType indices;
  IdListType list_ids;

  IdListImageType::Pointer id_image = IdListImageType::New();
  id_image->SetRegions( region );
  id_image->Allocate();
  id_image->FillBuffer( list_ids );

  IteratorType it1( input1, input1->GetLargestPossibleRegion() );
  IteratorType it2( input2, input2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();

  while( !it1.IsAtEnd() )
    {
    indices = it1.GetIndex();
    list_ids.clear();

    if( ( indices[0] < 5 ) && ( indices[1] < 5 ) )
      {
      list_ids.push_back( 1 );
      }

    if( ( indices[0] > 1 ) && ( indices[1] > 1 ) &&
        ( indices[0] < 8 ) && ( indices[1] < 8 ) )
      {
      list_ids.push_back( 2 );
      }

    id_image->SetPixel( indices, list_ids );

    it1.Set( vcl_sqrt(
             static_cast< float> ( ( indices[0] - 2 ) * ( indices[0] - 2 ) +
                                   ( indices[1] - 2 ) * ( indices[1] - 2 ) ) ) );

    it2.Set( vcl_sqrt(
             static_cast< float> ( ( indices[0] - 5 ) * ( indices[0] - 5 ) +
                                   ( indices[1] - 5 ) * ( indices[1] - 5 ) ) ) );
    ++it1;
    ++it2;
    }

  // Map of levelset bases
  std::map< itk::IdentifierType, LevelSetType::Pointer > level_set;
  level_set[1] = LevelSetType::New();
  level_set[1]->SetImage( input1 );

  level_set[2] = LevelSetType::New();
  level_set[2]->SetImage( input2 );

  // Insert the levelsets in a levelset container
  bool LevelSetNotYetAdded;
  LevelSetContainerType::Pointer lscontainer = LevelSetContainerType::New();
  LevelSetNotYetAdded = lscontainer->AddLevelSet( 0, level_set[1], false );
  if ( !LevelSetNotYetAdded )
  {
    return EXIT_FAILURE;
  }
  LevelSetNotYetAdded = lscontainer->AddLevelSet( 1, level_set[2], false );
  if ( !LevelSetNotYetAdded )
  {
    return EXIT_FAILURE;
  }
  std::cout << "Level set container created" << std::endl;

  // Define the Heaviside function
  HeavisideFunctionBaseType::Pointer heaviside = HeavisideFunctionBaseType::New();
  heaviside->SetEpsilon( 1.0 );

  // Create ChanAndVeseTerm
  ChanAndVeseTermType::Pointer cvTerm = ChanAndVeseTermType::New();
  cvTerm->SetHeaviside( heaviside );
  cvTerm->SetInput( input );
  cvTerm->SetCoefficient( 1.0 );
  cvTerm->SetCurrentLevelSet( 0 );
  cvTerm->SetLevelSetContainer( lscontainer );
  std::cout << "CV term created" << std::endl;

  TermContainerType::TermPointer temp = dynamic_cast< TermContainerType::TermType* >( cvTerm.GetPointer() );
  // Create Term Container
  TermContainerType::Pointer termContainer = TermContainerType::New();
  termContainer->AddTerm( 0, temp );
  std::cout << "Term container created" << std::endl;

  DomainMapImageFilterType::Pointer domainMapFilter = DomainMapImageFilterType::New();
  domainMapFilter->SetInput( id_image );
  domainMapFilter->Update();

  DomainMapImageFilterType::OutputImagePointer output = domainMapFilter->GetOutput();
  std::cout << "Domain partition computed" << std::endl;

  typedef DomainMapImageFilterType::OutputImageType       DomainOutputImageType;
  typedef DomainOutputImageType::LabelObjectType          LabelObjectType;
  typedef DomainOutputImageType::LabelObjectContainerType LabelObjectContainerType;

  LabelObjectContainerType label_objects = output->GetLabelObjectContainer();
  LabelObjectContainerType::const_iterator it = label_objects.begin();

  // 1-iterate on object
  while( it != label_objects.end() )
    {
    LabelObjectType* lo = it->second;

    LabelObjectType::LineContainerType lineContainer =
        lo->GetLineContainer();

    IdListType lout = lo->GetLabel();

    if( lout.empty() )
      {
      return EXIT_FAILURE;
      }

    // 2-iterate on lines
    for( LabelObjectType::LineContainerType::const_iterator
            lit = lineContainer.begin();
        lit != lineContainer.end(); ++lit )
      {
      DomainOutputImageType::IndexType firstIdx = lit->GetIndex();
      const DomainOutputImageType::OffsetValueType length = lit->GetLength();

      DomainOutputImageType::IndexValueType endIdx0 = firstIdx[0] + length;

      // 3-iterate on pixel (in a line)
      for ( DomainOutputImageType::IndexType idx = firstIdx;
            idx[0] < endIdx0;
            ++idx[0] )
        {
        std::cout << idx << " ";

        for( IdListType::iterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
          {
          std::cout << *lIt << " ";
          LevelSetType::Pointer levelSet = lscontainer->GetLevelSet( *lIt - 1);
          std::cout << levelSet->Evaluate( idx ) <<std::endl;

          if ( *lIt - 1 == 0 )
            {
            ChanAndVeseTermType::Pointer eqTerm = dynamic_cast< ChanAndVeseTermType* >(
                termContainer->GetTerm( *lIt - 1 ).GetPointer() );
            std::cout << eqTerm->Evaluate( idx ) <<std::endl;
            }
          }
        std::cout << std::endl;
        }
      }
    ++it;
    }


//  typedef std::map< itk::IdentifierType, DomainMapImageFilterType::NounToBeDefined >::iterator DomainMapIterator;
//  DomainMapIterator map_it = domainMapFilter->m_LevelSetMap.begin();
//  DomainMapIterator map_end = domainMapFilter->m_LevelSetMap.end();

//  LevelSetType::Pointer levelSet;
//  ChanAndVeseTermType::Pointer eqTerm;
//  while( map_it != map_end )
//    {
//    IdListImageType::RegionType temp_region = map_it->second.m_Region;

//    itk::ImageRegionConstIteratorWithIndex<IdListImageType >
//        temp_it( id_image, temp_region );
//    temp_it.GoToBegin();

//    while( !temp_it.IsAtEnd() )
//      {
//      std::cout << temp_it.GetIndex() << std::endl;
//      IdListType lout = map_it->second.m_List;

//      if( lout.empty() )
//        {
//        return EXIT_FAILURE;
//        }

//      for( IdListType::iterator lIt = lout.begin(); lIt != lout.end(); ++lIt )
//        {
//          std::cout << *lIt << " ";
//          levelSet = lscontainer->GetLevelSet( *lIt - 1);
//          std::cout << levelSet->Evaluate( temp_it.GetIndex() ) << std::endl;

//          if ( *lIt - 1 == 0 )
//          {
//            eqTerm =  dynamic_cast< ChanAndVeseTermType* >( termContainer->GetTerm( *lIt - 1 ).GetPointer() );
//            std::cout << eqTerm->Evaluate( temp_it.GetIndex() ) << std::endl;
//          }
//        }
//      std::cout << std::endl;
//      ++temp_it;
//      }
//    ++map_it;
//    }

// Update the term

  return EXIT_SUCCESS;
}
