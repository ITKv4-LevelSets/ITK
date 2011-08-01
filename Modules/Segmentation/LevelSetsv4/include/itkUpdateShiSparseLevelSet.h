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

#ifndef __itkUpdateShiSparseLevelSet_h
#define __itkUpdateShiSparseLevelSet_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"
#include "itkShiSparseLevelSetBase.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include <list>
#include "itkObject.h"

namespace itk
{
template< unsigned int VDimension, class TEquationContainer >
class UpdateShiSparseLevelSet : public Object
{
public:
  typedef UpdateShiSparseLevelSet       Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( UpdateShiSparseLevelSet, Object );

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef ShiSparseLevelSetBase< ImageDimension >      LevelSetType;
  typedef typename LevelSetType::Pointer               LevelSetPointer;
  typedef typename LevelSetType::InputType             LevelSetInputType;
  typedef typename LevelSetType::OutputType            LevelSetOutputType;
  typedef typename LevelSetType::OutputRealType        LevelSetOutputRealType;

  typedef typename LevelSetType::LabelMapType          LevelSetLabelMapType;
  typedef typename LevelSetType::LabelMapPointer       LevelSetLabelMapPointer;

  typedef typename LevelSetType::LabelObjectType       LevelSetLabelObjectType;
  typedef typename LevelSetType::LabelObjectPointer    LevelSetLabelObjectPointer;
  typedef typename LevelSetType::LabelObjectLengthType LevelSetLabelObjectLengthType;
  typedef typename LevelSetType::LabelObjectLineType   LevelSetLabelObjectLineType;

  typedef typename LevelSetType::LayerType             LevelSetLayerType;
  typedef typename LevelSetType::LayerIterator         LevelSetLayerIterator;
  typedef typename LevelSetType::LayerConstIterator    LevelSetLayerConstIterator;

  typedef typename LevelSetType::LayerMapType           LevelSetLayerMapType;
  typedef typename LevelSetType::LayerMapIterator       LevelSetLayerMapIterator;
  typedef typename LevelSetType::LayerMapConstIterator  LevelSetLayerMapConstIterator;

  typedef TEquationContainer                      EquationContainerType;
  typedef typename EquationContainerType::Pointer EquationContainerPointer;

  itkGetObjectMacro( OutputLevelSet, LevelSetType );

  void Update()
  {
    if( m_SparseLevelSet.IsNull() )
      {
      itkGenericExceptionMacro( <<"m_SparseLevelSet is NULL" );
      }

    if( this->m_SingleLevelSet != true )
      {
      itkWarningMacro( << "the current implementation does not allow to run multi object level sts");
      }

    m_OutputLevelSet->SetLayer( -1, m_SparseLevelSet->GetLayer( -1 ) );
    m_OutputLevelSet->SetLayer(  1, m_SparseLevelSet->GetLayer(  1 ) );

    m_OutputLevelSet->SetLabelMap( m_SparseLevelSet->GetLabelMap() );

    typedef LabelMapToLabelImageFilter<LevelSetLabelMapType, LabelImageType> LabelMapToLabelImageFilterType;
    typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput( m_SparseLevelSet->GetLabelMap() );
    labelMapToLabelImageFilter->Update();

    m_InternalImage = labelMapToLabelImageFilter->GetOutput();
    m_InternalImage->DisconnectPipeline();


    // neighborhood iterator
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

    // Step 2.1.1
    UpdateL_out();

    // Step 2.1.2 - for each point x in L_in
    LevelSetLayerType & list_in = m_OutputLevelSet->GetLayer( -1 );

    LevelSetLayerIterator nodeIt = list_in.begin();
    LevelSetLayerIterator nodeEnd = list_in.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool to_be_deleted = true;

      for( typename NeighborhoodIteratorType::Iterator
           i = neighIt.Begin();
           !i.IsAtEnd(); ++i )
        {
        if ( i.Get() > NumericTraits< LevelSetOutputType >::Zero )
          {
          to_be_deleted = false;
          break;
          }
        }
      if( to_be_deleted )
        {
//         std::cout << p.first << std::endl;
        LevelSetOutputType oldValue = -1;
        LevelSetOutputType newValue = -3;
        this->m_InternalImage->SetPixel( currentIndex, newValue );

        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        list_in.erase( tempIt );

        m_EquationContainer->GetEquation( m_CurrentLevelSetId )->UpdatePixel(
              currentIndex,
              static_cast< LevelSetOutputRealType >( oldValue ),
              static_cast< LevelSetOutputRealType >( newValue ) );
        }
      else
        {
        ++nodeIt;
        }
      }

    // Step 2.1.3 - for each point x in L_in
    UpdateL_in();

    // Step 2.1.4
    LevelSetLayerType & list_out = m_OutputLevelSet->GetLayer( 1 );

    nodeIt = list_out.begin();
    nodeEnd = list_out.end();

    while( nodeIt != nodeEnd )
      {
      LevelSetInputType currentIndex = nodeIt->first;

      neighIt.SetLocation( currentIndex );

      bool to_be_deleted = true;

      for( typename NeighborhoodIteratorType::Iterator
              i = neighIt.Begin();
          !i.IsAtEnd(); ++i )
        {
        if ( i.Get() < NumericTraits< LevelSetOutputType >::Zero )
          {
          to_be_deleted = false;
          break;
          }
        }
      if( to_be_deleted )
        {
        LevelSetOutputType oldValue = 1;
        LevelSetOutputType newValue = 3;
        this->m_InternalImage->SetPixel( currentIndex, newValue );

        LevelSetLayerIterator tempIt = nodeIt;
        ++nodeIt;
        list_out.erase( tempIt );

        m_EquationContainer->GetEquation( m_CurrentLevelSetId )->UpdatePixel(
              currentIndex,
              static_cast< LevelSetOutputRealType >( oldValue ),
              static_cast< LevelSetOutputRealType >( newValue )
              );
        }
      else
        {
        ++nodeIt;
        }
      }

    typedef LabelImageToLabelMapFilter< LabelImageType, LevelSetLabelMapType> LabelImageToLabelMapFilterType;
    typename LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();
    labelImageToLabelMapFilter->SetInput( m_InternalImage );
    labelImageToLabelMapFilter->SetBackgroundValue( 3 );
    labelImageToLabelMapFilter->Update();

    m_OutputLevelSet->GetLabelMap( )->Graft( labelImageToLabelMapFilter->GetOutput() );
  }

  // Set/Get the sparse levet set image
  itkSetObjectMacro( SparseLevelSet, LevelSetType );
  itkGetObjectMacro( SparseLevelSet, LevelSetType );

  itkGetMacro( RMSChangeAccumulator, LevelSetOutputRealType );

  // set the term container
  itkSetObjectMacro( EquationContainer, EquationContainerType );
  itkGetObjectMacro( EquationContainer, EquationContainerType );

  itkSetMacro( CurrentLevelSetId, IdentifierType );
  itkGetMacro( CurrentLevelSetId, IdentifierType );

  itkSetMacro( SingleLevelSet, bool );
  itkGetMacro( SingleLevelSet, bool );

protected:
  UpdateShiSparseLevelSet() :
    m_RMSChangeAccumulator( NumericTraits< LevelSetOutputRealType >::Zero ),
    m_SingleLevelSet( true )
    {
    m_OutputLevelSet = LevelSetType::New();
    }
  ~UpdateShiSparseLevelSet() {}

  // input
  LevelSetPointer   m_SparseLevelSet;

  // output
  LevelSetPointer   m_OutputLevelSet;

  IdentifierType           m_CurrentLevelSetId;
  LevelSetOutputRealType   m_RMSChangeAccumulator;
  EquationContainerPointer m_EquationContainer;
  bool                     m_SingleLevelSet;

  typedef Image< char, ImageDimension >     LabelImageType;
  typedef typename LabelImageType::Pointer  LabelImagePointer;

  LabelImagePointer m_InternalImage;

  typedef ShapedNeighborhoodIterator< LabelImageType > NeighborhoodIteratorType;

  // this is the same as Procedure 2
  // Input is a update image point m_UpdateImage
  // Input is also ShiSparseLevelSetBasePointer
  void UpdateL_out()
  {
    std::cout << "UpdateL_out" << std::endl;
    LevelSetLayerType & list_out  = m_SparseLevelSet->GetLayer(  1 );
    LevelSetLayerType & list_in   = m_SparseLevelSet->GetLayer( -1 );

    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

//     if( m_Update[1]->size() != list_out->size() )
//       {
//       itkGenericExceptionMacro( "m_Update[1]->size() != list_out->size()" );
//       }

    LevelSetLayerIterator nodeIt   = list_out.begin();
    LevelSetLayerIterator nodeEnd  = list_out.end();

    // for each point in Lz
    while( nodeIt != nodeEnd )
      {
      bool erased = false;
      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;

      // update the level set
      LevelSetOutputRealType update =
          m_EquationContainer->GetEquation( m_CurrentLevelSetId )->Evaluate( currentIndex );

//       std::cout << p.first << ' ' << int(p.second) << ' ' << update << std::endl;

      if( update < NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( Con( currentIndex, currentValue , update ) )
          {
          // CheckIn
          LevelSetOutputType oldValue = currentValue;
          LevelSetOutputType newValue = -1;

          this->m_InternalImage->SetPixel( currentIndex, newValue );

          m_EquationContainer->UpdatePixel( currentIndex, oldValue , newValue );

          list_in.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, newValue ) );

          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          list_out.erase( tempIt );
          erased = true;

          neighIt.SetLocation( currentIndex );

          for( typename NeighborhoodIteratorType::Iterator
              i = neighIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            if ( i.Get() == 3 )
              {
              LevelSetInputType tempIndex =
                  neighIt.GetIndex( i.GetNeighborhoodOffset() );

              oldValue = 3;
              newValue = 1;

              list_out.insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( tempIndex, newValue ) );

              this->m_InternalImage->SetPixel( tempIndex, newValue );

              m_EquationContainer->UpdatePixel( tempIndex, oldValue , newValue );
              }
            }
          }
        }
      if( !erased )
        {
        ++nodeIt;
        }
      }
    }

  void UpdateL_in()
  {
    std::cout << "UpdateL_in" << std::endl;
    LevelSetLayerType & list_out  = m_SparseLevelSet->GetLayer(  1 );
    LevelSetLayerType & list_in   = m_SparseLevelSet->GetLayer( -1 );

    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

    LevelSetLayerIterator nodeIt   = list_in.begin();
    LevelSetLayerIterator nodeEnd  = list_in.end();

    // for each point in Lz
    while( nodeIt != nodeEnd )
      {
      bool erased = false;
      LevelSetInputType   currentIndex = nodeIt->first;
      LevelSetOutputType  currentValue = nodeIt->second;

      // TODO
      // update for the current level set
      LevelSetOutputRealType update =
          m_EquationContainer->GetEquation( m_CurrentLevelSetId )->Evaluate( currentIndex );
//       std::cout << p.first << ' ' << int(p.second) << ' ' << update << std::endl;

      if( update > NumericTraits< LevelSetOutputRealType >::Zero )
        {
        if( Con( currentIndex, currentValue , update ) )
          {
          // CheckOut
          LevelSetOutputType oldValue = currentValue;
          LevelSetOutputType newValue = 1;

          this->m_InternalImage->SetPixel( currentIndex, newValue );

          m_EquationContainer->GetEquation( m_CurrentLevelSetId )->UpdatePixel(
                currentIndex,
                static_cast< LevelSetOutputRealType >( oldValue ),
                static_cast< LevelSetOutputRealType >( newValue )
                );

          list_out.insert(
                std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, newValue ) );

          LevelSetLayerIterator tempIt = nodeIt;
          ++nodeIt;
          list_out.erase( tempIt );

          erased = true;

          neighIt.SetLocation( currentIndex );

          for( typename NeighborhoodIteratorType::Iterator
              i = neighIt.Begin();
              !i.IsAtEnd(); ++i )
            {
            if ( i.Get() == -3 )
              {
              LevelSetInputType tempIndex =
                  neighIt.GetIndex( i.GetNeighborhoodOffset() );

              oldValue = -3;
              newValue = -1;

              list_in.insert(
                    std::pair< LevelSetInputType, LevelSetOutputType >( currentIndex, newValue ) );
              this->m_InternalImage->SetPixel( tempIndex, newValue );

              m_EquationContainer->GetEquation( m_CurrentLevelSetId )->UpdatePixel(
                    tempIndex,
                    static_cast< LevelSetOutputRealType >( oldValue ),
                    static_cast< LevelSetOutputRealType >( newValue )
                    );
              }
            }
          }
        }
      if( !erased )
        {
        ++nodeIt;
        }
      }
  }

  bool Con( const LevelSetInputType& iIdx,
            const LevelSetOutputType& iCurrentStatus,
            const LevelSetOutputRealType& iCurrentUpdate ) const
  {
    ZeroFluxNeumannBoundaryCondition< LabelImageType > sp_nbc;

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );

    NeighborhoodIteratorType neighIt( radius,
                                      m_InternalImage,
                                      m_InternalImage->GetLargestPossibleRegion() );

    neighIt.OverrideBoundaryCondition( &sp_nbc );

    typename NeighborhoodIteratorType::OffsetType sparse_offset;
    sparse_offset.Fill( 0 );

    for( unsigned int dim = 0; dim < ImageDimension; dim++ )
      {
      sparse_offset[dim] = -1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 1;
      neighIt.ActivateOffset( sparse_offset );
      sparse_offset[dim] = 0;
      }

    neighIt.SetLocation( iIdx );

    LevelSetOutputType opposite_status = ( iCurrentStatus == 1 ) ? -1 : 1;

    for( typename NeighborhoodIteratorType::Iterator
              i = neighIt.Begin();
          !i.IsAtEnd(); ++i )
      {
      if ( i.Get() == opposite_status )
        {
        LevelSetInputType tempIdx =
            neighIt.GetIndex( i.GetNeighborhoodOffset() );

        LevelSetOutputRealType neighborUpdate =
            m_EquationContainer->GetEquation( m_CurrentLevelSetId )->Evaluate( tempIdx );

        if ( neighborUpdate * iCurrentUpdate > NumericTraits< LevelSetOutputType >::Zero )
          {
          return true;
          }
        }
      }
   return false;
  }


private:
  UpdateShiSparseLevelSet( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkUpdateShiSparseLevelSet_h
