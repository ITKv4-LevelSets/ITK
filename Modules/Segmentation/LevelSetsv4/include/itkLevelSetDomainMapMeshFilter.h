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

#ifndef __itkLevelSetDomainMapMeshFilter_h
#define __itkLevelSetDomainMapMeshFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkIntTypes.h"

namespace itk
{
/**
  \class LevelSetDomainMapMeshFilter
  \tparam TInputMesh Image where the pixel type is a container of ids
  \tparam TOutputMesh Image where the pixel type is an integer to split the region
*/
template < class TInputMesh, class TOutputMesh >
class ITK_EXPORT LevelSetDomainMapMeshFilter : public
    QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
  public:
    typedef LevelSetDomainMapMeshFilter                      Self;
    typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh,
                                              TOutputMesh >  Superclass;
    typedef SmartPointer< Self >                             Pointer;
    typedef SmartPointer< const Self >                       ConstPointer;

    itkStaticConstMacro ( PointDimension, unsigned int,
                          TInputMesh::PointDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( LevelSetDomainMapMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

    typedef TInputMesh                            InputMeshType;
    typedef typename InputMeshType::ConstPointer  InputMeshConstPointer;
    typedef typename InputMeshType::PixelType     InputMeshPixelType;
    typedef typename InputMeshType::PointDataContainerPointer
                                                  InputPointDataContainerPointer;
    typedef typename InputMeshType::PointDataContainerIterator
                                                  InputPointDataContainerIterator;
    typedef typename InputMeshType::PointIdentifier
                                                  InputPointIdentifier;

    typedef TOutputMesh                           OutputMeshType;
    typedef typename OutputMeshType::Pointer      OutputMeshPointer;
    typedef typename OutputMeshType::PixelType    OutputMeshPixelType;
    typedef typename OutputMeshType::PointDataContainerPointer
                                                  OutputPointDataContainerPointer;
    typedef typename OutputMeshType::PointDataContainerIterator
                                                  OutputPointDataContainerIterator;
    typedef typename OutputMeshType::PointIdentifier
                                                  OutputPointIdentifier;

    struct NounToBeDefined // ~ kind of cache to speed up computations
      {
      NounToBeDefined() {}

      NounToBeDefined( const IdentifierType& id,
                       const InputPointIdentifier& iStart,
                       const InputPointIdentifier& iEnd,
                       const InputMeshPixelType& iList )
        : m_Id( id ), m_Start( iStart ), m_End( iEnd ), m_List( iList ) {}

      IdentifierType m_Id;
      InputPointIdentifier m_Start;
      InputPointIdentifier m_End;
      InputMeshPixelType m_List;
      };

    std::map< IdentifierType, NounToBeDefined > m_LevelSetMap;

  protected:
    LevelSetDomainMapMeshFilter();
    ~LevelSetDomainMapMeshFilter() {}

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

//    void ConsistencyCheck( bool& subRegionConsistent, InputImageRegionType& subRegion );
    void GenerateData();

  private:
    LevelSetDomainMapMeshFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkLevelSetDomainMapMeshFilter.txx"
#endif
