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

#ifndef __itkLevelSetQuadEdgeMeshBase_h
#define __itkLevelSetQuadEdgeMeshBase_h

#include "itkLevelSetBase.h"

namespace itk
{
template< class TMesh >
class LevelSetQuadEdgeMeshBase :
    public LevelSetBase<
      typename TMesh::PointIdentifier,
      TMesh::PointDimension,
      typename TMesh::PixelType,
      TMesh >
{
public:
  typedef TMesh                   MeshType;
  typedef typename TMesh::Pointer MeshPointer;

  typedef LevelSetQuadEdgeMeshBase        Self;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;
  typedef LevelSetBase< typename MeshType::PointIdentifier,
    MeshType::PointDimension,
    typename MeshType::PixelType,
    MeshType                    >         Superclass;

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( LevelSetQuadEdgeMeshBase, LevelSetBase );

  typedef typename Superclass::InputType    InputType;
  typedef typename Superclass::OutputType   OutputType;
  typedef typename Superclass::GradientType GradientType;
  typedef typename Superclass::HessianType  HessianType;

  itkSetObjectMacro( Mesh, MeshType );
  itkGetObjectMacro( Mesh, MeshType );

  OutputType Evaluate( const InputType& iP ) const
    {
    OutputType oValue = 0.;
    m_Mesh->GetPointData( iP, &oValue );
    return oValue;
    }

  GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }

  struct LevelSetDataType
    {
    OutputType Value;
    GradientType Gradient;
    HessianType Hessian;
    };

  LevelSetDataType* GetAllData( const InputType& iP ) const
    {
    LevelSetDataType* oData = new LevelSetDataType;
    oData->Value = m_Mesh->GetPixel( iP );

    return oData;
    }

protected:
  LevelSetQuadEdgeMeshBase()
    {
    m_Mesh = MeshType::New();
    }
  ~LevelSetQuadEdgeMeshBase() {}

  MeshPointer m_Mesh;

private:
  LevelSetQuadEdgeMeshBase( const Self& );
  void operator = ( const Self& );
};
}

#endif // __itkLevelSetQuadEdgeMeshBase_h
