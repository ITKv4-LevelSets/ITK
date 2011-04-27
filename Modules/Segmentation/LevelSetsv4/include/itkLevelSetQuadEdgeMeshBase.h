/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-11

 Copyright (c) 2009-11, President and Fellows of Harvard College.
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

#ifndef __itkLevelSetQuadEdgeMesh_h
#define __itkLevelSetQuadEdgeMesh_h

#include "itkLevelSetBase.h"

namespace itk
{
template< class TMesh >
class LevelSetQuadEdgeMeshBase :
    public LevelSetBase<
      typename TMesh::PointIdentifier,
      TMesh::PointDimension,
      typename TMesh::PixelType >
{
public:
  typedef TMesh MeshType;
  typedef typename TMesh::Pointer MeshPointer;

  typedef LevelSetQuadEdgeMeshBase Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef LevelSetBase< typename MeshType::PointIdentifier,
    MeshType::PointDimension,
    typename MeshType::PixelType > Superclass;

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

#endif // __itkLevelSetQuadEdgeMesh_h