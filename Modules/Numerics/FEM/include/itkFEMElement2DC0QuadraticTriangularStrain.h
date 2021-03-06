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
#ifndef __itkFEMElement2DC0QuadraticTriangularStrain_h
#define __itkFEMElement2DC0QuadraticTriangularStrain_h

#include "itkFEMElement2DC0QuadraticTriangular.h"
#include "itkFEMElement2DStrain.h"

namespace itk {
namespace fem {

/**
 * \class Element2DC0QuadraticTriangularStrain
 * \brief 3-noded finite element class in 2D space for linear elasticity problem.
 *
 * This element is combined from Element2DC0LinearTriangular and Element2DStrain.
 * \ingroup ITK-FEM
 */
class Element2DC0QuadraticTriangularStrain : public Element2DStrain<Element2DC0QuadraticTriangular>
{
FEM_CLASS(Element2DC0QuadraticTriangularStrain,Element2DStrain<Element2DC0QuadraticTriangular>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element2DC0QuadraticTriangularStrain();

  /**
   * Construct an element by specifying pointers to
   * 3 points and a material.
   */
  Element2DC0QuadraticTriangularStrain(
      NodeIDType n1_,
      NodeIDType n2_,
      NodeIDType n3_,
      NodeIDType n4_,
      NodeIDType n5_,
      NodeIDType n6_,
      Material::ConstPointer p_ );

}; // class Element2DC0QuadraticTriangularStrain

FEM_CLASS_INIT(Element2DC0QuadraticTriangularStrain)

}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement2DC0QuadraticTriangularStrain_h
