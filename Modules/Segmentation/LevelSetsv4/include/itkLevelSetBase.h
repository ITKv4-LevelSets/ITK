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

#ifndef __itkLevelSetBase_h
#define __itkLevelSetBase_h

#include "itkLevelSetBase.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkNumericTraits.h"
#include "itkConceptChecking.h"

namespace itk
{
template< class TInput,
          unsigned int VDimension,
          typename TOutput >
class LevelSetBase : public Object
{
public:
  typedef LevelSetBase               Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information */
  itkTypeMacro ( LevelSetBase, Object );

  typedef TInput                                           InputType;
  typedef TOutput                                          OutputType;
  typedef typename NumericTraits< OutputType >::RealType   OutputRealType;
  typedef Vector< OutputRealType, VDimension >             GradientType;
  typedef Matrix< OutputRealType, VDimension, VDimension > HessianType;

  virtual OutputType    Evaluate( const InputType& iP ) const = 0;
  virtual GradientType  EvaluateGradient( const InputType& iP ) const = 0;
  virtual HessianType   EvaluateHessian( const InputType& iP ) const = 0;

protected:
  LevelSetBase() {}
  virtual ~LevelSetBase() {}

private:
  LevelSetBase( const Self& );
  void operator = ( const Self& );

};
}

#endif // __itkLevelSetBase_h
