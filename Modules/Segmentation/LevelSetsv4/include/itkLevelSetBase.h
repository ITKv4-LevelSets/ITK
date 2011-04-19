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

#ifndef __itkLevelSetBase_h
#define __itkLevelSetBase_h

#include "itkLevelSetBase.h"
#include "itkVector.h"
#include "itkMatrix.h"
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

  typedef TInput                                       InputType;
  typedef TOutput                                      OutputType;
  typedef Vector< OutputType, VDimension >             GradientType;
  typedef Matrix< OutputType, VDimension, VDimension > HessianType;

  virtual OutputType    Evaluate( const InputType& iP ) const = 0;
  virtual GradientType  EvaluateGradient( const InputType& iP ) const = 0;
  virtual HessianType   EvaluateHessian( const InputType& iP ) const = 0;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */

  itkConceptMacro( DoubleConvertible,
                    ( Concept::Convertible< double, OutputType > ) );

  /** End concept checking */
#endif // ITK_USE_CONCEPT_CHECKING

protected:
  LevelSetBase() {}
  virtual ~LevelSetBase() {}

private:
  LevelSetBase( const Self& );
  void operator = ( const Self& );

};
}

#endif // __itkLevelSetBase_h
