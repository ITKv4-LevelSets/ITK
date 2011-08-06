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

#ifndef __itkLevelSetEquationContainerBase_h
#define __itkLevelSetEquationContainerBase_h

#include "itkObject.h"

namespace itk
{
template< class TTermContainer >
class LevelSetEquationContainerBase : public Object
{
public:
  typedef LevelSetEquationContainerBase Self;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;
  typedef Object                        Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationContainerBase, Object );

  typedef TTermContainer                            TermContainerType;
  typedef typename TermContainerType::Pointer       TermContainerPointer;

  typedef typename TermContainerType::InputImageType    InputImageType;
  typedef typename TermContainerType::InputImagePointer InputImagePointer;

  typedef typename TermContainerType::LevelSetOutputRealType  LevelSetOutputRealType;
  typedef typename TermContainerType::LevelSetInputIndexType  LevelSetInputIndexType;

  typedef typename TermContainerType::LevelSetIdentifierType    LevelSetIdentifierType;
  typedef typename TermContainerType::LevelSetContainerType     LevelSetContainerType;
  typedef typename TermContainerType::LevelSetContainerPointer  LevelSetContainerPointer;

  LevelSetContainerType* GetLevelSetContainer();

  void AddEquation( const LevelSetIdentifierType& iId,
                    TermContainerPointer iEquation );

  TermContainerType* GetEquation( const LevelSetIdentifierType& iId );

  void Update();

  void UpdatePixel( const LevelSetInputIndexType& iP,
                    const LevelSetOutputRealType & oldValue,
                    const LevelSetOutputRealType & newValue );

  void InitializeParameters();

  LevelSetOutputRealType GetCFLContribution();

  itkSetObjectMacro( Input, InputImageType );
  itkGetObjectMacro( Input, InputImageType );

protected:

  LevelSetEquationContainerBase();
  ~LevelSetEquationContainerBase();

  typedef std::map< LevelSetIdentifierType, TermContainerPointer >  MapContainerType;
  typedef typename MapContainerType::iterator                       MapContainerIterator;

  MapContainerType  m_Container;
  InputImagePointer m_Input;


private:
  LevelSetEquationContainerBase( const Self& );
  void operator = ( const Self& );

};
}

#include "itkLevelSetEquationContainerBase.hxx"
#endif // __itkLevelSetEquationContainerBase_h
