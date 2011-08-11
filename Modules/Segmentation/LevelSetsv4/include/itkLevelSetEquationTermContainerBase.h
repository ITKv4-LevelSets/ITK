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

#ifndef __itkLevelSetEquationTermContainerBase_h
#define __itkLevelSetEquationTermContainerBase_h

#include "itkLevelSetEquationTermBase.h"
#include "itkObject.h"

#include "itksys/hash_map.hxx"

#include <map>
#include <string>

namespace itk
{
template< class TInputImage,
          class TLevelSetContainer >
class LevelSetEquationTermContainerBase : public Object
{
public:
  typedef LevelSetEquationTermContainerBase Self;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;
  typedef Object                            Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  /** Run-time type information */
  itkTypeMacro( LevelSetEquationTermContainerBase,
                Object );

  typedef unsigned int                      TermIdType;

  typedef TInputImage                       InputImageType;
  typedef typename InputImageType::Pointer  InputImagePointer;

  typedef TLevelSetContainer                              LevelSetContainerType;
  typedef typename LevelSetContainerType::Pointer         LevelSetContainerPointer;
  typedef typename LevelSetContainerType::LevelSetIdentifierType
                                                          LevelSetIdentifierType;
//   typedef typename LevelSetContainerType::OutputPixelType LevelSetOutputPixelType;
  typedef typename LevelSetContainerType::OutputType      LevelSetOutputPixelType;
  typedef typename LevelSetContainerType::OutputRealType  LevelSetOutputRealType;
  typedef typename LevelSetContainerType::InputIndexType  LevelSetInputIndexType;
  typedef typename LevelSetContainerType::GradientType    GradientType;
  typedef typename LevelSetContainerType::HessianType     HessianType;

  typedef LevelSetEquationTermBase< InputImageType, LevelSetContainerType >
                                                                       TermType;
  typedef typename TermType::Pointer                                   TermPointer;

  itkSetObjectMacro( Input, InputImageType );
  itkGetObjectMacro( Input, InputImageType );

  void PushTerm( TermType* iTerm );

  void AddTerm( const TermIdType& iId, TermType* iTerm );

  TermType* GetTerm( const TermIdType& iId );
  TermType* GetTerm( const std::string& iName );

  void Initialize( const LevelSetInputIndexType& iP );

  void UpdatePixel( const LevelSetInputIndexType& iP,
                    const LevelSetOutputRealType & oldValue,
                    const LevelSetOutputRealType & newValue );

  void InitializeParameters();

  LevelSetOutputRealType Evaluate( const LevelSetInputIndexType& iP );

  void Update();

  LevelSetOutputRealType ComputeCFLContribution() const;

protected:
  LevelSetEquationTermContainerBase();

  virtual ~LevelSetEquationTermContainerBase();

  InputImagePointer     m_Input;

  struct hash_string
  {
    size_t operator()( const std::string& x ) const
    {
      return itksys::hash< const char* >()( x.c_str() );
    }
  };

  typedef itksys::hash_map< std::string,
                            TermPointer,
                            hash_string >                   HashMapStringTermContainerType;

  HashMapStringTermContainerType m_NameContainer;

  typedef std::map< TermIdType, TermPointer >           MapTermContainerType;
  typedef typename MapTermContainerType::iterator       MapTermContainerIteratorType;
  typedef typename MapTermContainerType::const_iterator MapTermContainerConstIteratorType;

  MapTermContainerType  m_Container;

  typedef std::map< TermIdType, LevelSetOutputRealType >  MapCFLContainerType;
  typedef typename MapCFLContainerType::iterator          MapCFLContainerIterator;
  typedef typename MapCFLContainerType::const_iterator    MapCFLContainerConstIterator;

  MapCFLContainerType   m_TermContribution;

private:
  LevelSetEquationTermContainerBase( const Self& );
  void operator = ( const Self& );
};

}
#include "itkLevelSetEquationTermContainerBase.hxx"
#endif // __itkLevelSetEquationTermContainerBase_h
