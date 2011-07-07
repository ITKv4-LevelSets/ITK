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

  typedef TInputImage                      InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;

  typedef TLevelSetContainer                              LevelSetContainerType;
  typedef typename LevelSetContainerType::Pointer         LevelSetContainerPointer;
  typedef typename LevelSetContainerType::IdentifierType  LevelSetIdentifierType;
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

  void AddTerm( const unsigned int& iId, TermPointer iTerm )
    {
    if ( iTerm.IsNotNull() )
      {
      if( iTerm->GetInput() == NULL )
        {
        if( m_Input.IsNotNull() )
          {
          iTerm->SetInput( m_Input );
          }
        else
          {
          itkGenericExceptionMacro( <<"m_Input and iTerm->GetInput are NULL" );
          }
        }
      m_Container[iId] = iTerm;
      m_TermContribution[iId] = NumericTraits< LevelSetOutputPixelType >::Zero;

      this->Modified();
      }
    else
      {
      itkGenericExceptionMacro( <<"Term supplied is null" );
      }
    }

  TermPointer GetTerm( const unsigned int& iId )
    {
    typename std::map< unsigned int, TermPointer >::iterator
        it = m_Container.find( iId );

    if( it != m_Container.end() )
      {
      return it->second;
      }
    else
      {
      itkGenericExceptionMacro( <<"this term does not exist" );
      return TermPointer();
      }
    }


  void Initialize( const LevelSetInputIndexType& iP )
  {
    typename std::map< unsigned int, TermPointer >::iterator
      term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
      term_end = m_Container.end();

    while( term_it != term_end )
      {
      ( term_it->second )->Initialize( iP );
      ++term_it;
      }
  }

  void UpdatePixel( LevelSetInputIndexType& iP, LevelSetOutputRealType & oldValue, LevelSetOutputRealType & newValue )
  {
    typename std::map< unsigned int, TermPointer >::iterator
      term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
      term_end = m_Container.end();

    while( term_it != term_end )
    {
      ( term_it->second )->UpdatePixel( iP, oldValue, newValue );
      ++term_it;
    }
  }

  void InitializeParameters()
  {
    typename std::map< unsigned int, TermPointer >::iterator
    term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
    term_end = m_Container.end();

    while( term_it != term_end )
    {
      ( term_it->second )->InitializeParameters();
      ++term_it;
    }
  }

  LevelSetOutputRealType Evaluate( const LevelSetInputIndexType& iP )
    {
    typename std::map< unsigned int, TermPointer >::iterator
        term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
        term_end = m_Container.end();

    LevelSetOutputRealType oValue = NumericTraits< LevelSetOutputRealType >::Zero;

    while( term_it != term_end )
      {
      LevelSetOutputRealType temp_val = ( term_it->second )->Evaluate( iP );

      m_TermContribution[ term_it->first ] =
          vnl_math_max( temp_val, m_TermContribution[ term_it->first ] );

      oValue += temp_val;
      ++term_it;
      }

    return oValue;
    }

  void Update()
    {
    typename std::map< unsigned int, TermPointer >::iterator
        term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
        term_end = m_Container.end();

    while( term_it != term_end )
      {
      ( term_it->second )->Update();
      m_TermContribution[ term_it->first ] =
          NumericTraits< LevelSetOutputPixelType >::Zero;
      ++term_it;
      }
    }

  LevelSetOutputRealType GetCFLContribution()
    {
    typename std::map< unsigned int, TermPointer >::iterator
        term_it = m_Container.begin();
    typename std::map< unsigned int, TermPointer >::iterator
        term_end = m_Container.end();

    LevelSetOutputRealType oValue = NumericTraits< LevelSetOutputRealType >::Zero;

    while( term_it != term_end )
      {
      LevelSetOutputRealType cfl = ( term_it->second )->GetCFLContribution();
      if( cfl == NumericTraits< LevelSetOutputRealType >::Zero )
        {
        cfl = m_TermContribution[ term_it->first ];
        }

      oValue += cfl;
      ++term_it;
      }

    return oValue;
    }

protected:
  LevelSetEquationTermContainerBase() : Superclass(), m_Input( NULL ) {}

  virtual ~LevelSetEquationTermContainerBase() {}

  std::map< unsigned int, TermPointer >              m_Container;
  std::map< unsigned int, LevelSetOutputRealType >   m_TermContribution;

  InputImagePointer                             m_Input;


private:
  LevelSetEquationTermContainerBase( const Self& );
  void operator = ( const Self& );
};

}
#endif // __itkLevelSetEquationTermContainerBase_h
