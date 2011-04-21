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

#ifndef __itkLevelSetContainerBase_h
#define __itkLevelSetContainerBase_h

#include <map>
#include "itkObject.h"

namespace itk
{
template< class TIdentifier,
          class TLevelSet >
class LevelSetContainerBase : public Object
{
public:
  typedef LevelSetContainerBase      Self;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object                     Superclass;

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( LevelSetContainerBase, Object );

  typedef TIdentifier IdentifierType;

  typedef TLevelSet                           LevelSetType;
  typedef typename LevelSetType::Pointer      LevelSetPointer;
  typedef typename LevelSetType::InputType    InputType;
  typedef typename LevelSetType::OutputType   OutputType;
  typedef typename LevelSetType::GradientType GradientType;
  typedef typename LevelSetType::HessianType  HessianType;

  typedef std::map< IdentifierType, LevelSetPointer >    LevelSetContainerType;
  typedef typename LevelSetContainerType::const_iterator LevelSetContainerConstIteratorType;
  typedef typename LevelSetContainerType::iterator       LevelSetContainerIteratorType;

  LevelSetPointer GetLevelSet( const IdentifierType& iId ) const
    {
    LevelSetContainerConstIteratorType it = m_Container.find( iId );

    if( it != m_Container.end() )
      {
      return it->second;
      }
    else
      {
      return NULL;
      }
    }

  bool AddLevelSet( const IdentifierType& iId,
                    LevelSetPointer iLevelSet,
                    const bool& iForce = true )
    {
    if( iForce )
      {
      m_Container[iId] = iLevelSet;
      this->Modified();
      return true;
      }
    else
      {
      LevelSetContainerIteratorType it = m_Container.find( iId );

      if( it != m_Container.end() )
        {
        return false;
        }
      else
        {
        m_Container.insert(
              std::pair< IdentifierType, LevelSetPointer >( iId, iLevelSet ) );
        return true;
        }
      }
    }

  bool RemoveLevelSet( const IdentifierType& iId )
    {
    LevelSetContainerIteratorType it = m_Container.find( iId );

    if( it != m_Container.end() )
      {
      it->second = NULL;
      m_Container.erase( it );

      this->Modified();

      return true;
      }
    else
      {
      return false;
      }
    }

protected:
  LevelSetContainerBase() {}
  ~LevelSetContainerBase() {}

  LevelSetContainerType m_Container;

private:
  LevelSetContainerBase( const Self & );
  void operator = ( const Self & );
};

}
#endif // __itkLevelSetContainerBase_h
