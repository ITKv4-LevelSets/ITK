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

#ifndef __itkLevelSetContainerBase_hxx
#define __itkLevelSetContainerBase_hxx

#include "itkLevelSetContainerBase.h"

namespace itk
{
template< class TIdentifier, class TLevelSet >
LevelSetContainerBase< TIdentifier, TLevelSet >
::LevelSetContainerBase() : Superclass(),
  m_Heaviside( NULL ), m_DomainMapFilter( NULL )
{}

template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::Iterator
LevelSetContainerBase< TIdentifier, TLevelSet >::Begin()
  {
  return Iterator( m_Container.begin() );
  }

template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::ConstIterator
LevelSetContainerBase< TIdentifier, TLevelSet >::Begin() const
  {
  return ConstIterator( m_Container.begin() );
  }

template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::Iterator
LevelSetContainerBase< TIdentifier, TLevelSet >::End()
  {
  return Iterator( m_Container.end() );
  }

template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::ConstIterator
LevelSetContainerBase< TIdentifier, TLevelSet >::End() const
  {
  return ConstIterator( m_Container.end() );
  }


template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::LevelSetIdentifierType
LevelSetContainerBase< TIdentifier, TLevelSet >::Size() const
  {
  return static_cast< LevelSetIdentifierType >( m_Container.size() );
  }

/** \brief Get the level set function given its id
  \param[in] iId
  \return the level set function if it is in the container, else NULL.
*/
template< class TIdentifier, class TLevelSet >
typename LevelSetContainerBase< TIdentifier, TLevelSet >::LevelSetPointer
LevelSetContainerBase< TIdentifier, TLevelSet >
::GetLevelSet( const LevelSetIdentifierType& iId ) const
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

/** \brief Add one level set function given its id.

  \param[in] iId id of the level set function
  \param[in] iLevelSet the level set function to be added
  \param[in] iForce if iForce is true (default) the level set function will be
  added to the container even if there is already one with the same id.

  \return true if the level set has been added.
*/
template< class TIdentifier, class TLevelSet >
bool LevelSetContainerBase< TIdentifier, TLevelSet >
::AddLevelSet( const LevelSetIdentifierType& iId,
               LevelSetPointer iLevelSet,
               const bool& iForce )
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
            std::pair< LevelSetIdentifierType, LevelSetPointer >( iId, iLevelSet ) );
      this->Modified();
      return true;
      }
    }
  }

/** \brief Remove one level set function given its id.
  \param[in] iId id of the level set function to be removed
  \return true if it has been removed, false if the id was not present in the
  container.
*/
template< class TIdentifier, class TLevelSet >
bool
LevelSetContainerBase< TIdentifier, TLevelSet >
::RemoveLevelSet( const LevelSetIdentifierType& iId )
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
}
#endif // __itkLevelSetContainerBase_hxx
