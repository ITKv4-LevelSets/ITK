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

#ifndef ITKLEVELSETCONTAINERBASE_H
#define ITKLEVELSETCONTAINERBASE_H

#include <map>
#include "itkObject.h"

namespace itk
{
template< class TIdentifier,
          class TLevelSet >
class LevelSetContainerBase : public Object
{
public:
  typedef LevelSetContainerBase Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef Object Superclass;

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( LevelSetContainerBase, Object );

  typedef TIdentifier IdentifierType;

  typedef TLevelSet   LevelSetType;
  typedef typename LevelSetType::Pointer LevelSetPointer;
  typedef typename LevelSetType::InputType InputType;
  typedef typename LevelSetType::OutputType OutputType;
  typedef typename LevelSetType::GradientType GradientType;
  typedef typename LevelSetType::HessianType HessianType;

  LevelSetPointer GetLevelSet( const IdentifierType& iId ) const
    {
    typename std::map< IdentifierType, LevelSetPointer >::const_iterator
        it = m_Container.find( iId );

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
      typename std::map< IdentifierType, LevelSetPointer >::iterator
          it = m_Container.find( iId );

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
    typename std::map< IdentifierType, LevelSetPointer >::iterator
        it = m_Container.find( iId );

    if( it != m_Container.end() )
      {
      m_Container.erase( it );
      it->second = NULL;
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

  std::map< IdentifierType, LevelSetPointer > m_Container;


private:
  LevelSetContainerBase( const Self & );
  void operator = ( const Self & );
};

}
#endif // ITKLEVELSETCONTAINERBASE_H
