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
#include "itkHeavisideStepFunctionBase.h"
#include "itkLevelSetDomainMapImageFilter.h"

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

  /** typedefs related to the type of level set*/
  typedef TLevelSet                             LevelSetType;
  typedef typename LevelSetType::Pointer        LevelSetPointer;
  typedef typename LevelSetType::ImageType      LevelSetImageType;
  typedef typename LevelSetType::InputType      InputIndexType;
  typedef typename LevelSetType::OutputType     OutputType;
  typedef typename LevelSetType::OutputRealType OutputRealType;
  typedef typename LevelSetType::GradientType   GradientType;
  typedef typename LevelSetType::HessianType    HessianType;

  /** IdentifierType */
  typedef TIdentifier IdentifierType;

  typedef std::map< IdentifierType, LevelSetPointer >    LevelSetContainerType;
  typedef typename LevelSetContainerType::const_iterator LevelSetContainerConstIteratorType;
  typedef typename LevelSetContainerType::iterator       LevelSetContainerIteratorType;

  typedef HeavisideStepFunctionBase< OutputRealType, OutputRealType >
                                            HeavisideType;
  typedef typename HeavisideType::Pointer   HeavisidePointer;

  itkStaticConstMacro ( ImageDimension, unsigned int,
                       LevelSetImageType::ImageDimension );

  typedef std::list< IdentifierType >                    IdListType;
  typedef typename IdListType::iterator                  IdListIterator;
  typedef Image< IdListType, ImageDimension >            IdListImageType;
  typedef Image< short, ImageDimension >                 CacheImageType;
  typedef LevelSetDomainMapImageFilter< IdListImageType, CacheImageType >
                                                         DomainMapImageFilterType;

  typedef typename DomainMapImageFilterType::Pointer         DomainMapImageFilterPointer;
  typedef typename DomainMapImageFilterType::NounToBeDefined NounToBeDefined;

  typedef typename std::map< IdentifierType, NounToBeDefined >::iterator DomainIteratorType;

  LevelSetContainerIteratorType Begin()
    {
    return m_Container.begin();
    }

  LevelSetContainerConstIteratorType Begin() const
    {
    return m_Container.begin();
    }

  LevelSetContainerIteratorType End()
    {
    return m_Container.end();
    }

  LevelSetContainerConstIteratorType End() const
    {
    return m_Container.end();
    }

  /** \brief Get the level set function given its id
    \param[in] iId
    \return the level set function if it is in the container, else NULL.
  */
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

  /** \brief Add one level set function given its id.

    \param[in] iId id of the level set function
    \param[in] iLevelSet the level set function to be added
    \param[in] iForce if iForce is true (default) the level set function will be
    added to the container even if there is already one with the same id.

    \return true if the level set has been added.
  */
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

  /// \warning why
  itkSetObjectMacro( Heaviside, HeavisideType );
  itkGetObjectMacro( Heaviside, HeavisideType );

  // set the domain map image filter
  itkSetObjectMacro( DomainMapFilter, DomainMapImageFilterType );
  itkGetObjectMacro( DomainMapFilter, DomainMapImageFilterType );

protected:
  /** \brief Default Constructor */
  LevelSetContainerBase() {}

  /** \brief Default Destructor */
  ~LevelSetContainerBase() {}

  HeavisidePointer            m_Heaviside;
  DomainMapImageFilterPointer m_DomainMapFilter;
  LevelSetContainerType       m_Container;

private:
  LevelSetContainerBase( const Self & );
  void operator = ( const Self & );
};

}
#endif // __itkLevelSetContainerBase_h
