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
#ifndef __itkResourceProbesCollectorBase_h
#define __itkResourceProbesCollectorBase_h


#include "itkResourceProbe.h"
#include "itkMemoryUsageObserver.h"

namespace itk
{
/** \class ResourceProbesCollectorBase
 *  \brief Class for aggregating a set of probes.
 *
 *  This class defines a set of ResourceProbes and assign names to them.
 *  The user can start and stop each one of the probes by addressing them by name.
 *
 *  \sa ResourceProbe
 *
 * \ingroup ITK-Common
 */
template< class TProbe >
class ITK_EXPORT ResourceProbesCollectorBase
{
public:
  typedef std::string                IdType;
  typedef std::map< IdType, TProbe > MapType;

  /** destructor */
  virtual ~ResourceProbesCollectorBase();

  /** Start a probe with a particular name. If the time probe does not
   * exist, it will be created */
  virtual void Start(const char *name);

  /** Stop a time probe identified with a name */
  virtual void Stop(const char *name);

  /** Report the summary of results from the probes */
  virtual void Report(std::ostream & os = std::cout) const;

  /** Destroy the set of probes. New probes can be created after invoking this
    method. */
  virtual void Clear(void);

protected:
  MapType m_Probes;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResourceProbesCollectorBase.txx"
#endif

#endif //__itkResourceProbesCollectorBase_h
