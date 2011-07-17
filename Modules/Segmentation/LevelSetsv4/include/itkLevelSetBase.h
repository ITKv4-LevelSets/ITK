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

#include "itkVector.h"
#include "itkMatrix.h"
#include "itkNumericTraits.h"
#include "itkDataObject.h"

namespace itk
{
template< class TInput,
          unsigned int VDimension,
          typename TOutput,
          class TDomain >
class LevelSetBase : public DataObject
{
public:
  typedef LevelSetBase               Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information */
  itkTypeMacro ( LevelSetBase, DataObject );

  typedef TInput                                           InputType;
  typedef TOutput                                          OutputType;
  typedef typename NumericTraits< OutputType >::RealType   OutputRealType;
  typedef Vector< OutputRealType, VDimension >             GradientType;
  typedef Matrix< OutputRealType, VDimension, VDimension > HessianType;

  typedef TDomain DomainType;

  /** Type used to define Regions */
  typedef long RegionType;

  virtual OutputType    Evaluate( const InputType& iP ) const = 0;
  virtual GradientType  EvaluateGradient( const InputType& iP ) const = 0;
  virtual HessianType   EvaluateHessian( const InputType& iP ) const = 0;

  /** Get the maximum number of regions that this data can be
   * separated into. */
  itkGetConstMacro(MaximumNumberOfRegions, RegionType);

  virtual void Initialize()
    {
    Superclass::Initialize();
    }

  /** Methods to manage streaming. */
  virtual void UpdateOutputInformation()
    {
    if( this->GetSource() )
      {
      this->GetSource()->UpdateOutputInformation();
      }

    // Now we should know what our largest possible region is. If our
    // requested region was not set yet, (or has been set to something
    // invalid - with no data in it ) then set it to the largest
    // possible region.
    if ( m_RequestedRegion == -1 && m_RequestedNumberOfRegions == 0 )
      {
      this->SetRequestedRegionToLargestPossibleRegion();
      }
    }

  virtual void SetRequestedRegionToLargestPossibleRegion()
    {
    m_RequestedNumberOfRegions     = 1;
    m_RequestedRegion           = 0;
    }

  virtual void CopyInformation(const DataObject *data)
    {
    const LevelSetBase *LevelSet = NULL;

    try
      {
      LevelSet = dynamic_cast< const LevelSetBase * >( data );
      }
    catch ( ... )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::LevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( LevelSetBase * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::LevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( LevelSetBase * ).name() );
      }

    m_MaximumNumberOfRegions = LevelSet->GetMaximumNumberOfRegions();

    m_NumberOfRegions = LevelSet->m_NumberOfRegions;
    m_RequestedNumberOfRegions = LevelSet->m_RequestedNumberOfRegions;
    m_BufferedRegion  = LevelSet->m_BufferedRegion;
    m_RequestedRegion = LevelSet->m_RequestedRegion;
    }

  virtual void Graft(const DataObject *data)
    {
    // Copy Meta Data
    this->CopyInformation(data);

    const Self *LevelSet = NULL;

    try
      {
      LevelSet = dynamic_cast< const Self * >( data );
      }
    catch ( ... )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::LevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }

    if ( !LevelSet )
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::LevelSetBase::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( Self * ).name() );
      }
    }

  virtual bool RequestedRegionIsOutsideOfTheBufferedRegion()
    {
    if ( m_RequestedRegion != m_BufferedRegion
         || m_RequestedNumberOfRegions != m_NumberOfRegions )
      {
      return true;
      }

    return false;
    }

  virtual bool VerifyRequestedRegion()
    {
    bool retval = true;

    // Are we asking for more regions than we can get?
    if ( m_RequestedNumberOfRegions > m_MaximumNumberOfRegions )
      {
      itkExceptionMacro(<< "Cannot break object into "
                        << m_RequestedNumberOfRegions << ". The limit is "
                        << m_MaximumNumberOfRegions);
      }

    if ( m_RequestedRegion >= m_RequestedNumberOfRegions
         || m_RequestedRegion < 0 )
      {
      itkExceptionMacro(<< "Invalid update region " << m_RequestedRegion
                        << ". Must be between 0 and "
                        << m_RequestedNumberOfRegions - 1);
      }

    return retval;
    }

  /** Set the requested region from this data object to match the requested
   * region of the data object passed in as a parameter.  This method
   * implements the API from DataObject. The data object parameter must be
   * castable to a PointSet. */
  virtual void SetRequestedRegion(DataObject *data)
    {
    Self *LevelSet = dynamic_cast< Self * >( data );

    if ( LevelSet )
      {
      // only copy the RequestedRegion if the parameter is another PointSet
      m_RequestedRegion = LevelSet->m_RequestedRegion;
      m_RequestedNumberOfRegions = LevelSet->m_RequestedNumberOfRegions;
      }
    }

  /** Set/Get the Requested region */
  virtual void SetRequestedRegion(const RegionType & region)
    {
    if ( m_RequestedRegion != region )
      {
      m_RequestedRegion = region;
      }
    }

  itkGetConstMacro(RequestedRegion, RegionType);

  /** Set/Get the Buffered region */
  virtual void SetBufferedRegion(const RegionType & region)
    {
    if ( m_BufferedRegion != region )
      {
      m_BufferedRegion = region;
      this->Modified();
      }
    }

  itkGetConstMacro(BufferedRegion, RegionType);


protected:
  LevelSetBase() : Superclass() {}
  virtual ~LevelSetBase() {}

  // If the RegionType is ITK_UNSTRUCTURED_REGION, then the following
  // variables represent the maximum number of region that the data
  // object can be broken into, which region out of how many is
  // currently in the buffered region, and the number of regions and
  // the specific region requested for the update. Data objects that
  // do not support any division of the data can simply leave the
  // MaximumNumberOfRegions as 1. The RequestedNumberOfRegions and
  // RequestedRegion are used to define the currently requested
  // region. The LargestPossibleRegion is always requested region = 0
  // and number of regions = 1;
  RegionType m_MaximumNumberOfRegions;
  RegionType m_NumberOfRegions;
  RegionType m_RequestedNumberOfRegions;
  RegionType m_BufferedRegion;
  RegionType m_RequestedRegion;

private:
  LevelSetBase( const Self& );
  void operator = ( const Self& );

};
}

#endif // __itkLevelSetBase_h
