/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkPointSetToListAdaptor.txx
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPointSetToListAdaptor_txx
#define _itkPointSetToListAdaptor_txx

#include "itkPointSetToListAdaptor.h"

namespace itk{ 
namespace Statistics{

template < class TPointSet >
void
PointSetToListAdaptor< TPointSet >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "PointSet: " << m_PointSet << std::endl;
}

template < class TPointSet >
void
PointSetToListAdaptor< TPointSet >
::SetPointSet(TPointSet* pointSet) 
{ 
  m_PointSet = pointSet ; 
  m_PointsContainer = pointSet->GetPoints() ;
}

template < class TPointSet >
TPointSet*
PointSetToListAdaptor< TPointSet >
::GetPointSet() 
{
  return m_PointSet.GetPointer() ; 
}  

/** returns the number of measurement vectors in this container*/
template < class TPointSet >
unsigned int
PointSetToListAdaptor< TPointSet >
::Size() const
{
  return m_PointsContainer->Size() ;
}

template < class TPointSet >
inline typename PointSetToListAdaptor< TPointSet >::MeasurementVectorType
PointSetToListAdaptor< TPointSet >
::GetMeasurementVector(const InstanceIdentifier &id)
{
  m_PointSet->GetPoint(id, &m_TempPoint) ;
  return (MeasurementVectorType&) m_TempPoint ;
}

template < class TPointSet >
inline typename PointSetToListAdaptor< TPointSet >::FrequencyType
PointSetToListAdaptor< TPointSet >
::GetFrequency(const InstanceIdentifier &) const 
{
  return 1 ;
}

template < class TPointSet >
typename PointSetToListAdaptor< TPointSet >::FrequencyType
PointSetToListAdaptor< TPointSet >
::GetTotalFrequency() const
{ 
  return this->Size() ; 
}

} // end of namespace Statistics 
} // end of namespace itk

#endif



