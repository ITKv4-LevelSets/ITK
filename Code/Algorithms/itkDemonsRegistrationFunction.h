/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkDemonsRegistrationFunction.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDemonsRegistrationFunction_h_
#define _itkDemonsRegistrationFunction_h_

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

namespace itk {

/**
 * \class DemonsRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration 
 * algorithm. It is used by DemonsRegistrationFilter to compute the 
 * output deformation field which will map a reference image onto a
 * a target image.
 *
 * Non-integer reference image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetReferenceInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the Reference image type, Target image type
 * and the deformation field type.
 *
 * \warning This filter assumes that the reference type, target type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa DemonsRegistrationFilter
 * \ingroup Functions
 */
template<class TReference, class TTarget, class TDeformationField>
class ITK_EXPORT DemonsRegistrationFunction : 
  public PDEDeformableRegistrationFunction< TReference,
    TTarget, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef DemonsRegistrationFunction    Self;
  typedef PDEDeformableRegistrationFunction< TReference,
    TTarget, TDeformationField >    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( DemonsRegistrationFunction, 
    PDEDeformableRegistrationFunction );

  /** Reference image type. */
  typedef typename Superclass::ReferenceType     ReferenceType;
  typedef typename Superclass::ReferencePointer  ReferencePointer;

  /** Target image type. */
  typedef typename Superclass::TargetType     TargetType;
  typedef typename Superclass::TargetPointer  TargetPointer;
  typedef typename TargetType::IndexType      IndexType;
  typedef typename TargetType::SizeType       SizeType;
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   
    DeformationFieldTypePointer;

  /** Inherit some enums from the superclass. */
  enum{ ImageDimension = Superclass::ImageDimension };

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType     PixelType;
  typedef typename Superclass::RadiusType    RadiusType;
  typedef typename Superclass::NeighborhoodType    NeighborhoodType;
  typedef typename Superclass::BoundaryNeighborhoodType    
                   BoundaryNeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType TimeStepType;

  /** Interpolator type. */
  typedef InterpolateImageFunction<ReferenceType> InterpolatorType;
  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorType::PointType       PointType;
  typedef LinearInterpolateImageFunction<ReferenceType>
    DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double,ImageDimension> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<TargetType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;

  /** Set the reference interpolator. */
  void SetReferenceInterpolator( InterpolatorType * ptr )
    { m_ReferenceInterpolator = ptr; }

  /** Get the reference interpolator. */
  InterpolatorPointer GetReferenceInterpolator()
    { return m_ReferenceInterpolator; }

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void *GlobalData) const
    { return m_TimeStep; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer() const
    {
    GlobalDataStruct *global = new GlobalDataStruct();
    return global;
    }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
    { delete (GlobalDataStruct *) GlobalData;  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  virtual PixelType  ComputeUpdate(const NeighborhoodType &neighborhood,
                     void *globalData,
                     const FloatOffsetType &offset = m_ZeroOffset) const;

  /** This method is called by a finite difference solver image filter at
   * each pixel that lies on the data set boundary. */
  virtual PixelType  ComputeUpdate(const BoundaryNeighborhoodType
                     &neighborhood, void *globalData,
                     const FloatOffsetType &offset = m_ZeroOffset) const;

protected:
  DemonsRegistrationFunction();
  ~DemonsRegistrationFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Target image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<TargetType> TargetNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the target image. */
  struct GlobalDataStruct
   {
   TargetNeighborhoodIteratorType   m_TargetIterator;
   };

private:
  DemonsRegistrationFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  /** Cache target image information. */
  const double *                  m_TargetSpacing;
  const double *                  m_TargetOrigin;

  /** Function to compute derivatives of the target image. */
  GradientCalculatorPointer       m_TargetGradientCalculator;

  /** Function to interpolate the reference image. */
  InterpolatorPointer             m_ReferenceInterpolator;

  /** The global timestep. */
  TimeStepType                    m_TimeStep;

  /** Constant used to guard against division by zero. */
  double                          m_EpsilonDenominator;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsRegistrationFunction.txx"
#endif

#endif
