/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkImageToImageAffineMutualInformationGradientDescentRegistration.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToImageAffineMutualInformationGradientDescentRegistration_h
#define __itkImageToImageAffineMutualInformationGradientDescentRegistration_h

#include "itkRegistrationMethod.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkGradientDescentOptimizer.h"
#include "itkImageMapper.h"
#include "itkAffineTransform.h"
#include "itkPoint.h"

namespace itk
{

/**
 *  Traits class that defines the different types to be
 *  used by this registration method
 */
template <class TReference, class TTarget>
class ITK_EXPORT 
  ImageToImageAffineMutualInformationGradientDescentRegistrationTraits 
{
public:
  /**  Type of the reference. */
  typedef TReference  ReferenceType;

  /**  Type of the target. */
  typedef TTarget TargetType;

  /** Image dimensions. */
  enum { ImageDimension = ReferenceType::ImageDimension }; 

  /**  Type of the transformation. */
  typedef AffineTransform< double, ImageDimension >  TransformationType;

  /**  Type of the parameters. */
  typedef typename TransformationType::ParametersType  ParametersType;

  /** Parameters dimensions. */
  enum { ParametersDimension = TransformationType::ParametersDimension }; 

  /**  Type of the mapper. */
  typedef ImageMapper<ReferenceType,
    TransformationType> MapperType;

  /**  Type of the metric. */
  typedef MutualInformationImageToImageMetric<TargetType, 
    MapperType> MetricType;

  /**  Type of the optimizer.  */
  typedef GradientDescentOptimizer<MetricType> OptimizerType;
};

/** \class ImageToImageAffineMutualInformationGradientDescentRegistration
 * \brief Register two images using mutual information.
 *
 * ImageToImageAffineMutualInformationGradientDescentRegistration performs an
 * affine registration
 * between a target and reference image using mutual information.
 * It uses the optimization method of Viola and Wells to find the
 * best affine transform to register the reference image onto the target
 * image. The two images are defined via methods SetTarget()
 * and SetReference().
 *
 * The mutual information value and its derivatives are estimated
 * using spatial sampling. The performance
 * of the registration depends on good choices of the parameters
 * used to estimate the mutual information. Refer to the documentation
 * for MutualInformationImageToImageMetric for details on these
 * parameters and how to set them.
 *
 * The registration uses a simple stochastic gradient ascent scheme. Steps
 * are repeatedly taken that are proportional to the approximate
 * deriviative of the mutual information with respect to the affine
 * transform parameters. The stepsize is governed by the LearningRate
 * parameter. The default is 1.0. The LearningRate is set via the method
 * SetLearningRate().
 *
 * Since the parameters of the linear part is different in magnitude
 * to the parameters in the offset part, scaling is required
 * to improve convergence. The scaling can set via the optimizer.
 *
 * NB: In the Viola and Wells paper, the scaling is specified by
 * using different learning rates for the linear and offset part.
 * The following formula translate their scaling parameters to
 * those used in this framework:
 *
 * LearningRate = lambda_R
 * TranslationScale = sqrt( lambda_T / lambda_R );
 *
 * In the optimizer's scale transform set the scaling for
 * all the translation parameters to TranslationScale^{-2}.
 * Set the scale for all other parameters to 1.0.
 *
 * Note: the optimization performance can be improved by 
 * setting the image origin to center of mass of the image.
 *
 * Implementaton of this class is based on:
 * Viola, P. and Wells III, W. (1997).
 * "Alignment by Maximization of Mutual Information"
 * International Journal of Computer Vision, 24(2):137-154
 *
 * This class is templated over the reference image type and the
 * the target image type.
 *
 * \sa MutualInformationImageToImageMetric
 * \sa AffineTransform
 *
 *  \ingroup AffineImageRegistration
 */
template <class TReference, class TTarget>
class ITK_EXPORT 
ImageToImageAffineMutualInformationGradientDescentRegistration 
: public RegistrationMethod< 
   ImageToImageAffineMutualInformationGradientDescentRegistrationTraits
   <TReference, TTarget>  >
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageAffineMutualInformationGradientDescentRegistration Self;
  typedef RegistrationMethod< 
          ImageToImageAffineMutualInformationGradientDescentRegistrationTraits<
          TReference,TTarget> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageAffineMutualInformationGradientDescentRegistration,
               RegistrationMethod);

  /**  Type of the reference. */
  typedef TReference  ReferenceType;

  /**  Type of the target. */
  typedef TTarget TargetType;

  /** Image dimensions. */
  enum {ImageDimension = ReferenceType::ImageDimension};

  /** Parmeters dimensions. */
  enum { ParametersDimension = ImageDimension * ( ImageDimension + 1 ) };

  /** Parameters type. */
  typedef typename Superclass::ParametersType   ParametersType;

  /**  Type of the transformation. */
  typedef typename Superclass::TransformationType    TransformationType;

  /** Point type. */
  typedef typename TransformationType::InputPointType  PointType;

  /**  Type of the mapper. */
  typedef typename Superclass::MapperType   MapperType;

  /**  Type of the metric. */
  typedef typename Superclass::MetricType   MetricType;

  /**  Type of the optimizer.  */
  typedef typename  Superclass::OptimizerType  OptimizerType;

  /**  Pointer type for the optimizer.  */
  typedef typename OptimizerType::Pointer  OptimizerPointer;

  /**  Pointer type for the reference. */
  typedef typename ReferenceType::ConstPointer ReferenceConstPointer;

  /**  Pointer type for the target. */
  typedef typename TargetType::ConstPointer TargetConstPointer;

  /**  Pointer type for the transformation. */
  typedef  typename TransformationType::Pointer TransformationPointer;

  /**  Pointer type for the metric. */
  typedef  typename MetricType::Pointer  MetricPointer;

  /**  Pointer type for the mapper. */
  typedef typename MapperType::Pointer  MapperPointer;

  /** Method that initiates the registration. */
  void StartRegistration( void );

  /** Set the transform parameters. */
  void SetParameters( const ParametersType&  params )
    { m_Parameters = params; }

  /** Get the transform parameters. */
  const ParametersType& GetParameters( void ) const
    { return m_Parameters; }

  /** Set the learning rate. This is used in the optimization scheme.
   * Typically, the learning rate needs to be annealed (decreased over
   * time). Default is 1.0. */
  itkSetClampMacro( LearningRate, double, 0.0,
                    NumericTraits<double>::max() );

  /** Get the learning rate. */
  itkGetMacro( LearningRate, double );

  /** Set the number of iterations. This determines the number of
   * iterations performed in the steepest descent optimization.
   * Default is 1000. */
  itkSetMacro( NumberOfIterations, unsigned int );

  /** Get the number of iterations. */
  itkGetMacro( NumberOfIterations, unsigned int );

protected:
  ImageToImageAffineMutualInformationGradientDescentRegistration();
  virtual ~ImageToImageAffineMutualInformationGradientDescentRegistration();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  ImageToImageAffineMutualInformationGradientDescentRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  ParametersType             m_Parameters;

  /** Optimization related variables. */
  double                     m_LearningRate;
  unsigned int               m_NumberOfIterations;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToImageAffineMutualInformationGradientDescentRegistration.txx"
#endif

#endif


