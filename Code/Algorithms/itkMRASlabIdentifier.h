/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkMRASlabIdentifier.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMRASlabIdentifier_h
#define __itkMRASlabIdentifier_h

#include "itkObject.h"
#include "itkImage.h"
#include <vector>


namespace itk
{

/** \class MRASlabIdentifier
 * \brief identifies slab in MR images comparing minimum intensity averages
 *
 * This class is templated over the type of image. 
 * In many cases, a 3D MR image is constructed by merging smaller 3D 
 * blocks (slabs) which were acquired with different settings such as magnetic
 * settings and patient positions. Therefore, stripe like patterns with slabs 
 * can be present in the resulting image. Such artifacts are called "slab 
 * boundary" artifacts or "venetian blind" artifacts.
 *
 * With the slab boundary artifacts in an image, even a same tissue class's
 * intensity values might vary significantly along the borders of slabs.
 * Such rough value changes are not appropriate for some image processing
 * methods. For example, MRIBiasFieldCorrectionFilter assumes a smooth bias 
 * field. However, with the slab boundary artifacts, the bias field estimation
 * scheme that MRIBiasFieldCorrectionFilter uses might not adopt well.
 * So, the MRIBiasFieldCorrectionFilter creates regions for slabs using the 
 * MRASlabIdentifier and then apply its bias correction scheme to each slab. 
 *
 * For this identifier, a slice means 2D image data which is extracted from
 * the input image along one of three axes (x, y, z). Users can specify
 * the slicing axis using the SetSlicingDirection(int dimension) member.
 * (0 - x, 1 - y, 2 - z).
 *
 * The identification scheme used here is very simple.  
 * 1) Users should specify how many pixels per slice the identifier
 *    will sample. 
 * 2) For each slice, the identifier searches the specified number of pixels 
 *    of which intensity values are greater than 0 and less than those 
 *    of the other pixels in the slice
 * 3) The identifier calculates the average for each slice and the overall
 *    average using the search results.
 * 4) For each slice, it subtracts the overall average from the slice average.
 *    If the sign of the subtraction result changes, then it assumes that a 
 *    slab ends and another slab begins.
 */
template <class TInputImage>            
class ITK_EXPORT MRASlabIdentifier : public Object 
{
public:
  /** Standard class typedefs. */
  typedef MRASlabIdentifier      Self;
  typedef Object  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Object);

  /** Type definition for the input image. */
  typedef TInputImage  ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer  ImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType ImagePixelType;
  typedef typename TInputImage::IndexType ImageIndexType;
  typedef typename TInputImage::SizeType ImageSizeType;
  typedef typename TInputImage::RegionType ImageRegionType ;
  typedef std::vector<ImageRegionType> SlabRegionVectorType ; 

  /** Set the input image. */
  void SetImage(ImagePointer image) ;

  /** Set the number of minimum intensity pixels per slice. */
  void SetNumberOfMinimumsPerSlice(int no) ;

  /** Get the number of minimum intensity pixels per slice. */
  int GetNumberOfMinimumsPerSlice(void) 
    {return m_NumberOfMinimumsPerSlice;}

  /** Set/Get the direction of slicing/
   * 0 - x axis, 1 - y axis, 2 - z axis */
  void SetSlicingDirection(int dimension)
    { m_SlicingDirection = dimension ; }
  int GetSlicingDirection()
    { return m_SlicingDirection ; }

  /** Compute the average values of miminum intensity pixels for each slice and
   * compare the average values with overall averages. */
  void GenerateSlabRegions(void);

  /** Get slab regions. */
  SlabRegionVectorType GetSlabRegionVector(void) ;

protected:
  MRASlabIdentifier() ;
  virtual ~MRASlabIdentifier() {} 
  
private:
  MRASlabIdentifier(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  /** Target image pointer that MRASlabIdentifier will use. */
  ImagePointer m_Image ;

  /** The number of pixels per slice which will be included 
   * for average calculation. In a sense, it's sampling size per slice. */
  int m_NumberOfMinimumsPerSlice ;
  
  int m_SlicingDirection ;

  SlabRegionVectorType m_Slabs ;
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRASlabIdentifier.txx"
#endif

#endif /* __itkMRASlabIdentifier_h */
