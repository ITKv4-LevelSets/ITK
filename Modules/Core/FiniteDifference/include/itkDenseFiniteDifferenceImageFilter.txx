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
#ifndef __itkDenseFiniteDifferenceImageFilter_txx
#define __itkDenseFiniteDifferenceImageFilter_txx
#include "itkDenseFiniteDifferenceImageFilter.h"

#include <list>
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
void
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::CopyInputToOutput()
{
  typename TInputImage::ConstPointer input  = this->GetInput();
  typename TOutputImage::Pointer output = this->GetOutput();

  if ( !input || !output )
    {
    itkExceptionMacro(<< "Either input and/or output is NULL.");
    }

  // Check if we are doing in-place filtering
  if ( this->GetInPlace() && this->CanRunInPlace() )
    {
    typename TInputImage::Pointer tempPtr =
      dynamic_cast< TInputImage * >( output.GetPointer() );
    if ( tempPtr && tempPtr->GetPixelContainer() == input->GetPixelContainer() )
      {
      // the input and output container are the same - no need to copy
      return;
      }
    }

  ImageRegionConstIterator< TInputImage > in( input, output->GetRequestedRegion() );
  ImageRegionIterator< TOutputImage >     out( output, output->GetRequestedRegion() );

  while ( !out.IsAtEnd() )
    {
    out.Value() =  static_cast< PixelType >( in.Get() );  // Supports input
                                                          // image adaptors only
    ++in;
    ++out;
    }
}

template< class TInputImage, class TOutputImage >
void
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::AllocateUpdateBuffer()
{
  // The update buffer looks just like the output.
  typename TOutputImage::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetOrigin( output->GetOrigin() );
  m_UpdateBuffer->SetSpacing( output->GetSpacing() );
  m_UpdateBuffer->SetDirection( output->GetDirection() );
  m_UpdateBuffer->SetLargestPossibleRegion( output->GetLargestPossibleRegion() );
  m_UpdateBuffer->SetRequestedRegion( output->GetRequestedRegion() );
  m_UpdateBuffer->SetBufferedRegion( output->GetBufferedRegion() );
  m_UpdateBuffer->Allocate();
}

template< class TInputImage, class TOutputImage >
void
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::ApplyUpdate(TimeStepType dt)
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;

  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ApplyUpdateThreaderCallback,
                                            &str);
  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitely call Modified on GetOutput here
  // since ThreadedApplyUpdate changes this buffer
  // through iterators which do not increment the
  // output timestamp
  this->GetOutput()->Modified();
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::ApplyUpdateThreaderCallback(void *arg)
{
  DenseFDThreadStruct *str;
  int                  total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

  str = (DenseFDThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if ( threadId < total )
    {
    str->Filter->ThreadedApplyUpdate(str->TimeStep, splitRegion, threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template< class TInputImage, class TOutputImage >
typename
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >::TimeStepType
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::CalculateChange()
{
  int          threadCount;
  TimeStepType dt;

  // Set up for multithreaded processing.
  DenseFDThreadStruct str;

  str.Filter = this;
  str.TimeStep = NumericTraits< TimeStepType >::Zero;  // Not used during the
  // calculate change step.
  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->CalculateChangeThreaderCallback,
                                            &str);

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  threadCount = this->GetMultiThreader()->GetNumberOfThreads();
  str.TimeStepList = new TimeStepType[threadCount];
  str.ValidTimeStepList = new bool[threadCount];
  for ( int i = 0; i < threadCount; ++i )
    {
    str.ValidTimeStepList[i] = false;
    }

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  dt = this->ResolveTimeStep(str.TimeStepList, str.ValidTimeStepList, threadCount);
  delete[] str.TimeStepList;
  delete[] str.ValidTimeStepList;

  // Explicitely call Modified on m_UpdateBuffer here
  // since ThreadedCalculateChange changes this buffer
  // through iterators which do not increment the
  // update buffer timestamp
  this->m_UpdateBuffer->Modified();

  return dt;
}

template< class TInputImage, class TOutputImage >
ITK_THREAD_RETURN_TYPE
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::CalculateChangeThreaderCallback(void *arg)
{
  DenseFDThreadStruct *str;
  int                  total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

  str = (DenseFDThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;

  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if ( threadId < total )
    {
    str->TimeStepList[threadId] =
      str->Filter->ThreadedCalculateChange(splitRegion, threadId);
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_VALUE;
}

template< class TInputImage, class TOutputImage >
void
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType & regionToProcess,
                      int)
{
  ImageRegionIterator< UpdateBufferType > u(m_UpdateBuffer,    regionToProcess);
  ImageRegionIterator< OutputImageType >  o(this->GetOutput(), regionToProcess);

  u = u.Begin();
  o = o.Begin();

  while ( !u.IsAtEnd() )
    {
    o.Value() += static_cast< PixelType >( u.Value() * dt );  // no adaptor
                                                              // support here
    ++o;
    ++u;
    }
}

template< class TInputImage, class TOutputImage >
typename
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >::TimeStepType
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::ThreadedCalculateChange(const ThreadRegionType & regionToProcess, int)
{
  typedef typename OutputImageType::RegionType                    RegionType;
  typedef typename OutputImageType::SizeType                      SizeType;
  typedef typename OutputImageType::IndexType                     IndexType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType NeighborhoodIteratorType;

  typedef ImageRegionIterator< UpdateBufferType > UpdateIteratorType;

  typename OutputImageType::Pointer output = this->GetOutput();

  TimeStepType timeStep;
  void *       globalData;

  // Get the FiniteDifferenceFunction to use in calculations.
  const typename FiniteDifferenceFunctionType::Pointer df =
    this->GetDifferenceFunction();
  const SizeType radius = df->GetRadius();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because input has been copied to output.
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< OutputImageType >
  FaceCalculatorType;

  typedef typename FaceCalculatorType::FaceListType FaceListType;

  FaceCalculatorType faceCalculator;

  FaceListType faceList = faceCalculator(output, regionToProcess, radius);
  typename FaceListType::iterator fIt = faceList.begin();

  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  globalData = df->GetGlobalDataPointer();

  // Process the non-boundary region.
  NeighborhoodIteratorType nD(radius, output, *fIt);
  UpdateIteratorType       nU(m_UpdateBuffer,  *fIt);
  nD.GoToBegin();
  while ( !nD.IsAtEnd() )
    {
    nU.Value() = df->ComputeUpdate(nD, globalData);
    ++nD;
    ++nU;
    }

  // Process each of the boundary faces.

  NeighborhoodIteratorType bD;
  UpdateIteratorType       bU;
  for ( ++fIt; fIt != faceList.end(); ++fIt )
    {
    bD = NeighborhoodIteratorType(radius, output, *fIt);
    bU = UpdateIteratorType  (m_UpdateBuffer, *fIt);

    bD.GoToBegin();
    bU.GoToBegin();
    while ( !bD.IsAtEnd() )
      {
      bU.Value() = df->ComputeUpdate(bD, globalData);
      ++bD;
      ++bU;
      }
    }

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  timeStep = df->ComputeGlobalTimeStep(globalData);
  df->ReleaseGlobalDataPointer(globalData);

  return timeStep;
}

template< class TInputImage, class TOutputImage >
void
DenseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
