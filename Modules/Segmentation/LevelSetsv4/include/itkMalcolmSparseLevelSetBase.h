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

#ifndef __itkMalcolmSparseLevelSetBase_h
#define __itkMalcolmSparseLevelSetBase_h

#include "itkImage.h"
#include "itkLevelSetImageBase.h"

namespace itk
{
template< unsigned int VDimension >
class MalcolmSparseLevelSetBase :
    public LevelSetImageBase< Image< char, VDimension > >
{
public:
  typedef char                                    OutputType;
  typedef Image< OutputType, VDimension >         ImageType;
  typedef typename ImageType::Pointer             ImagePointer;

  typedef MalcolmSparseLevelSetBase               Self;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;
  typedef LevelSetImageBase< ImageType >          Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MalcolmSparseLevelSetBase, LevelSetImageBase);

  typedef typename Superclass::InputType      InputType;
  typedef typename Superclass::OutputRealType OutputRealType;
  typedef typename Superclass::GradientType   GradientType;
  typedef typename Superclass::HessianType    HessianType;

  typedef std::pair< InputType, OutputType >        NodePairType;
  typedef std::list< NodePairType >                 NodeListType;
  typedef typename NodeListType::iterator           NodeListIterator;
  typedef typename NodeListType::const_iterator     NodeListConstIterator;

  /*
  GradientType EvaluateGradient( const InputType& iP ) const
    {
    return GradientType();
    }

  HessianType EvaluateHessian( const InputType& iP ) const
    {
    return HessianType();
    }*/

  NodeListType* GetListNode()
    {
    return &m_List;
    }

protected:

  MalcolmSparseLevelSetBase() : Superclass()
  {}

  ~MalcolmSparseLevelSetBase()    {}

  NodeListType        m_List;

private:
  MalcolmSparseLevelSetBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkMalcolmSparseLevelSetBase_h
