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

#ifndef __itkListPixel_h
#define __itkListPixel_h

#include <list>
#include <ostream>
#include "itkNumericTraits.h"

namespace itk
{
template< typename T >
class ListPixel : public std::list< T >
  {
public:
  typedef ListPixel       Self;
  typedef std::list< T >  Superclass;

  typedef typename Superclass::iterator       Iterator;
  typedef typename Superclass::const_iterator ConstIterator;

  ListPixel() : Superclass() {}
  ListPixel( const T& ) {}

  ~ListPixel() {}

  friend std::ostream & operator<< ( std::ostream & os, const Self & c)
  {
    os << "[ ";
    for( ConstIterator it = c.begin();
          it != c.end();
         ++it )
      {
      os << *it << " ";
      }
    os << "]" <<std::endl;
    return os;
  }

  };

template< typename T >
class NumericTraits< ListPixel< T > >
{
private:

  typedef typename NumericTraits< T >::AbsType        ElementAbsType;
  typedef typename NumericTraits< T >::AccumulateType ElementAccumulateType;
  typedef typename NumericTraits< T >::FloatType      ElementFloatType;
  typedef typename NumericTraits< T >::PrintType      ElementPrintType;
  typedef typename NumericTraits< T >::RealType       ElementRealType;
public:

  /** Return the type of the native component type. */
  typedef T ValueType;

  typedef ListPixel< T > Self;

  /** Unsigned component type */
  typedef ListPixel< ElementAbsType > AbsType;

  /** Accumulation of addition and multiplication. */
  typedef ListPixel< ElementAccumulateType > AccumulateType;

  /** Typedef for operations that use floating point instead of real precision
    */
  typedef ListPixel< ElementFloatType > FloatType;

  /** Return the type that can be printed. */
  typedef ListPixel< ElementPrintType > PrintType;

  /** Type for real-valued scalar operations. */
  typedef ListPixel< ElementRealType > RealType;

  /** Type for real-valued scalar operations. */
  typedef                                    ElementRealType ScalarRealType;

  /** Measurement vector type */
  typedef Self MeasurementVectorType;

  /** Component wise defined elements
   *
   * \note minimum value for floating pointer types is defined as
   * minimum positive normalize value.
   */
  static const Self max(const Self &)
  {
    return Self( NumericTraits< T >::max() );
  }

  static const Self min(const Self &)
  {
    return Self( NumericTraits< T >::min() );
  }

  static const Self max()
  {
    return Self( NumericTraits< T >::max() );
  }

  static const Self min()
  {
    return Self( NumericTraits< T >::min() );
  }

  static const Self NonpositiveMin()
  {
    return Self ( NumericTraits< ValueType >::NonpositiveMin() );
  }

  static const Self ZeroValue()
  {
    return Self(NumericTraits< T >::Zero);
  }

  static const Self OneValue()
  {
    return Self(NumericTraits< T >::One);
  }

  /** \note: the functions are prefered over the member variables as
   * they are defined for all partial specialization
   */
  static const Self ITKCommon_EXPORT Zero;
  static const Self ITKCommon_EXPORT One;
};

template< typename T >
const typename NumericTraits< ListPixel< T > >::Self
NumericTraits< ListPixel< T > >::Zero = Self();

template< typename T >
const typename NumericTraits< ListPixel< T > >::Self
NumericTraits< ListPixel< T > >::One = Self();
}
#endif // __itkListPixel_h
