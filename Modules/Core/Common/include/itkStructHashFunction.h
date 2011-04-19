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
#ifndef __itkStructHashFunction_h
#define __itkStructHashFunction_h

namespace itk
{
/** \class StructHashFunction
 *
 *  \brief Generic hash function for an arbitrary struct (or class).
 *
 *  This regards the input key as a string of bytes, and applies a
 *  hash function similar to one that has been used in perl.  If the
 *  data in the input key has pointers to other classes, then the
 *  function will be hashing on the value of the pointer itself, not
 *  on the data it points to.  It is imagined that this function will
 *  be used more for structs (with fully exposed data) than for
 *  general classes.
 * \ingroup ITK-Common
 */
template< class TInput >
class ITK_EXPORT StructHashFunction
{
public:

  /** Standard class typedefs. */
  typedef StructHashFunction Self;

  /** Input type */
  typedef TInput InputType;

  unsigned int operator()(const InputType & key) const;
};

template< class TInput >
inline unsigned int
StructHashFunction< TInput >
::operator()(const InputType & key) const
{
  unsigned int len = static_cast< unsigned int >( sizeof( InputType ) );
  const char * p = reinterpret_cast< const char * >( &key );
  unsigned int hash = 0;
  while ( len-- )
    {
    hash = hash * 65 + *p++;
    }
  hash += ( hash >> 5 );
  return hash;
}
}

#endif  // ndef itkStructHashFunction_h
