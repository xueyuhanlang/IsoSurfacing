/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *      Copyright (C) 2001-2011 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openflipper.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenFlipper.                                        *
 *                                                                           *
 *  OpenFlipper is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenFlipper is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenFlipper. If not,                                  *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 10745 $                                                       *
 *   $Author: moebius $                                                      *
 *   $Date: 2011-01-26 10:23:50 +0100 (Mi, 26. Jan 2011) $                   *
 *                                                                           *
\*===========================================================================*/




//=============================================================================
//
//  CLASS NumLimitsT
//
//=============================================================================

#ifndef ACG_NUMLIMITS_HH
#define ACG_NUMLIMITS_HH

//#ifndef WIN32
#  undef min
#  undef max
//#endif


//== INCLUDES =================================================================

#include <limits.h>
#include <float.h>

//== NAMESPEACES ==============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================


/** \class NumLimitsT NumLimitsT.hh <ACG/Utils/NumLimitsT.hh

    This class provides the maximum and minimum values a certain
    scalar type (\c int, \c float, or \c double) can store. You can
    use it like this:
    \code
    #include <ACG/Utils/NumLimitsT.hh>

    int   float_min   = ACG::NumLimitsT<float>::min();
    float double_max  = ACG::NumLimitsT<double>::max();
    \endcode

    \note This functionality should be provided by
    std::numeric_limits.  This template does not exist on gcc <=
    2.95.3. The class template NumLimitsT is just a workaround.
**/
template <typename Scalar>
class NumLimitsT
{
public:
  /// Return the smallest \em absolte value a scalar type can store.
  static inline Scalar min() { return 0; }
  /// Return the maximum \em absolte value a scalar type can store.
  static inline Scalar max() { return 0; }
};


template<> inline int NumLimitsT<int>::min() { return INT_MIN; }
template<> inline int NumLimitsT<int>::max() { return INT_MAX; }

template<> inline float NumLimitsT<float>::min() { return FLT_MIN; }
template<> inline float NumLimitsT<float>::max() { return FLT_MAX; }

template<> inline double NumLimitsT<double>::min() { return DBL_MIN; }
template<> inline double NumLimitsT<double>::max() { return DBL_MAX; }


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NUMLIMITS_HH defined
//=============================================================================

