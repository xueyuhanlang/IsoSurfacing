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
 *   $Revision: 14455 $                                                       *
 *   $Author: bommes $                                                      *
 *   $Date: 2012-04-27 18:48:14 +0200 (Fr, 27 Apr 2012) $                   *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS VectorT
//
//=============================================================================


#ifndef ACG_VECTOR_HH
#define ACG_VECTOR_HH


//== INCLUDES =================================================================


#include <OpenMesh/Core/Geometry/VectorT.hh>


//== NAMESPACES ===============================================================


namespace ACG {

//=============================================================================


using OpenMesh::VectorT;

/** 1-byte signed vector */
typedef VectorT<signed char,1> Vec1c;
/** 1-byte unsigned vector */
typedef VectorT<unsigned char,1> Vec1uc;
/** 1-short signed vector */
typedef VectorT<signed short int,1> Vec1s;
/** 1-short unsigned vector */
typedef VectorT<unsigned short int,1> Vec1us;
/** 1-int signed vector */
typedef VectorT<signed int,1> Vec1i;
/** 1-int unsigned vector */
typedef VectorT<unsigned int,1> Vec1ui;
/** 1-float vector */
typedef VectorT<float,1> Vec1f;
/** 1-double vector */
typedef VectorT<double,1> Vec1d;

/** 2-byte signed vector */
typedef VectorT<signed char,2> Vec2c;
/** 2-byte unsigned vector */
typedef VectorT<unsigned char,2> Vec2uc;
/** 2-short signed vector */
typedef VectorT<signed short int,2> Vec2s;
/** 2-short unsigned vector */
typedef VectorT<unsigned short int,2> Vec2us;
/** 2-int signed vector */
typedef VectorT<signed int,2> Vec2i;
/** 2-int unsigned vector */
typedef VectorT<unsigned int,2> Vec2ui;
/** 2-float vector */
typedef VectorT<float,2> Vec2f;
/** 2-double vector */
typedef VectorT<double,2> Vec2d;

/** 3-byte signed vector */
typedef VectorT<signed char,3> Vec3c;
/** 3-byte unsigned vector */
typedef VectorT<unsigned char,3> Vec3uc;
/** 3-short signed vector */
typedef VectorT<signed short int,3> Vec3s;
/** 3-short unsigned vector */
typedef VectorT<unsigned short int,3> Vec3us;
/** 3-int signed vector */
typedef VectorT<signed int,3> Vec3i;
/** 3-int unsigned vector */
typedef VectorT<unsigned int,3> Vec3ui;
/** 3-float vector */
typedef VectorT<float,3> Vec3f;
/** 3-double vector */
typedef VectorT<double,3> Vec3d;
/** 3-bool vector */
typedef VectorT<bool,3> Vec3b;

/** 4-byte signed vector */
typedef VectorT<signed char,4> Vec4c;
/** 4-byte unsigned vector */
typedef VectorT<unsigned char,4> Vec4uc;
/** 4-short signed vector */
typedef VectorT<signed short int,4> Vec4s;
/** 4-short unsigned vector */
typedef VectorT<unsigned short int,4> Vec4us;
/** 4-int signed vector */
typedef VectorT<signed int,4> Vec4i;
/** 4-int unsigned vector */
typedef VectorT<unsigned int,4> Vec4ui;
/** 4-float vector */
typedef VectorT<float,4> Vec4f;
/** 4-double vector */
typedef VectorT<double,4> Vec4d;

/** 6-byte signed vector */
typedef VectorT<signed char,6> Vec6c;
/** 6-byte unsigned vector */
typedef VectorT<unsigned char,6> Vec6uc;
/** 6-short signed vector */
typedef VectorT<signed short int,6> Vec6s;
/** 6-short unsigned vector */
typedef VectorT<unsigned short int,6> Vec6us;
/** 6-int signed vector */
typedef VectorT<signed int,6> Vec6i;
/** 6-int unsigned vector */
typedef VectorT<unsigned int,6> Vec6ui;
/** 6-float vector */
typedef VectorT<float,6> Vec6f;
/** 6-double vector */
typedef VectorT<double,6> Vec6d;


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_VECTOR_HH defined
//=============================================================================
