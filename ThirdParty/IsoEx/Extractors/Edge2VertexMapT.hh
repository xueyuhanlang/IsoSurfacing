/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

//=============================================================================
//
//  CLASS Edge2VertexMapT
//
//=============================================================================


#ifndef ISOEX_EDGE2VERTEXMAPT_HH
#define ISOEX_EDGE2VERTEXMAPT_HH


//== INCLUDES =================================================================

#include <map>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

/** \class Edge2VertexMapT Edge2VertexMapT.hh <IsoEx/Extractors/Edge2VertexMapT.hh>

    This class map edges (referenced by a point index and an axis) to
    vertex handles. Both PointIdx and VertexHandle are template
    arguments.  Internally this is implemented using a std::map, the
    key is a combination of point and axis index, the value is a
    VertexHandle.
    \note PointIdx must be comparable, i.e. have operator<().
*/	      

template <class PointIdx, class VertexHandle>
class Edge2VertexMapT
{
public:
   
  /// Constructor
  Edge2VertexMapT() {}


  /// clear the map
  void clear() { map_.clear(); }

  /// Store vertex in map
  void insert(PointIdx _p0, PointIdx _p1, VertexHandle _vhnd)
  {
    map_[EdgeKey(_p0, _p1)] = _vhnd;
  }

  /// Get vertex handle from map. Returns invalid handle if not found.
  VertexHandle find(PointIdx _p0, PointIdx _p1) const 
  {
    MyMapIterator it = map_.find(EdgeKey(_p0, _p1));
    if (it != map_.end())  return it->second;
    else return VertexHandle();
  }


private:

  class EdgeKey 
  {
  public:

    EdgeKey(PointIdx _p0, PointIdx _p1) {
      if (_p0 < _p1)  { p0_ = _p0;  p1_ = _p1; }
      else            { p0_ = _p1;  p1_ = _p0; }
    }

    bool operator<(const EdgeKey& _rhs) const 
    {
      if (p0_ != _rhs.p0_)
	return (p0_ < _rhs.p0_);
      else
	return (p1_ < _rhs.p1_);
    }

  private:
    PointIdx p0_, p1_;
  };


  typedef std::map<EdgeKey, VertexHandle>  MyMap;
  typedef typename MyMap::const_iterator   MyMapIterator;

  MyMap  map_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_EDGE2VERTEXMAPT_HH defined
//=============================================================================

