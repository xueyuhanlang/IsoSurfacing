#ifndef ISOEX_EDGEKEY_HH
#define ISOEX_EDGEKEY_HH


//== INCLUDES =================================================================


//== NAMESPACES ===============================================================

namespace IsoEx
{

  template <class PointIdx>
  class EdgeKey
  {
  public:
    EdgeKey(){}
    ~EdgeKey(){}

    EdgeKey(PointIdx _p0, PointIdx _p1)
    {
      if (_p0 < _p1)  { p0_ = _p0;  p1_ = _p1; }
      else            { p0_ = _p1;  p1_ = _p0; }
    }


    enum Location {
      NONE,
      SAME,
      TOP,
      BOTTOM,
      LEFT,
      RIGHT,
      TOPLEFT,
      TOPRIGHT,
      BOTTOMLEFT,
      BOTTOMRIGHT
    };

    const PointIdx &p0() const {return p0_;}
    const PointIdx &p1() const {return p1_;}

    bool is_on_z_layer (unsigned int &_x_res,
                        unsigned int &_y_res,
                        unsigned int &_z_res,
                        unsigned int &_layer)
    {
      PointIdx idx = p0_;
      idx /= _x_res;
      idx /= _y_res;;
      unsigned int z0 = idx;

      idx = p1_;
      idx /= _x_res;
      idx /= _y_res;;
      unsigned int z1 = idx;

      if (z0 == z1)
      {
        _layer = z0;
        return true;
      }
      else return false;

    }

    bool is_x_edge ()
    {
      if ((p1_ - p0_) == 1) return true;
      else return false;
    }

    //     bool x_neighbor( const EdgeKey & _ne,
    //                      const unsigned int &_x_res,
    //                      const unsigned int &_y_res,
    //                      const unsigned int &_z_res,
    //                      int &_location )
    //     {
    //       //points of e
    //       OpenMesh::Vec2i p0;
    //       OpenMesh::Vec2i p1;
    //
    //       //points of ne
    //       OpenMesh::Vec2i np0;
    //       OpenMesh::Vec2i np1;
    //
    //       //directions of e and ne
    //       OpenMesh::Vec2i de;
    //       OpenMesh::Vec2i dne;
    //
    //       PointIdx idx = p0_;
    //       p0[0] = idx % _x_res; idx /= _x_res;
    //       p0[1] = idx % _y_res; idx /= _y_res;
    //
    //       idx = p1_;
    //       p1[0] = idx % _x_res; idx /= _x_res;
    //       p1[1] = idx % _y_res; idx /= _y_res;
    //
    //       idx = _ne.p0();
    //       np0[0] = idx % _x_res; idx /= _x_res;
    //       np0[1] = idx % _y_res; idx /= _y_res;
    //
    //       idx = _ne.p1();
    //       np1[0] = idx % _x_res; idx /= _x_res;
    //       np1[1] = idx % _y_res; idx /= _y_res;
    //
    //       de = p1 - p0;
    //       dne = np1 - np0;
    //
    //       _location = UNKNOWN;
    //       bool neighbor = true;
    //
    //       if (de == OpenMesh::Vec2i(1,0) )
    //       {
    //         if (dne == OpenMesh::Vec2i(1,0))
    //         {
    //           if (abs(p0_ - _ne.p0()) == 1 ) _location = SAME;
    //           else
    //           {
    //             if (p0[1] < np0[1])
    //             {
    //               neighbor = false;
    //               _location = TOP;
    //             }
    //             else if (p0[1] > np0[1])
    //             {
    //               neighbor = false;
    //               _location = BOTTOM;
    //             }
    //             else
    //             {
    //               neighbor = false;
    //               _location = SAME;
    //             }
    //           }
    //         }
    //         if (dne == OpenMesh::Vec2i(0,1))
    //         {
    //           if (np0 == p0)
    //           {
    //             _location = TOP;
    //             neighbor = false;
    //           }
    //           else if (np0 == p1)
    //           {
    //             _location = TOP;
    //             neighbor = false;
    //           }
    //           else if (np1 == p0)
    //           {
    //             _location = BOTTOM;
    //             neighbor = false;
    //           }
    //           else if (np1 == p1)
    //           {
    //             _location = BOTTOM;
    //             neighbor = false;
    //           }
    //           else
    //           {
    //             neighbor = false;
    //             if (np0[1] > p0[1]) _location = TOP;
    //             else _location = BOTTOM;
    //           }
    //         }
    //       }
    //       else  neighbor = false;
    //
    //       return neighbor;
    //     }

    void x_neighbor( const EdgeKey & _ne,
                     const unsigned int &_x_res,
                     const unsigned int &_y_res,
                     const unsigned int &_z_res,
                     int &_location )
    {
      //points of e
      OpenMesh::Vec2i p0;
      OpenMesh::Vec2i p1;

      //points of ne
      OpenMesh::Vec2i np0;
      OpenMesh::Vec2i np1;

      //directions of e and ne
      OpenMesh::Vec2i de;
      OpenMesh::Vec2i dne;

      PointIdx idx = p0_;
      p0[0] = idx % _x_res; idx /= _x_res;
      p0[1] = idx % _y_res; idx /= _y_res;

      idx = p1_;
      p1[0] = idx % _x_res; idx /= _x_res;
      p1[1] = idx % _y_res; idx /= _y_res;

      idx = _ne.p0();
      np0[0] = idx % _x_res; idx /= _x_res;
      np0[1] = idx % _y_res; idx /= _y_res;

      idx = _ne.p1();
      np1[0] = idx % _x_res; idx /= _x_res;
      np1[1] = idx % _y_res; idx /= _y_res;

      de = p1 - p0;
      dne = np1 - np0;

      _location = NONE;

      if (de == OpenMesh::Vec2i(1,0) )
      {
        if (dne == OpenMesh::Vec2i(1,0))
        {
          if ((int)p0_ - (int)_ne.p0() == 1 ) _location = LEFT;
          else if ((int)p0_ - (int)_ne.p0() == -1) _location = RIGHT;
          else if ((int)p0_ - (int)_ne.p0() == 0) _location = SAME;
          else
          {
            if (p0[1] - np0[1] == -1)
            {
              _location = TOP;
            }
            else if (p0[1] - np0[1] == 1)
            {
              _location = BOTTOM;
            }
          }
        }
        else if (dne == OpenMesh::Vec2i(0,1))
        {
          if (np0 == p0)
          {
            _location = TOPLEFT;
          }
          else if (np0 == p1)
          {
            _location = TOPRIGHT;
          }
          else if (np1 == p0)
          {
            _location = BOTTOMLEFT;
          }
          else if (np1 == p1)
          {
            _location = BOTTOMRIGHT;
          }
        }
      } 
    }


    bool operator<(const EdgeKey& _rhs) const
    {
      if (p0_ != _rhs.p0())
        return (p0_ < _rhs.p0());
      else
        return (p1_ < _rhs.p1());
    }

  private:
    PointIdx p0_, p1_;
  };


  //=============================================================================
} // namespace IsoEx


//=============================================================================
#endif // ISOEX_EDGEKEY_HH defined
//=============================================================================
