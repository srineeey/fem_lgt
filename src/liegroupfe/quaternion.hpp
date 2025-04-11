#ifndef FILE_QUATERNION
#define FILE_QUATERNION

#include <bla.hpp>

namespace ngsolve
{
  using namespace ngbla;
  
  template <typename T = double>
  class Quaternion{
  protected:
      T scal;
      Vec<3, T> vec;

  public:
      Quaternion(T _scal, Vec<3, T> _vec)
              : scal(_scal), vec(_vec) { ; }

      Quaternion(Vec<4, T> _vec)
              : scal(_vec[0]), vec{_vec[1], _vec[2], _vec[3]} { ; }

      auto getScal() const { return scal; }

      auto getVec() const { return vec; }

      //addition of quaternion (quaternion+quaternion=quaternion)
      Quaternion operator+(Quaternion<T> q2) { return Quaternion<T>{scal + q2.getScal(), vec + q2.getVec()}; }

      //TODO:TEST
      //quaternion product (quaternion*quaternion=quaternion)
      Quaternion<T> operator*(Quaternion<T> q2) const {
          return Quaternion<T>{
                  //getScal() * q2.getScal() - getVec()[0] * q2.getVec()[0] - getVec()[1] * q2.getVec()[1] - getVec()[2] * q2.getVec()[2],
                  getScal() * q2.getScal() - ( InnerProduct(getVec(), q2.getVec()) ),
                  /*
                  Vec<3, T>{
                          getScal() * q2.getVec()[0] + getVec()[0] * q2.getScal() + getVec()[1] * q2.getVec()[2] - getVec()[2] * q2.getVec()[1],
                          getScal() * q2.getVec()[1] + getVec()[1] * q2.getScal() - getVec()[0] * q2.getVec()[2] + getVec()[2] * q2.getVec()[0],
                          getScal() * q2.getVec()[2] + getVec()[2] * q2.getScal() + getVec()[0] * q2.getVec()[1] - getVec()[1] * q2.getVec()[0],
                  }
                   */
                  //TODO: check sign of cross product!
                  getScal()*q2.getVec() + q2.getScal()*getVec() + Cross(getVec(), q2.getVec())
          };
      };

      Quaternion<T> conjugate() const {
          return Quaternion<T>{getScal(), -1.*getVec()};
      };

  };


  template <typename T = double>
  class UnitQuaternion : public Quaternion<T>  {
  public:

      UnitQuaternion(T _scal, Vec<3, T> _vec) : Quaternion<T>{_scal, _vec}
      {
          //TODO: normalize values?
          //(*this) = normalize(Quaternion<T>{_scal, _vec});
      }

      UnitQuaternion(T _scal, Vec<3, T> _vec, bool tonormalize) : Quaternion<T>{_scal, _vec}
      {
        if (tonormalize){
          //TODO: normalize values?
          (*this) = normalize(Quaternion<T>{_scal, _vec});
        }
      }

      UnitQuaternion(Vec<4, T> _vec) : Quaternion<T>{_vec}
      {
          //TODO: normalize values?
          //(*this) = normalize(Quaternion<T>{this->getScal(), this->getVec()});
      }

      UnitQuaternion(Vec<4, T> _vec, bool tonormalize) : Quaternion<T>{_vec}
      {
        if (tonormalize){
          //TODO: normalize values?
          (*this) = normalize(Quaternion<T>{this->getScal(), this->getVec()});
        }
      }

      UnitQuaternion(Quaternion<T> q) : Quaternion<T>(q.getScal(), q.getVec())
      {
          //TODO: normalize values?
          //(*this) = normalize(q);
      }

      UnitQuaternion(Quaternion<T> q, bool tonormalize) : Quaternion<T>(q.getScal(), q.getVec())
      {
        if (tonormalize){
          //TODO: normalize values?
          (*this) = normalize(q);
        }
      }

      //TODO: is Quaternion() overload necessary?
      operator Quaternion<T>() const {return Quaternion<T>{UnitQuaternion<T>::scal, UnitQuaternion<T>::vec};};

      //TODO:TEST
      //TODO: is this automatically chosen for UnitQuaternions?
      //quaternion product (quaternion*quaternion=quaternion)
      UnitQuaternion<T> operator*(UnitQuaternion<T> q2) const {
          return UnitQuaternion<T>{Quaternion<T>{(*this)}*Quaternion<T>{q2}};
      };

      UnitQuaternion<T> conjugate() const {
          return UnitQuaternion<T>{Quaternion<T>::conjugate()};
      };

  };

template  <typename T>
//scalar multiplication of quaternion (real*quaternion=quaternion)
auto operator*(T c, Quaternion<T> q) { return Quaternion<T>{c * q.getScal(), c * q.getVec()}; };

template  <typename T>
//quaternion dot product = R4 dot product (quaternion*quaternion=real)
T InnerProduct(Quaternion<T> q1, Quaternion<T> q2) {
    return ( q1.getScal() * q2.getScal() + q1.getVec()[0] * q2.getVec()[0] + q1.getVec()[1] * q2.getVec()[1] + q1.getVec()[2] * q2.getVec()[2] );
    //return ( q1.getScal() * q2.getScal() + InnerProduct(q1.getVec(), q2.getVec()) );
};

template  <typename T>
//normalize the quaternion to get unit quaternion
UnitQuaternion<T> normalize(Quaternion<T> q) {
    auto norm = sqrt(InnerProduct(q,q));
    //TODO: what if q = 0?
    auto invnorm = 1. / norm;
    return UnitQuaternion<T>{invnorm*q.getScal(), invnorm*q.getVec(), false};
};
    


//calculate unit quaternion representing a 3D rotation
template  <typename T>
UnitQuaternion<T> GetRotQuaternion(T angle, Vec<3, T> axis)
{
    return UnitQuaternion<T>{cos(0.5*angle), (sin(0.5*angle)/sqrt(InnerProduct(axis,axis)))*axis};
};

//slerp interpolation quaternion
//(1/sin(theta)) * [ sin((1-t)theta) * q_0 + sin(t theta) * q_0 ]
template  <typename T>
UnitQuaternion<T> SLERPInterpolate (double t, UnitQuaternion<T> a, UnitQuaternion<T> b)
{
    //get the relative 4D angle theta between a and b
    //TODO: are a and b normalized? should be by design if SO3
    T theta = InnerProduct(a, b);

    UnitQuaternion<T> _inter_q = sin((1.- t)*theta) * a + sin(t*theta) * b;
    _inter_q = (1./sin(theta)) * _inter_q;

    return _inter_q;   // tbd
};



template <typename T>
inline ostream & operator<< (ostream & ost, const Quaternion<T> & q)
{
  ost << q.getScal()
    << ", " << q.getVec()(0)
    << ", " << q.getVec()(1)
    << ", " << q.getVec()(2);

return ost;
};

}




#endif
