#ifndef FILE_LIEGROUPS
#define FILE_LIEGROUPS


#include "quaternion.hpp"

namespace ngsolve
{
  using namespace ngbla;
  
  template <typename T = double>
  class SO3
  {
    UnitQuaternion<T> q;

  public:
    SO3 (UnitQuaternion<T> _q) : q(_q) { }

    UnitQuaternion<T> getq() const {return q;}

    //get SO(3) matrix expression
    Mat<3,3,T> GetMatrix() const {
        //TODO: TEST
        T q0 = q.getScal();
        Vec<3,T> qvec{q.getVec()};

        Mat<3,3,T> _qmat;

        for(int k = 0; k <=2 ; k++)
        {
            for(int l = 0; l <=2 ; l++)
            {
                // _qmat(k, l) += qvec[k]*qvec[l];
                _qmat(k, l) = qvec[k]*qvec[l];
            }
        }

        //TODO: calculate cross product by recasting Mat<k,.> into Vec<.>?
        //cross product terms
        _qmat(0, 1) = _qmat(0, 1) - q0*qvec[2];
        _qmat(1, 0) = _qmat(1, 0) + q0*qvec[2];

        _qmat(2, 0) = _qmat(2, 0) - q0*qvec[1];
        _qmat(0, 2) = _qmat(0, 2) + q0*qvec[1];

        _qmat(1, 2) = _qmat(1, 2) - q0*qvec[0];
        _qmat(2, 1) = _qmat(2, 1) + q0*qvec[0];

        auto qvecqvec = qvec[0]*qvec[0] + qvec[1]* qvec[1] + qvec[2] *qvec[2];

        for(int k = 0; k <=2 ; k++)
        {
            _qmat(k, k) = _qmat(k, k) - qvecqvec;
            //_qmat(k, k) = _qmat(k, k) - InnerProduct(qvec, qvec);
        }

        _qmat = (2./(q.getScal()*q.getScal() + qvecqvec)) * _qmat;
        //_qmat = (2./InnerProduct(q, q)) * _qmat;

        //add identity
        _qmat(0, 0) = _qmat(0, 0) + 1.;
        _qmat(1, 1) = _qmat(1, 1) + 1.;
        _qmat(2, 2) = _qmat(2, 2) + 1.;
        //_qmat = _qmat + Id<3>;
        return _qmat;
    } // tbd

    //TODO: return type?
    //Mat<3,3> ApplyMatrix(Vec<3,T> v) const {
    Vec<3,T> ApplyMatrix(Vec<3,T> v) const {
        //R(v) = qvq^*

        //embed v into quaternion space
        Quaternion<T> q_v{0., v};
        //Rotate in quaternion space
        q_v = q * q_v * (q.conjugate());
        //q_v = q * q_v ;
        //q_v = q_v * (q.conjugate());
        //retract/extract 3D vector
        return q_v.getVec();
      } // tbd

    SO3<T> Inverse() const {
        return SO3<T>(getq().conjugate());
    }

    SO3<T> operator* (SO3<T> rot2) const {return SO3<T>(getq()*rot2.getq());}

  };

  template <typename T = double>
  class so3
    {
        Vec<3,T> sigma;

    public:
        so3(Vec<3,T> _sigma) : sigma(_sigma) { }

        so3(T _x, T _y, T _z) : sigma(Vec<3,T>{_x, _y, _z}) { }

        Vec<3,T> getvec() const {return sigma;}

        //get so(3) matrix expression
        Mat<3,3,T> GetMatrix() const {
            //TODO: implement matrix formula
            //TODO: initialized with 0?
            //Mat<3,3,T> sigma_mat{0.0};
            Mat<3,3,T> sigma_mat;
            //mat_bc = -sigma_a eps_abc
            for(size_t i=0; i<=2 ;i++)
            {
                Vec<3,T> _ei_vec(0., 0., 0.);
                _ei_vec(i) = -1.;
                Vec<3,T> _matvec = Cross(getvec(),_ei_vec);
                for(size_t j=0; j<=2 ;j++)
                {
                    sigma_mat(i,j) = _matvec(j);
                }
            }
            return sigma_mat;
        } // tbd

        Mat<3,3,T> ApplyMatrix(Vec<3,T> v) const {
            return Cross(getvec(),v);
        } // tbd

      so3 operator+ (so3 sigma2){return so3{this->getvec()+sigma2.getvec()};}

      so3 LieBracket(so3 sigma2)
      {
            //TODO: implement Lie bracket
            //structure constants correspond to 3D vector cross product
            return so3{ Cross(getvec(),sigma2.getvec())};
      }
    };

  template <typename T = double>
  so3<T> operator* (T s, so3<T> sigma){return so3<T>{s*sigma.getvec()};}
  

    //slerp interpolation quaternion
    //(1/sin(theta)) * [ sin((1-t)theta) * q_0 + sin(t theta) * q_0 ]
  template <typename T = double>
  inline SO3<T> Interpolate (T t, SO3<T> a, SO3<T> b)
  {
    return SO3<T>{SLERPInterpolate(t, a.getq(), b.getq())};   // tbd
  }

    //TODO: return Lie algebra element, not vector
    template <typename T = double>
    so3<T> Log(SO3<T> g){

        //TODO: what if q = 1?
        if( std::fabs((1-g.getq().getScal())) < std::numeric_limits<T>::epsilon())
        {
            return so3{Vec<3,T>{0.,0.,0.}};
        }
        T _phir = 2*acos(g.getq().getScal());
        Vec<3,T> lievec{g.getq().getVec()};
        return so3<T>{_phir*( 1./(sqrt(1-g.getq().getScal()*g.getq().getScal())) )*lievec};
    }

    //TODO: take Lie algebra element, not vector
    //TODO: Log and Exp do not match up exactly
    template <typename T = double>
    SO3<T> Exp(so3<T> sigma){

        //length of a Lie Algebra element = real rotation angle
        Vec<3,T> vec = sigma.getvec();
        T phir = sqrt(InnerProduct(vec, vec));

        /*
        if( std::fabs(phir) < std::numeric_limits<double>::epsilon())
        {
            return SO3{ UnitQuaternion<double>{1., Vec<3,double>{0.,0.,0.}}};
        }
         */
        //convert to quaternion (vec/phir normalizes vec!)
        //TODO: what if sigma=0?
        UnitQuaternion<T> uq{cos(0.5*phir), (sin(0.5*phir)/phir)*vec};

        return SO3<T>{uq};
    }


    template <typename T = double>
    inline ostream & operator<< (ostream & ost, const SO3<T> & g)
    {
        auto gmat = g.GetMatrix();
        ost << gmat(0,0) << " " << gmat(0,1) << " " << gmat(0,2) << "\n"
            << gmat(1,0) << " " << gmat(1,1) << " " << gmat(1,2) << "\n"
            << gmat(2,0) << " " << gmat(2,1) << " " << gmat(2,2) << "\n";
        return ost;
    };

    template <typename T = double>
    inline ostream & operator<< (ostream & ost, const so3<T> & sigma)
    {
        auto sigmamat = sigma.GetMatrix();
        ost << sigmamat(0,0) << " " << sigmamat(0,1) << " " << sigmamat(0,2) << "\n"
            << sigmamat(1,0) << " " << sigmamat(1,1) << " " << sigmamat(1,2) << "\n"
            << sigmamat(2,0) << " " << sigmamat(2,1) << " " << sigmamat(2,2) << "\n";
        return ost;
    };

  
}

#endif
