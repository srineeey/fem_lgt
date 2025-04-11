//using namespace std::complex<T>_literals;

/// the number of coefficients used to represent a group element:
/// 4 because quaternion, even though su(2) is 3 dimensional
/// Lie algebra basis chosen: t_k = 1/i sigma_k k=x,y,z

/// sin(x)/x function that avoids division by 0
template<typename T = double>
T sinc(T x, T eps)
{
    if (x*x < eps*eps)
    {
        //std::cout << "using taylor at 0 for x = "<< x << std::endl;
        //return 1. - (x*x/ (T(6.)) ) + (x*x*x*x/ (T(120.)) );
        /// TODO: is this precise enough?
        return 1.;
    }
    else
    {
        return sin(x)/x;
    }
}

template<typename T = double>
T sinc(T x)
{
    return sinc(x, T(0.0001) );
}


/// complex datatype that can be vectorized
/// noticed issues with std::complex<std::vector> 
template<typename T = double>
class MyComplex
{
    private:
        std::pair<T,T> z;

    public:
        MyComplex<T>(): z(T(0.), T(0.)){};

        //MyComplex<T>(MyComplex<T> _zc): z(_zc.GetPair()){};
        MyComplex<T>(const MyComplex<T>& _zc): z(_zc.GetPair()){};
        
        //MyComplex<T>(std::pair<T,T> _z): z(_z){};
        MyComplex<T>(const std::pair<T,T> & _z): z(_z){};
        MyComplex<T>(std::pair<T,T> && _z): z(_z){};

        //MyComplex(T _re, T _im): z(_re, _im){};
        MyComplex(const T & _re, const T & _im): z(_re, _im){};
        MyComplex(T && _re, T && _im): z(_re, _im){};

        MyComplex(const T & _re): z(_re, T(0.)){
            //std::get<0>(z) = _re;
            //std::get<1>(z) = T(0.);
        };

        MyComplex<T> MyComplexPolar(T _r, T _phi){
            return MyComplex(_r*cos(_phi), _r*sin(_phi));
            };

        std::pair<T,T> GetPair() const {return z;};

        //T Re() const {return std::get<0>(z);};
        //T Im() const {return std::get<1>(z);};
        //T mod() const {return sqrt(Re()*Re() + Im()*Im());};
        //T arg() const {return atan2(Im(), Re());};

        T Re() const {return std::get<0>(z);};
        T Im() const {return std::get<1>(z);};
        T mod() const {return sqrt(Re()*Re() + Im()*Im());};
        T arg() const {return atan2(Im(), Re());};


        MyComplex<T> operator*(const MyComplex<T>& zl) const {            
            return MyComplex(Re()*zl.Re() - Im()*zl.Im(), Re()*zl.Im() + Im()*zl.Re());
            };

        MyComplex<T> Conj() const {
            return MyComplex(Re(), -1.*Im());
            };
        
        MyComplex<T> operator+(const MyComplex<T>& zl) const {
            return MyComplex(Re() + zl.Re(), Im() + zl.Im());
            };

        MyComplex<T> operator-(const MyComplex<T>& zl) const {
            return MyComplex(Re() - zl.Re(), Im() - zl.Im());
            };

        MyComplex<T> operator+=(const MyComplex<T>& zl) const {
            return MyComplex((*this) + zl);
            };

        MyComplex<T> operator-=(const MyComplex<T>& zl) const {
            return MyComplex((*this) - zl);
            };

        std::string toString() const{
            std::stringstream s;
            s << *this;
            return s.str();
        }

};

template <typename T = double>
MyComplex<T> operator*(T s, MyComplex<T> &zl)
{
            return MyComplex(s*zl.Re(), s*zl.Re());
};


template <typename T = double>
std::ostream& operator<<(std::ostream& out, const MyComplex<T>& z)
{

    out << z.Re();
    out << "+";
    out << z.Im();
    out << "i";
    return out;
};


using namespace ngbla;



template <typename T = double>
class su2
{
    /// coefficients in the chosen su2 basis (quaternion basis! t_k = 1/i sigma_k k=x,y,z)
    Vec<3,T> A;

public:
    su2() : A{T(0.),T(0.),T(0.)} { }
    
    su2(Vec<3,T> _A) : A(_A) { }

    su2(T _x, T _y, T _z) : A(Vec<3,T>{_x, _y, _z}) { }

    Vec<3,T> GetVec() const {return A;}


    Matrix<MyComplex<T>> GetMatrix() const 
    {
        Matrix<MyComplex<T>> A_mat(2,2);

        auto _vec = this->GetVec();

        /// Matrix = sum_k A[k] (1/i) sigma_k
        /// watch out for real and imaginary parts!
        // A_mat(0,0) = MyComplex(_vec[2], 0.);
        // A_mat(1,1) = MyComplex(-1.*_vec[2], 0.);
        // A_mat(0,1) = MyComplex(_vec[1], -1.*_vec[0]);
        // A_mat(1,0) = MyComplex(-1.*_vec[1], -1.*_vec[0]);

        A_mat(0,0) = MyComplex<T>( T(0.) , T(-1.)*_vec[2]);
        A_mat(1,1) = MyComplex<T>( T(0.), _vec[2]);
        A_mat(0,1) = MyComplex<T>(T(-1.)*_vec[1], -1.*_vec[0]);
        A_mat(1,0) = MyComplex<T>(_vec[1], -1.*_vec[0]);

        return A_mat;
    }

    /// Conj() is equivalent to inverse for SU2!
    su2<T> Conj() const {

        return su2<T>( T(-1.)*this->GetVec() );
    }

    Vec<3,T> Trta() const 
    {
        /// Trace(sigma_k sigma_l) = 2 delta_{kl} -> this implies for t_k = 1/i sigma_k k=x,y,z
        return T(-2.)*this->GetVec();
    }


    su2 operator+ (su2 A2) const
    {
        return su2{this->GetVec() + A2.GetVec()};
    }

    
};

template <typename T = double>
su2<T> LieBracket(su2<T> A1, su2<T> A2)
    {
        //structure constants correspond to 3D vector cross product
        // because effectively the quaternion basis is chosen
        // [sigma_k, sigma_l] = 2 i \eps_{klm} sigma_m -> this implies for t_k = 1/i sigma_k k=x,y,z
        return su2{ T(2.)*Cross(A1.GetVec(),A2.GetVec())};
    }

template <typename T = double>
su2<T> operator* (T s, su2<T> A){
    return su2<T>{s*(A.GetVec())};
    }

template <typename T = double>
T gInnerProduct(su2<T> A, su2<T> B)
{
    T IP{T(0.)};
    /// CAUTION: the inner product daggers the second entry, but we are working with REAL coefficients (vecs)
    for(int c = 0; c < 3; c++)
    {
        /// see su2::Trta() to understand formula ...
        IP += A.GetVec()[c] * B.GetVec()[c];
    }
    return T(2.) * IP;
}
  
template <typename T = double>
class SU2
{
    Matrix<MyComplex<T>> U{2,2};

public:
    SU2 () : U(Exp(su2<T>(Vec<3,T>(0., 0., 0.))).GetMatrix()){ }
    //SU2 (const SU2 &Uo): U(Uo.GetMatrix()){ }
    SU2 (const Matrix<MyComplex<T>>& _U) : U(_U){ }
    SU2 (const SU2<T>& _U) : U(_U.GetMatrix()){ }
    SU2 (const Vec<3,T>& _A) : U(Exp(su2<T>(_A)).GetMatrix()) { }
    SU2 (const su2<T>& _A) : U(Exp(_A).GetMatrix()) { }


    /// take four quaternion (!) coefficients
    /// _vec needs to be normalized
    SU2 (const Vec<4,T>& _vec)
    {
        U(0,0) = MyComplex<T>(_vec[0], -1.*_vec[3]);
        U(1,1) = MyComplex<T>(_vec[0], _vec[3]);
        U(0,1) = MyComplex<T>(-1.*_vec[2], -1.*_vec[1]);
        U(1,0) = MyComplex<T>(_vec[2], -1.*_vec[1]);
    }

    /// take four quaternion (!) coefficients
    /// need to be normalized
    SU2 (const T& u0, const T& u1, const T& u2, const T& u3)
    {  
        U(0,0) = MyComplex<T>(u0, -1.*u3);
        U(1,1) = MyComplex<T>(u0, u3);
        U(0,1) = MyComplex<T>(-1.*u2, -1.*u1);
        U(1,0) = MyComplex<T>(u2, -1.*u1);
    }

    Matrix<MyComplex<T>> GetMatrix() const {
        return U;
    }


    /// Conj() is equivalent to inverse for SU2!
    SU2<T> Conj() const {
        Matrix<MyComplex<T>> Uconjmat{2,2};

        Uconjmat(1,0) = U(0,1).Conj();
        Uconjmat(0,1) = U(1,0).Conj();
        Uconjmat(0,0) = U(0,0).Conj();
        Uconjmat(1,1) = U(1,1).Conj();

        return SU2<T>(Uconjmat);
    }

    MyComplex<T> Tr() const {   
        return (GetMatrix())(0,0) + (GetMatrix())(1,1);
    }

    T ReTr() const {   
        return ( (GetMatrix())(0,0) + (GetMatrix())(1,1) ).Re();
    }

    /// Opt for manual calculation of Trace(ta*U)
    /// instead of working with coefficients (would require taking logarithm!) 
    Vec<3,T> Trta() const {

        Vec<3,T> TrtaU{T(0.),T(0.),T(0.)};

        for(int a = 1; a <= 3; a++)
        {
            Vec<4,T> vecta{T(0.),T(0.),T(0.),T(0.)};
            vecta[a] = T(1.); 
            SU2<T> taU_ = SU2<T>(vecta) * (*this);
            TrtaU[a-1] = taU_.ReTr();
        }

        // Vec<4,T> UCoeff = GetCoeff();
        // TrtaU[0] = T(-2.)*UCoeff[1];
        // TrtaU[1] = T(-2.)*UCoeff[2];
        // TrtaU[2] = T(-2.)*UCoeff[3];
        
        return TrtaU;
    }

    /// extract quaternion (!) coefficients
    /// NOT Lie algebra element A which fulfills exp(A) = U!
    Vec<4,T> GetCoeff() const {

        Vec<4,T> UCoeff{T(0.),T(0.),T(0.),T(0.)};
        /// factor 0.5 in fundamental representation (trace of 2D identity)
        UCoeff[0] = T(0.5)*(this->ReTr());

        for(int a = 1; a <= 3; a++)
        {
            Vec<4,T> vecta{0., 0., 0., 0.};
            vecta[a] = T(1.); 
            SU2<T> taU_ = SU2<T>(vecta) * (*this);
            /// factor -0.5 because of quaternion basis t_k = 1/i sigma_k k=x,y,z
            UCoeff[a] = T(-0.5)*taU_.ReTr();
        }

        return UCoeff;
    }

    /// manually implement matrix multiplication
    SU2<T> operator* (const SU2<T>& U2) const {

        Matrix<MyComplex<T>> U1mat = this->GetMatrix();
        Matrix<MyComplex<T>> U2mat = U2.GetMatrix();
        Matrix<MyComplex<T>> U1U2mat(2,2);
        U1U2mat(0,0) = T(0.);
        U1U2mat(0,1) = T(0.);
        U1U2mat(1,0) = T(0.);
        U1U2mat(1,1) = T(0.);

        for(int r = 0; r<=1; r++)
        {
            for(int c = 0; c<=1; c++)
            {
                for(int i = 0; i<=1; i++)
                {
                    // += not working?
                    //U1U2mat(r,c) += U1mat(r,i)*U2mat(i,c);
                    U1U2mat(r,c) = U1U2mat(r,c) + U1mat(r,i)*U2mat(i,c);
                    //auto UU_ = U1mat(r,i)*U2mat(i,c);
                    //std::cout << U1mat(r,i) << " * " << U2mat(i,c) << " = " << UU_ << std::endl << std::flush;
                }   
            }   
        }

        return SU2<T>( Matrix<MyComplex<T>>(U1U2mat) );
    }

    std::string toString() const{
        std::stringstream s;
        s << *this;
        return s.str();
    }

};

/// << operator to print matrix values 
template<typename T>
std::ostream& operator<< (std::ostream& out, const SU2<T>& U){
    for(int r = 0; r<=1; r++)
    {
        for(int c = 0; c<=1; c++)
        {
            out << (U.GetMatrix())(r,c) << " ";   
        }   
        out << std::endl;
    }

    return out;
}


/// the functions below are respective "inverses"

/// expects: 0 < q0 < 1
template <typename T = double>
su2<T> Log(SU2<T> U)
{

    //Vec<3,T> lievec = T(-0.5) * U.Trta();

    Vec<4,T> lievec = U.GetCoeff();
    T q0 = lievec[0];
    T Amod = arccos(q0);
    T fac = Amod/sqrt(1-q0*q0);

    return su2<T>(Vec<3,T>(fac*lievec[0],fac*lievec[1],fac*lievec[2]));
}


template <typename T = double>
SU2<T> Exp(su2<T> A)
{

    //length of a Lie Algebra element = real rotation angle
    Vec<3,T> vec = A.GetVec();

    T phi = sqrt(InnerProduct(vec, vec));
    //Vec<3,T> SU2vec =  (sin(phi)/phi)*vec;
    Vec<3,T> SU2vec =  sinc(phi)*vec;

    return SU2{Vec<4,T>{cos(phi), SU2vec[0], SU2vec[1], SU2vec[2]}};


}
