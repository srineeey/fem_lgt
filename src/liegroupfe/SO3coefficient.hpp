/*********************************************************************/
/* File:   SO3coefficient.cpp                                        */
/* Author: Srinath Bulusu                                            */
/* Date:   21. Oct. 2022                                             *///
/*********************************************************************/

/*
   SO3 Finite Element Coefficient Function
*/

#include <core/register_archive.hpp>
#include <fem.hpp>
// #include <coefficient_impl.hpp>

// #include <../ngstd/evalfunc.hpp>
#include <algorithm>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp> // for shallow archive
#endif // NGS_PYTHON

#include "quaternion.hpp"
#include "liegroups.hpp"
#include "liegroupfe.hpp"

#ifndef LIEGROUPFEM_SO3COEFFICIENT_H
#define LIEGROUPFEM_SO3COEFFICIENT_H



namespace ngfem
{

    shared_ptr<CoefficientFunction> VecCrossMatCF (shared_ptr<CoefficientFunction> coef);
    shared_ptr<CoefficientFunction> SO3CF (shared_ptr<CoefficientFunction> coef);


    class VecCrossMatCoefficientFunction : public T_CoefficientFunction<VecCrossMatCoefficientFunction>
    {
        shared_ptr<CoefficientFunction> veccf;
        using BASE = T_CoefficientFunction<VecCrossMatCoefficientFunction>;
        public:
        VecCrossMatCoefficientFunction() = default;
        VecCrossMatCoefficientFunction (shared_ptr<CoefficientFunction> veccf_in)
            : T_CoefficientFunction<VecCrossMatCoefficientFunction>(3*3, veccf_in->IsComplex()), veccf(veccf_in)
        {
            this->SetDimensions (ngstd::INT<2> (3,3));
        }

        void DoArchive(Archive &ar) override
        {
            CoefficientFunction::DoArchive(ar);
            ar.Shallow(veccf);
        }

        virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override
        {
            veccf->TraverseTree(func);
            func(*this);
        }

        using T_CoefficientFunction<VecCrossMatCoefficientFunction>::Evaluate;
        template <typename MIR, typename T, ORDERING ORD>
        void T_Evaluate (const MIR & mir,
                        BareSliceMatrix<T,ORD> result) const
        {
            // int hd = Dimensions()[0]; //should be 3
            int hd = 3;
            STACK_ARRAY(T, hmem, mir.Size()*hd);
            FlatMatrix<T,ORD> tempvec(hd, mir.Size(), &hmem[0]);
            veccf->Evaluate (mir, tempvec);

            for (size_t i = 0; i < mir.Size(); i++)
            {  
                result(0*hd+0, i) = 0.;
                result(1*hd+1, i) = 0.;
                result(2*hd+2, i) = 0.;

                result(0*hd+1, i) = tempvec(2,i);
                result(1*hd+0, i) = -1.*tempvec(2,i);

                result(0*hd+2, i) = -1.*tempvec(1,i);
                result(2*hd+0, i) = tempvec(1,i);

                result(1*hd+2, i) = tempvec(0,i);
                result(2*hd+1, i) = -1.*tempvec(0,i);
            }
        }  

        template <typename MIR, typename T, ORDERING ORD>
        void T_Evaluate (const MIR & ir,
                        FlatArray<BareSliceMatrix<T,ORD>> input,                       
                        BareSliceMatrix<T,ORD> values) const
        {
            int hd = Dimensions()[0];
            
            auto in0 = input[0];

            for (size_t i = 0; i < ir.Size(); i++)
            {  
                values(0*hd+0, i) = 0.;
                values(1*hd+1, i) = 0.;
                values(2*hd+2, i) = 0.;

                values(0*hd+1, i) = in0(2,i);
                values(1*hd+0, i) = -1.*in0(2,i);

                values(0*hd+2, i) = -1.*in0(1,i);
                values(2*hd+0, i) = in0(1,i);

                values(1*hd+2, i) = in0(0,i);
                values(2*hd+1, i) = -1.*in0(0,i);
            }

        }

        shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                                shared_ptr<CoefficientFunction> dir) const override
        {
            if (this == var) return dir;
            return VecCrossMatCF(veccf->Diff(var, dir));
        }
    
    };
  

    class SO3CoefficientFunction : public T_CoefficientFunction<SO3CoefficientFunction>
    {
        //generate SO3 CF from R4 CF qcf
        shared_ptr<CoefficientFunction> r4cf;
        size_t matdim = 3;

    public:
        SO3CoefficientFunction() = default;
        SO3CoefficientFunction (shared_ptr<CoefficientFunction> r4cf_in)
                : T_CoefficientFunction<SO3CoefficientFunction>(3*3, r4cf_in->IsComplex()), r4cf(r4cf_in), matdim(3)
        {
            this->SetDimensions (ngstd::INT<2> (3,3));
        }


        using T_CoefficientFunction<SO3CoefficientFunction>::Evaluate;
        virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
        {
            throw Exception ("SO3CF:: scalar evaluate for matrix called");
        }

        void DoArchive(Archive &ar) override
        {
            CoefficientFunction::DoArchive(ar);
            ar.Shallow(r4cf);
        }

        virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override
        {
            r4cf->TraverseTree(func);
            func(*this);
        }


        template <typename MIR, typename T, ORDERING ORD>
        void T_Evaluate (const MIR & mir,
                        BareSliceMatrix<T,ORD> result) const
        {
            //r4cf->Evaluate (mir, result);

            STACK_ARRAY(T, hmem, mir.Size()*4);
            FlatMatrix<T,ORD> tempq(4, mir.Size(), &hmem[0]);
            r4cf->Evaluate (mir, tempq);

            for (size_t i = 0; i < mir.Size(); i++)
            {
                Vec<4,T> _q;
                for (size_t k = 0; k < 4; k++)
                    // _q(k) = result(k, i);
                    _q(k) = tempq(k, i);
                
                auto _SO3mat = (ngsolve::SO3<T>(ngsolve::UnitQuaternion<T>(_q, true))).GetMatrix();
                for (size_t j = 0; j < matdim; j++)
                    for (size_t k = 0; k < matdim; k++)
                        result(j*matdim+k, i) = _SO3mat(j,k);
            }

        }  

        template <typename MIR, typename T, ORDERING ORD>
        void T_Evaluate (const MIR & ir,
                        FlatArray<BareSliceMatrix<T,ORD>> input,                       
                        BareSliceMatrix<T,ORD> values) const
        {
            auto in0 = input[0];
            for (size_t i = 0; i < ir.Size(); i++)
            {
                Vec<4,T> _q;
                for (size_t k = 0; k < 4; k++)
                    _q(k) = in0(k, i);
                
                auto _SO3mat = ngsolve::SO3<T>(ngsolve::UnitQuaternion<T>(_q, true)).GetMatrix();
                for (size_t j = 0; j < matdim; j++)
                    for (size_t k = 0; k < matdim; k++)
                        values(j*matdim+k, i) = _SO3mat(j,k);
            }
        }

        shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                              shared_ptr<CoefficientFunction> dir) const override
        {
            //return (dq q Ad_q + Ad_q (dq q)^*) as a matrix operation
            //simple formula: s = dq q -> diff = 2 (s_0*Id + svec cross)R
            if (this == var)
            return  dir;

            auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());

            auto qs = MakeComponentCoefficientFunction(r4cf, 0);
            auto qv = MakeSubTensorCoefficientFunction(r4cf, 1, Array<int>{3}, Array<int>{1}); //size or list argument?
            auto ps = MakeComponentCoefficientFunction(dir, 0);
            auto pv = MakeSubTensorCoefficientFunction(dir, 1, Array<int>{3}, Array<int>{1});
            auto ss = qs*ps + InnerProduct(qv, pv); 
            auto sv = qs*pv - ps*qv - CrossProduct(pv, qv);
            auto fac = 2.* (ss * IdentityCF(3) + VecCrossMatCF(sv));
            return fac * thisptr;
        }

    
    };



}

shared_ptr<ngfem::CoefficientFunction> ngfem::VecCrossMatCF (shared_ptr<CoefficientFunction> coef)
    {
        if (coef->IsZeroCF())
            return ngfem::ZeroCF(Array<int>());

        return make_shared<ngfem::VecCrossMatCoefficientFunction> (coef);
    };

shared_ptr<ngfem::CoefficientFunction> ngfem::SO3CF (shared_ptr<CoefficientFunction> coef)
    {
        return make_shared<ngfem::SO3CoefficientFunction> (coef);
    };

#endif //LIEGROUPFEM_SO3COEFFICIENT_H
