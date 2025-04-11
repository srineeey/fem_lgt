#include <algorithm>

#include <comp.hpp>
#include <python_comp.hpp>

#include "comp.hpp"


using namespace ngcomp;



extern "C" void lgt_module(py::object & res) {
  cout << "called lgt_module" << endl;
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    

  // BEGIN USER DEFINED CODE
 
 /// TODO: calculate mass matrix entry for weights
  auto CalcMassMat = [](std::shared_ptr<MeshAccess> ma)
  {
    /// TODO: get H1space of order 1
    //shared_ptr<H1HighOrderFESpace> fesh1(1);
    Flags flags;
    flags.SetFlag ("order", int(1) );
    
    //shared_ptr<H1HighOrderFESpace> fesh1_ptr = make_shared<H1HighOrderFESpace> (fesh1);
    auto fesh1 = make_shared<H1HighOrderFESpace>(ma, flags);
    
    auto u = fesh1->GetTrialFunction();
    auto up = fesh1->GetTestFunction();

    /// TODO: define bilinear form for the mass matrix
    Flags flags_massbf;
    auto massbf = make_shared<T_BilinearFormSymmetric<double>> (fesh1, "massh1", flags_massbf);

    /// TODO: add mass term (BilinearFormIntegrator)
    massbf->AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(u * up, VOL, VOL));


    /// TODO: (Re)Assemble()
    /// How large should LocalHeap be?
    
    //size_t heap_size = sizeof(double)*ma->GetNV()*ma->GetNV();
    size_t heap_size = sizeof(double)*10;
    LocalHeap lh(heap_size, "massbf", true);
    
    /// TODO: python kernel dies...
    massbf->Assemble(lh);


    /// TODO: GetMatrix() - compiler error: can not convert BaseMatrix to SparseMatrix
    
    SparseMatrix<double> massmat = massbf->GetMatrix();
    
    //auto& tbonemassmat = massbf->GetMatrix();
    
    /// Vnum = row number
    for(size_t Vnum = 0; Vnum < ma->GetNV(); Vnum++)
    {
        //std::cout <<  "Vnum " << Vnum << std::endl << std::flush;
        /// TODO: find way to index matrix object

        /// j = number of non-zero entry in row
        // for(size_t j = massmat.GetRowIndices[Vnum]; j < massmat.GetRowIndices[Vnum+1]; j++)
        // {
        //   if (massmat.col[j] == Vnum)
        //   {
        //     massmatvalue = tbonemassmat.GetRowValues(j);
        //     break;
        //   }
        // }
        

        // double massmatvalue = massmat(Vnum, Vnum);
        //std::cout <<  "massmatvalue " << massmatvalue << std::endl << std::flush;
    }
  };


m.def("CalcMassMat", CalcMassMat);

// END USER DEFINED CODE

res = m;  

}