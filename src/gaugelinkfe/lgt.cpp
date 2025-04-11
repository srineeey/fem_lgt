#include <algorithm>

#include <comp.hpp>
#include <python_comp.hpp>

#include "comp.hpp"

#include "su2.hpp"
#include "MeshVars.hpp"

using namespace ngcomp;


/// PYBIND EXPORTS
/// START

void ExportSU2(py::module &m){

    //std::cout << "export SU2" << std::endl;
    py::class_<SU2<double>, std::shared_ptr<SU2<double>>>
        (m, "SU2", "docu.")
        .def("__init__",
            [](SU2<double> * instance)
            {
                new (instance) SU2<double>();
            }                
        )
        .def("__init__",
            [](SU2<double> * instance, double u0, double u1, double u2, double u3)
            {
                new (instance) SU2<double>{u0, u1, u2, u3};
            }, py::arg("u0")=1., py::arg("u1")=0., py::arg("u2")=0., py::arg("u3")=0., 
            "u0 to u3 values in a quaternion representation"    
        )
        .def("__init__",
            [](SU2<double> * instance, double A1, double A2, double A3, bool exp)
            {
                new (instance) SU2<double>(su2<double>(A1, A2, A3));
            }, py::arg("A1")=0., py::arg("A2")=0., py::arg("A3")=0., py::arg("expbool")=true, 
            "A values to exponentiate"    
        )
        .def("GetMatrix", [](SU2<double> * instance){
            return instance->GetMatrix();}
            )
        .def("Conj", [](SU2<double> * instance){
            return instance->Conj();}
            )
        .def("GetCoeff", [](SU2<double> * instance, int a){
            return (instance->GetCoeff())[a];}
            )
        .def("Trta", [](SU2<double> * instance){
            return instance->Trta();}
            )
        .def("__mul__", &SU2<double>::operator *)
        //.def("__mul__", &SU2<double>::operator *)
        .def("__str__", &SU2<double>::toString)
        ;
    
}

void ExportMyComplex(py::module &m){

    std::cout << "export MyComplex" << std::endl;
    py::class_<MyComplex<double>, std::shared_ptr<MyComplex<double>>>
        (m, "MyComplex", "docu.")
        .def("__init__",
            [](MyComplex<double> * instance)
            {
                new (instance) MyComplex<double>();
            }                
        )
        .def("__init__",
            [](MyComplex<double> * instance, double re, double im)
            {
                new (instance) MyComplex<double>(re, im);
            }, py::arg("re")=0., py::arg("im")=0., 
            "real part re and imaginary part im"    
        )
        .def("Conj", [](MyComplex<double> * instance){
            return instance->Conj();}
            )
        .def("__str__", &MyComplex<double>::toString)
        //.def_property("imag", &myComplex<double>::getImag, &myComplex<double>::setImag)
        ;
    
}


void ExportMeshVars(py::module &m){

    std::cout << "export MeshVars" << std::endl;
    py::class_<MeshVars, std::shared_ptr<MeshVars>>
        (m, "MeshVars", "docu.")
        .def("__init__",
            [](MeshVars * instance, std::shared_ptr<MeshAccess> _ma, double _dt)
            {
                new (instance) MeshVars(_ma, _dt);
            }, py::arg("_ma")=nullptr, py::arg("_dt")=1.,
            "shared pointer to mesh (MeshAccess) and time step"    
        )
        .def("GetMeshAccess",
            [](MeshVars * instance)
            {
                return instance->GetMeshAccess();
            }
        )
        .def("EvaluateFacetVorb",
            [](MeshVars * instance)
            {
                return instance->EvaluateFacetVorb();
            }
        )
        .def("FacetVorb",
            [](MeshVars * instance, size_t Inum)
            {
                return instance->FacetVorb(Inum);
            }
        )
        .def("EvaluateVertexVorb",
            [](MeshVars * instance)
            {
                return instance->EvaluateVertexVorb();
            }
        )
        .def("VertexVorb",
            [](MeshVars * instance, size_t Vnum)
            {
                return instance->VertexVorb(Vnum);
            }
        )
        .def("CalcMassDiagEntry",
            [](MeshVars * instance, size_t Vnum)
            {
                return instance->CalcMassDiagEntry(Vnum);
            }
        )
	.def("ElLinkOrientation",
            [](MeshVars * instance, const int& Inum, const int& Elnum)
            {
                return instance->ElLinkOrientation(Inum, Elnum);
            }
        )
        .def("CalcVertexOrientedEnums",
            [](MeshVars * instance, const size_t& Vnum)
            {
                return instance->CalcVertexOrientedEnums(Vnum);
            }
        )
        .def("GetOrientedEdgesOfVertex",
            [](MeshVars * instance, const size_t& Vnum)
            {
                return instance->GetOrientedEdgesOfVertex(Vnum);
            }
        )
        .def("GettBoneWeight",
            [](MeshVars * instance, const size_t& Vnum)
            {
                return instance->GettBoneWeight(Vnum);
            }
        )
        .def("GetBoneWeight",
            [](MeshVars * instance, const size_t& Inum)
            {
                return instance->GetBoneWeight(Inum);
            }
        )
        .def("RefreshBoneWeights",
            [](MeshVars * instance)
            {
                return instance->RefreshBoneWeights();
            }
        )
        .def("RefreshtBoneWeights",
            [](MeshVars * instance)
            {
                return instance->RefreshtBoneWeights();
            }
        )
        ;
    
}


extern "C" void lgt_module(py::object & res) {
  cout << "called lgt_module" << endl;
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    


ExportSU2(m);
ExportMyComplex(m);
ExportMeshVars(m);


/// PYBIND EXPORTS
/// STOP

// BEGIN USER DEFINED CODE


/*
Depending on where the Wilson loops starting point and direction, the oriented Edges need to be permuted and or inverted.
this function does exactly that.
dir == true: Link at Enum is multiplied FIRST
dir == false: Link at Enum is multiplied LAST
*/
auto ArrangeOrientedEdgesOfVertex = [](std::vector<std::tuple<int, bool>>& E_ors, const size_t& Vnum, const size_t& Enum, const bool& dir)
{

  // find index of Enum and return iterator
  auto eqEnum = [Enum](std::tuple<int, bool> _E_or){ return (std::get<0>(_E_or) == Enum); };
  auto Enumit = std::find_if(E_ors.begin(), E_ors.end(), eqEnum); // != vec.end()

  if ( dir == true )
  {
    /// permute Enum to the first position in E_ors
    std::rotate(E_ors.begin(), Enumit, E_ors.end());  
  }
  else if ( dir == false )
  {
    /// permute Enum to the last position of E_ors
    //std::rotate(E_ors.begin(), Enumit+1, E_ors.end());  

    /// permute Enum to the first position of E_ors
    std::rotate(E_ors.begin(), Enumit, E_ors.end());  
    
    /// then reverse vector
    std::reverse(E_ors.begin(), E_ors.end());

    /// and flip the direction bools!
    for (auto _E_or : E_ors)
    {
      std::get<1>(_E_or) = !(std::get<1>(_E_or));
    }
  }

  return E_ors;
};



/// GridFunction has 4 entries for each link
/// obtain link value at edge Inum
auto GetLink = [](std::shared_ptr<GridFunction> gfU, const int& Inum, bool forward)
{
  FlatVector<> gfUvec = gfU->GetVector().FV<double>();
  Array<DofId> dnums;
  gfU->GetFESpace()->GetEdgeDofNrs (Inum, dnums);
  Vec<4,double> vecU{gfUvec[dnums[0]], gfUvec[dnums[1]], gfUvec[dnums[2]], gfUvec[dnums[3]]};
  if (forward == false)
  {
      for (int a = 1; a <=3; a++)
      {
        vecU[a] *= -1.;
      }
  }

  //std::cout << vecU << std::endl << std::flush;
  //std::cout << "with norm:" << std::endl << std::flush;
  //std::cout << sqrt(vecU[0]*vecU[0]+vecU[1]*vecU[1]+vecU[2]*vecU[2]+vecU[3]*vecU[3]) << std::endl << std::flush;

  return SU2(vecU);
};


/// GridFunction has 4 entries for each link
/// set link value at edge Inum
auto SetLink = [](std::shared_ptr<GridFunction> gfU, const int& Inum, SU2<double> U, bool forward)
{
  FlatVector<> gfUvec = gfU->GetVector().FV<double>();
  Array<DofId> dnums;
  gfU->GetFESpace()->GetEdgeDofNrs (Inum, dnums);

  Vec<4,double> qUvec = U.GetCoeff(); 

  if (forward == false)
  {
      for (int a = 1; a <=3; a++)
      {
        qUvec[a] *= -1.;
      }
  }

  for(int i = 0; i <=3; i++)
  {
    gfUvec(dnums[i]) = qUvec[i];  
  }
};


/// useful links:
///https://github.com/NGSolve/ngsolve/blob/master/comp/python_comp.cpp
///https://github.com/NGSolve/ngsolve/blob/master/comp/gridfunction.cpp
///https://github.com/NGSolve/ngsolve/blob/master/fem/coefficient.cpp
///https://github.com/NGSolve/ngsolve/blob/master/fem/intrule.cpp
///https://docu.ngsolve.org/latest/i-tutorials/unit-9.2-C++Assemble/cppassembling.html?highlight=cppassembling

//GetElementTrafo
//IntegrationRule
//IntegrationPoint


/// WORKS!
//auto Evaluatef = [](MeshVars &meshvars, std::shared_ptr<GridFunction>& f, Vec<2, double> ElC)
auto Evaluatef = [](MeshVars &meshvars, std::shared_ptr<GridFunction>& f, size_t Elnum)
{
  // from __call__ in python_comp.cpp
  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  Vec<2, double> ElC = meshvars.ElBarycenter(Elnum);

  /// (global?) IntegrationPoint
  IntegrationPoint ip(ElC, 1.);

  ElementId ElId(VOL, Elnum);

  /// TODO: how large LocalHeap?
  //LocalHeap lh;
  size_t heap_size = sizeof(double)*100;
  LocalHeap lh(heap_size, "Evaluate", true);
  auto space = f->GetFESpace();
  auto evaluator = space->GetEvaluator();
  //std::cout << "evaluator->Dim(): " << evaluator->Dim() << std::endl;
  const FiniteElement & fel = space->GetFE(ElId, lh);
  Array<int> dnums(fel.GetNDof(), lh);
  space->GetDofNrs(ElId, dnums);
  auto & trafo = ma->GetTrafo(ElId, lh);
  Vector<> elvec(fel.GetNDof()*space->GetDimension());
  Vector<> values(evaluator->Dim());
  f->GetElementVector(dnums, elvec);

  // for (auto& value : elvec)
  // { 
  //   value = 0.;
  // }

  for (auto& value : values)
  { 
    value = 0.;
  }
  
  evaluator->Apply(fel, trafo(ip, lh), elvec, values, lh);
  return values;
};




//auto IntegratedualA = [](MeshVars &meshvars, std::shared_ptr<CoefficientFunction>& A, size_t glob_Enum, size_t glob_Elnum)
auto IntegratedualA = [](MeshVars &meshvars, std::shared_ptr<GridFunction>& A, size_t glob_Enum, size_t glob_Elnum)
{

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  ElementId ElId(VOL, glob_Elnum);

  /// Check if edge hangs at element
  Array<int> ElInums( ma->GetElEdges(ElId) );
  bool edge_at_el = false;
  /// iterate over all element interfaces
  for(int Inum : ElInums)
  {
    //std::cout << "Checking Edge " << Inum << std::endl;
    /// check if interface hangs at element
    if (Inum == glob_Enum)
    {
      edge_at_el = true;
      break;
    }
  
  }

  /// set edge number in reference triangle
  size_t ref_Enum = 0;
  if (edge_at_el == true)
  {
    /// this is a way to reverse-engineer the local edge number from the global one
    /// by iterating over all edges at an element
    /// the order in ElEnums reflects their respective local numbers
    //ma->GetElEdges(glob_Elnum, ElEnums);
    //Array<int> ElEnums = ma->GetElEdges(ElId, ElEnums);
    Array<int> ElEnums (ma->GetElEdges(ElId)) ;
    ref_Enum = 0;
    for (auto& ElEnum : ElEnums)
    {
      if (ElEnum == glob_Enum)
      {
        break;
      }
      ref_Enum++;
    }

  }
  /// if the interface does not belong to this element ...
  if (edge_at_el != true)
  {
    /// ... raise exception
    /// (or pick some edge)
    std::cout << "WARNING:" << std::endl;
    std::cout << "Edge " << glob_Enum << " NOT at element " << glob_Elnum << std::endl;
  }
  
  //std::cout << "Setting ref_Enum to " << ref_Enum << std::endl;


  /// integration rule order required
  size_t n = A->GetFESpace()->GetOrder();

  /// number of A components (only color - excluding spatial ones?)
  size_t n_A = A->GetNComponents();
  //const size_t n_A = 3;
  if(n_A == 0)
  {
    n_A = 1;
  }
  //std::cout << "order n: " << n << std::endl;
  //std::cout << "components n_A: " << n_A << std::endl;


  /// define reference element dual lattice points
  /// TODO: customize for different elements: trigs, quads, ...
  ELEMENT_TYPE El_type = ma->GetElType(ElId);
  /// order of edges consistent with ElementTopology
  /// See ElementTopology::GetNormals()
  
  /// default element typ is ET_TRIG
  Vec<3,double> ref_Elc_p{1./3., 1./3., 0.};
  //Vec<3, double> ref_Edgec_ps[3] = { {0., 0.5, 0.},  {0.5, 0.5, 0.},  {0.5, 0., 0.} };
  Vec<3, double> ref_Edgec_ps[4] = { {0.5, 0., 0.}, {0., 0.5, 0.}, {0.5, 0.5, 0.}, {0., 0., 0.}};

  if (El_type == ET_TRIG)
  {
    //std::cout << "Element " << glob_Elnum << " ET_TRIG" << std::endl;
  }
  else if (El_type == ET_QUAD)
  {
    //std::cout << "Element " << glob_Elnum << " ET_QUAD" << std::endl;
    ref_Elc_p = Vec<3,double>{0.5, 0.5, 0.};
    //ref_Edgec_ps = Vec<3, double>[4] { {0.5, 0., 0.}, {0.5, 1., 0.}, {0., 0.5, 0.},  {1., 0.5, 0.} };
    ref_Edgec_ps[0] = Vec<3, double>{0.5, 0., 0.};
    ref_Edgec_ps[1] = Vec<3, double>{0.5, 1., 0.};
    ref_Edgec_ps[2] = Vec<3, double>{0., 0.5, 0.};
    ref_Edgec_ps[3] = Vec<3, double>{1., 0.5, 0.};
  }
  
  
  /// define starting and stopping points of dual edge (on reference element)
  Vec<3,double> refpstart = ref_Elc_p;
  Vec<3,double> refpstop = ref_Edgec_ps[ref_Enum];

  //std::cout << "refpstart: " << refpstart << std::endl;
  //std::cout << "refpstop: " << refpstop << std::endl;

  /// evaluate the direction unit vector (reference element)
  // Vec<2,double> dirvec;
  // double ref_length = 0.;
  // for (size_t k = 0; k < 2; k++)
  // {
  //   dirvec[k] = (refpstop - refpstart)[k];
  //   ref_length += dirvec[k] * dirvec[k];
  // }
  // ref_length = sqrt(ref_length);
  // dirvec = (1./ref_length) * dirvec;
  //std::cout << "dirvec " << dirvec << std::endl;
  

  /// evaluate the direction unit vector (mapped element)
  Vec<2,double> dirvec;
  auto pstart = meshvars.ElBarycenter(glob_Elnum);
  auto pstop = meshvars.EBarycenter(glob_Enum);

  //std::cout << "pstart: " << pstart << std::endl;
  //std::cout << "pstop: " << pstop << std::endl;

  double length = 0.;
  for (size_t k = 0; k < 2; k++)
  {
    dirvec[k] = (pstop - pstart)[k];
    length += dirvec[k] *  dirvec[k];
  }
  length = sqrt(length);
  dirvec = (1./length) * dirvec;
  //std::cout << "dirvec " << dirvec << std::endl;
  //std::cout << "length " << length << std::endl;



  /// get segment integration rule for order n
  IntegrationRule seg_ir(ET_SEGM, n);

  /// define IntegrationRule (on reference triangle)
  /// IntegrationRule :: IntegrationRule (size_t asize, double (*pts)[3], double * weights)
  /// keep it empty for now:
  IntegrationRule dualedge_ir;

  /// add IntegrationPoints
  for (IntegrationPoint seg_ip : seg_ir)
  {
    /// obtain weight of segment integration point
    double weight = seg_ip.Weight();
    //std::cout << "weight: " << weight << std::endl;
    /// obtain distance s to the starting point
    //std::cout << "seg_ip: " << seg_ip << std::endl;
    double s = seg_ip.Point()[0];
    
    /// calculate point on the dual reference edge via a convex combination
    Vec<3,double> ipvec = (1.-s)*refpstart + s*refpstop;
    //std::cout << "ipvec: " << ipvec << std::endl;
    double ipcoord[3] = {ipvec[0], ipvec[1], ipvec[2]};

    IntegrationPoint ip(ipcoord, weight);

    /// add integration point
    /// INLINE void AddIntegrationPoint (const IntegrationPoint & ip)
    dualedge_ir.AddIntegrationPoint(ip);
  }

  /// (IntegrationPoints in __call__ function in python.comp)
  /// from __call__ in python_comp.cpp:
  /// TODO: how large LocalHeap?
  size_t heap_size = sizeof(double)*100;
  LocalHeap lh(heap_size, "IntegratedualA", true);
  auto space = A->GetFESpace();
  auto evaluator = space->GetEvaluator();
  //std::cout << "evaluator->Dim(): " << evaluator->Dim() << std::endl;
  const FiniteElement & fel = space->GetFE(ElId, lh);
  Array<int> dnums(fel.GetNDof(), lh);
  space->GetDofNrs(ElId, dnums);
  auto & trafo = ma->GetTrafo(ElId, lh);
  Vector<> elvec(fel.GetNDof()*space->GetDimension());
  Vector<> values(evaluator->Dim());
  A->GetElementVector(dnums, elvec);
  
  // for (auto& value : elvec)
  // { 
  //   value = 0.;
  // }

  for (auto& value : values)
  { 
    value = 0.;
  }

  /// get BaseMappedIntegrationRule = Integration rule on physical triangle
  //auto & dualedge_mir = trafo(dualedge_ir, lh);


  /// integrate A (including all space and color directions) using quadrature weights
  Vector<double> Acdirint(evaluator->Dim());

  for (auto & Acdir : Acdirint)
  {
    Acdir = 0.;
  }

  //std::cout << "number of IPs: " << dualedge_ir.GetNIP() << std::endl;
  for (int i = 0 ; i < dualedge_ir.GetNIP(); i++)
  {
    //MappedIntegrationPoint<2,2> mip(dualedge_ir[i], trafo);
    //Vec<2*n_A,double> value = A->Evaluate(mip);

    evaluator->Apply(fel, trafo(dualedge_ir[i], lh), elvec, values, lh);

    //auto value = A->Evaluate(mip);
    //std::cout << "typeof(value): "  << typeid(value) << std::endl;
    //std::cout << "integration point "  << i << ":" << std::endl;
    //std::cout << "reference point "  << dualedge_ir[i].Point() << ":" << std::endl;
    //std::cout << "weight "  << dualedge_ir[i].Weight() << ":" << std::endl;
    //std::cout << "values: "  << values << std::endl;
    //std::cout << "elvec: "  << elvec << std::endl;
    //std::cout << std::endl;

    Acdirint += dualedge_ir[i].Weight() * values;
  }

  //std::cout << "Acdirint: "  << Acdirint << std::endl;


  /// for each color component Ac
  /// take inner product with dual edge direction vector
  /// three components c for three colors
  Vec<3,double> Aintvec{0., 0., 0.};
  for(size_t c = 0; c < n_A; c++)
  {
    /// take the inner product of the A vector with the integration direction
    //double Acint = Acdirint[2*c]*dirvec[0] + Acdirint[2*c + 1]*dirvec[1];
    double Acint = 0;
    for (size_t k = 0; k < 2; k++)
    {
      Acint += Acdirint[2*c + k]*dirvec[k];
    }

    /// collect those values
    Aintvec[c] = Acint;
  }

  /// need to scale with real dual mesh length
  Aintvec = length*Aintvec;
  //std::cout << std::endl;
  //std::cout << "Aintvec: " << Aintvec << std::endl;

  /// return values (as vectors)
  return Aintvec;
};

/// TODO: use CoefficientFunction instead
//auto InterpolateAtoU = [&IntegratedualA, &SetLink](MeshVars &meshvars, std::shared_ptr<GridFunction>& A, std::shared_ptr<GridFunction>& gfU, const int & order)
auto InterpolateAtoU = [&IntegratedualA, &SetLink](MeshVars &meshvars, std::shared_ptr<GridFunction>& A, std::shared_ptr<GridFunction>& gfU)
{
  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  //std::cout << "Number of Elements: " << ma->GetNE() << std::endl;
  //std::cout << "Number of Edges: " << ma->GetNEdges() << std::endl;
  for(int Enum = 0; Enum < ma->GetNEdges(); Enum++)
  {

    SU2<double> expAEint{1., 0., 0., 0.};

    /// obtain elements hanging at edge
    Array<int> EElnums;
    ma->GetEdgeElements(Enum, EElnums);


    

    /// special case of border elements - only one link
    if(EElnums.Size() == 1)
    {
      su2<double> Aint0 = su2<double> ( IntegratedualA(meshvars, A, Enum, EElnums[0]) );
      if ( meshvars.ElLinkOrientation(Enum, EElnums[0]) )
      {
        expAEint = Exp(Aint0);
      }
      else
      {
        expAEint = Exp(Aint0.Conj());
      }
    }
    else if(EElnums.Size() == 2)
    {
      su2<double> Aint0 = su2<double> ( IntegratedualA(meshvars, A, Enum, EElnums[0]) );
      su2<double> Aint1 = su2<double> ( IntegratedualA(meshvars, A, Enum, EElnums[1]) );

      if ( meshvars.ElLinkOrientation(Enum, EElnums[0]) )
      {
        expAEint = Exp(Aint0) * Exp(Aint1.Conj());
      }
      else
      {
        expAEint = Exp(Aint1) * Exp(Aint0.Conj());
      }
    }
    else
    {
      std::cout << "WARNING: Edge" << Enum << " has " << EElnums.Size() << " elements"  << std::endl;
    }

    SetLink(gfU, Enum, expAEint, true);
  }
};

auto GetJFlux = [](std::shared_ptr<GridFunction> gfj, const int& Inum, bool forward)
{
  FlatVector<> gfjvec = gfj->GetVector().FV<double>();
  Array<DofId> dnums;
  gfj->GetFESpace()->GetEdgeDofNrs (Inum, dnums);
  Vec<3,double> vecj{gfjvec[dnums[0]], gfjvec[dnums[1]], gfjvec[dnums[2]]};
  if (forward == false)
  {
    //vecj *= -1;
    for (int a = 0; a < 3; a++)
    {
      vecj[a] *= -1.;
    }
  }

  //return vecj;
  return su2<double>(vecj);
};


/// TODO: find correct function to integrate rho over element
//auto GetRho = [](std::shared_ptr<GridFunction> gfrho, const int& Elnum, bool forward)
//{
//  return integrated charge over element
//  return Integrate(gfrho, mesh(Elnum));
//}
/// the above function should replace the one below

/// atm below: a vector entry is picked out
/// and integration is performed in python code
auto GetRho = [](std::shared_ptr<GridFunction> gfrhoint, const int& Elnum, bool forward)
{
  FlatVector<> gfvecrhoint = gfrhoint->GetVector().FV<double>();
  Array<DofId> dnums;
  //gfrhoint->GetFESpace()->GetDofNrs (Elnum, dnums);
  gfrhoint->GetFESpace()->GetDofNrs (ElementId(VOL, Elnum), dnums);
  //Array<DofId> dnums (gfrhoint->GetFESpace()->GetDofNrs (Elnum) );
  Vec<3,double> vecrhoint{gfvecrhoint[dnums[0]], gfvecrhoint[dnums[1]], gfvecrhoint[dnums[2]]};
  if (forward == false)
  {
    //vecrhoint *= -1;
    for (int a = 0; a < 3; a++)
    {
      vecrhoint[a] *= -1.;
    }
  }

  return su2<double>(vecrhoint);
};


/*
Calculate the Wilson loop around bone Bnum starting with gauge link at Interface Inum and with direction dir.
*/



auto GetWilsonLoop = [&ArrangeOrientedEdgesOfVertex, &GetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, const int& Bnum, int Inum, const bool& dir)
{
  SU2<double> W_B{1.,0.,0.,0.}; 

  // Bnum: bone number
  // 2D: vertex Vnum
  // 3D: edge

  // Inum: (Inter)Facet number
  // 2D: edge Enum
  // 3D: face

  //std::shared_ptr<MeshAccess> ma = gfU->GetMeshAccess();
  //std::shared_ptr<MeshAccess> ma = meshvars->GetMeshAccess();
  /// instead of recalculating edges hanging at vertex every time ...
  //std::vector<std::tuple<int, bool>> _I_ors = CalcVertexOrientedEnums(ma, Bnum);
  //std::vector<std::tuple<int, bool>> _I_ors = GetOrientedEdgesOfVertex(ma, Bnum);
  /// ... retrieve values from table/array
  std::vector<std::tuple<int, bool>> _I_ors = meshvars.GetOrientedEdgesOfVertex(Bnum);


  //std::cout << "before sorting" << std::endl << std::flush;
  //for(auto _I_or : _I_ors)
  //{
  //  std::cout << " Enum " << std::get<0>(_I_or) << " or " << std::get<1>(_I_or) << std::endl << std::flush;
  //}
  
  // ArrangeOrientedEdgesOfVertex(_I_ors, Bnum, Inum, dir);
  
  /// If the Wilson loop should start with the link at the interface Inum
  if(Inum >= 0)
  {
    /// then the sequence of interfaces needs to be cycled
    // such that the provided interface Inum is the first in order
  }
  else
  {
    /// if Inum < 0 is specified, an arbitrary starting link can be chosen
    /// still the direction needs to be flipped if dir==false
    Inum = std::get<0>(_I_ors[0]);
  }
  _I_ors = ArrangeOrientedEdgesOfVertex(_I_ors, Bnum, Inum, dir);  


  // //std::cout << "after sorting" << std::endl << std::flush;
  // for(auto _I_or : _I_ors)
  // {
  //   //std::cout << " Enum " << std::get<0>(_I_or) << " or " << std::get<1>(_I_or) << std::endl << std::flush;
  // }

  //std::cout << "starting multiplication loop" << std::endl << std::flush;
  /// pick out gauge links for each interface
  /// CAUTION: the variable _or takes the GLOBAL orientation into account (see GetOrientedEdgesOfVertex(int Bnum))
  for(auto _I_or : _I_ors)
  {
    auto _Inum = std::get<0>(_I_or);
    auto _or = std::get<1>(_I_or);

    SU2 U = GetLink(gfU, _Inum, _or);
    //std::cout << "U:" << std::endl << std::flush;
    //std::cout << U << std::endl << std::flush;

    /// LGT MULTIPLICATION CONVENTION
    /// and sequentially group multiply
    W_B = W_B*U;
    //W_B = U*W_B;
    //std::cout << "W_B:" << std::endl << std::flush;
    //std::cout << W_B << std::endl << std::flush;
  }

  /// return group element
  return W_B;
};



auto GettWilsonLoop = [&GetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, const int& Bnum, const bool& dir)
{

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  /// obtain elements hanging at edge
  Array<int> EElnums;
  ma->GetEdgeElements(Bnum, EElnums);

  int EElnum = EElnums[0];

  /// to evaluate orientation
  bool _or = dir;
  if (dir == true)
  {
    _or = meshvars.ElLinkOrientation(Bnum, EElnum);
  }
  else
  {
    _or = !( meshvars.ElLinkOrientation(Bnum, EElnum) );
  }

  SU2 Unew = GetLink(gfUnew, Bnum, _or);
  SU2 Uold = GetLink(gfUold, Bnum, !_or);

  SU2<double> W_B = Unew * Uold; 

  /// return group element
  return W_B;
};



/// The function that implements a timestep update
/// according to the spatial variation of the mass lumped simplicial wilson action
auto EOMLinkUpdate = [&GetLink, &GetWilsonLoop, &GetJFlux, &SetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfjhat, const int& Inum, const double& dt)
{
  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  /// get bones hanging at interface
  /// 2 vertices for each edge in 2D
  INT<2> EVnums = ma->GetEdgePNums(Inum);

  /// to get current flux (upwind!)
  /// evaluate u*n
  /// choose correct side of the element
  /// extract integration point from edge id
  /// get link values
  /// or: do this in python beforehand (this is done here)

  SU2<double> U = GetLink(gfU, Inum, true);
  SU2<double> Uoldinv = GetLink(gfUold, Inum, false);



  su2<double> jflux = GetJFlux(gfjhat, Inum, true);


  /// there are two spatial loops for each edge
  /// with opposite orientations
  //SU2<double> W_v0 = GetWilsonLoop(gfU, EVnums[0], Inum, true);
  /// this one goes first along Inum
  SU2<double> W_v0 = GetWilsonLoop(meshvars, gfU, EVnums[0], Inum, true);

  /// this one goes LAST along Inum
  /// therefore needs a minus sign in the update
  //SU2<double> W_v1 = GetWilsonLoop(gfU, EVnums[1], Inum, false);
  SU2<double> W_v1 = GetWilsonLoop(meshvars, gfU, EVnums[1], Inum, false);
  
  /// one old temporal loop for each edge
  SU2<double> oldP = U*Uoldinv;

  /// TODO: check signs and orientation, add correct weights and dt factors - is the overall formula correct?


  /// instead of fixing weights ...
  //double Ct0 = 0.00111111111;
  /// ... retrieve weights from table
  /// ~ inverse mesh mass*dt^2/area^2 essentially
  double Ct0 = meshvars.GettBoneWeight(size_t(EVnums[0]));

  //double Ct1 = tbone_weights[size_t(EVnums[1])];
  double Ct1 = meshvars.GettBoneWeight(size_t(EVnums[1]));
  
  /// ~ inverse mesh mass*length^2/dual_length^2 essentially
  //double Ce = bone_weights[Inum];
  double Ce = meshvars.GetBoneWeight(size_t(Inum));

  double Ceinv = 1./Ce;


  //Vec<3,double> newP_evec = Ct * ( (W_v0.Trta()) - (W_v1.Trta()) );
  //Vec<3,double> newP_evec =  (Ct0 *  (W_v0.Trta())) - (Ct1 * (W_v1.Trta()));
  
  Vec<3,double> newP_evec =  (Ct0 * (W_v0.Trta())) - (Ct1 * (W_v1.Trta())) + jflux.Trta();
  
  /// newP = Unew*Uinv
  
  /// difference to python code (=vec_from_q?) is -1/2
  /// CAUTION: variable holding new coefficients should be in the QUATERNION representation
  /// formula: quaternion coefficients W_a = -1/2 * Tr(t_a W) where W is the corresponding SU2 matrix
  /// therefore the factor -1/2 is needed in front of Trta
  /// to obtain quaternion coefficients of the new Polyakov loop newP
  /// (see SU2 implementation of GetCoeff)
  newP_evec = -0.5 * newP_evec;


  /// negative prefactor for space-time metric in Ct weights!
  newP_evec = -1.*newP_evec;

  /// not necessary if dt is already included in Meshvars weights
  //newP_evec = dt*dt*newP_evec;

  newP_evec = Ceinv*newP_evec;
  

  Vec<3,double> oldP_evec =  -0.5 * oldP.Trta();
  newP_evec += oldP_evec;

  //std::cout << "update vec:" << std::endl << std::flush;
  //std::cout << newP_evec << std::endl << std::flush;

  double newP_e0 = 0.;
  double newP_evec_norm2 = newP_evec[0]*newP_evec[0] + newP_evec[1]*newP_evec[1] + newP_evec[2]*newP_evec[2];
  
  //std::cout << "with norm:" << std::endl << std::flush;
  //std::cout << newP_evec_norm2 << std::endl << std::flush;
  //std::cout << std::endl << std::flush;

  if (newP_evec_norm2 < 1.)
  {
    newP_e0 = sqrt(1. - newP_evec_norm2);
  }
  else
  {
    std::cout << "invalid SU2 values with norm= " << newP_evec_norm2 << std::endl << std::flush;
    std::cout << newP_evec << std::endl << std::flush;
    newP_e0 = 0.;
    std::cout << "capping u0 to " << newP_e0 << std::endl << std::flush;
  }
  //std::cout << newP_e0 << std::endl << std::flush;


  /// LGT MULTIPLICATION CONVENTION
  SU2<double> newP(newP_e0, newP_evec[0], newP_evec[1], newP_evec[2]);
  SU2<double> newU = newP * U;

  SetLink(gfUnew, Inum, newU, true);
};

/// global EOM update
/// currently fixed boundary conditions
auto EOMUpdate = [&GetWilsonLoop, &EOMLinkUpdate](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfjhat, const double& dt)
{
  /// iterate over all interfaces
  /// edges in 2D
  /// and perform the EOM update

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  /// can be parallelized
  for(int Inum = 0; Inum < ma->GetNEdges(); Inum++)
  {
    //std::cout << "edge: " << Inum << std::endl;
    if (meshvars.FacetVorb(Inum) == true)
    //if (true)
    {
      //std::cout << Inum << " is in VOL" << std::endl;
      EOMLinkUpdate(meshvars, gfU, gfUold, gfUnew, gfjhat, Inum, dt);  
    }
    else
    {
      //std::cout << Inum << " is in BND" << std::endl;
    }
  }
};




auto GaussElementVecError = [&GetRho, &GetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfrhoint, const int& Elnum, const double& eps)
{
  //std::cout << "GaussElementCheck for Element: " << Elnum << std::endl;
  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  SU2<double> rho_el = GetRho(gfrhoint, Elnum, true);

  /// get interfaces hanging at elements
  //Array<int> ElInums;
  //ma->GetElEdges(ElementId(VOL, Elnum), ElInums);
  
  Array<int> ElInums( ma->GetElEdges(ElementId(VOL,Elnum)) );

  Vec<3,double> gauss_flux{0., 0., 0.};

  /// iterate over all element interfaces
  for(int Inum : ElInums)
  {
    /// evaluate the contribution of this interface to the overall gauss flux out of this element

    //std::cout << "interface number: " << Inum << std::endl;
    double Ce = meshvars.GetBoneWeight(size_t(Inum));
    //std::cout << "bone weight: " << Ce << std::endl;
  


    SU2<double> Unew_el{1., 0., 0., 0.};
    SU2<double> Uold_el{1., 0., 0., 0.};
    
    /// flip Us if necessary
    /// the Wilson loops for each edge involved in the gauss constraint
    /// have a specific orientation with respect to the element
    bool ElLinkOr = meshvars.ElLinkOrientation(Inum, Elnum);
    if ( ElLinkOr == true )
    {
      //std::cout << "not flipping W" << std::endl;
      Unew_el = GetLink(gfUnew, Inum, true);
      Uold_el = GetLink(gfUold, Inum, false);
    }
    else
    {
      //std::cout << "flipping W" << std::endl;
      Unew_el = GetLink(gfUnew, Inum, false);
      Uold_el = GetLink(gfUold, Inum, true);
    }

    /// temporal Wilson loop (~ exp{\int E dt})
    SU2<double> W = Unew_el*Uold_el;

    /// add the contribution
    gauss_flux += Ce * W.Trta();
  }

  /// subtract rho sources from the field flux gauss_flux
  Vec<3,double> gauss_diff = 2 * gauss_flux  - rho_el.Trta();

  return gauss_diff;
};



auto GaussElementCheck = [&GetRho, &GetLink, &GaussElementVecError](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfrhoint, const int& Elnum, const double& eps)
{
  Vec<3,double> gauss_diff = GaussElementVecError(meshvars, gfUold, gfUnew, gfrhoint, Elnum, eps);

  double gauss_err2 = 0.;
  for(int c = 0; c < 3; c++)
  {
    gauss_err2 += gauss_diff[c]*gauss_diff[c];
  }

  // if(gauss_err2 < eps*eps)
  // {
  //   return true;
  // }
  // else
  // {
  //   return false;
  // }

  return gauss_err2;
};


/// global Gauss check
auto GaussCheck = [&GaussElementCheck](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfrhoint, const double& eps)
{
  /// iterate over all interfaces
  /// edges in 2D
  /// and perform the Gauss update

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  double gauss_err2 = 0.;
  //std::cout << "Number of Elements: " << ma->GetNE() << std::endl;
  for(int Elnum = 0; Elnum < ma->GetNE(); Elnum++)
  {
    //std::cout << Elnum << std::endl;
    gauss_err2 += GaussElementCheck(meshvars, gfUold, gfUnew, gfrhoint, Elnum, eps);
    //std::cout << "gauss constraint error: " << gauss_err2 << std::endl;
    if (gauss_err2 > eps*eps)
    {
      //std::cout << "element: " << Elnum << " gauss_err2: " << gauss_err2 << std::endl;
      //return false;
    }
  }

  //return true;
  return gauss_err2;

};




/// an iterative method of how to generate two field configurations that fulfill the gauss constraint
/// by variation of the gauss constraint violation
/// and a group multiplicative update
auto GaussLinkUpdate = [&GetLink, &SetLink, &GetRho, &GaussElementVecError](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfrhoint, const int& Inum, const double& eps)
{
  //std::cout << "GausslinkUpdate for Edge: " << Inum << std::endl;
  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  double Ce = meshvars.GetBoneWeight(size_t(Inum));
  //std::cout << "bone weight: " << Ce << std::endl;

  /// get elements hanging at interface
  /// 2 triangles for each edge in 2D
  Array<int> EElnums;
  ma->GetEdgeElements(Inum, EElnums);
  //std::cout << "Edge Elements: " << EElnums << std::endl;

  su2<double> var_dir[3];
  var_dir[0] = su2<double>{1.,0.,0.};
  var_dir[1] = su2<double>{0.,1.,0.};
  var_dir[2] = su2<double>{0.,0.,1.};


  /// the updated link value
  SU2<double> Uupdated = GetLink(gfUnew, Inum, true);

  for(int var_i = 0; var_i < 3; var_i++)
  {
    /// variation of gauß constraint violation
    double var_err = 0.;

    /// go through all elements connected to interface
    // the element loop is necessary to ensure proper W orientations (always out of element!) in variation terms
    for(int Elnum : EElnums)
    {

      //std::cout << "Edge Element: " << Elnum << std::endl;

      SU2<double> Unew_el{1., 0., 0., 0.};
      SU2<double> Uold_el{1., 0., 0., 0.};

      /// flip Us if necessary
      /// the Wilson loops for each edge involved in the gauss constraint
      /// have a specific orientation with respect to the element
      
      bool ElLinkOr = meshvars.ElLinkOrientation(Inum, Elnum);
      if ( ElLinkOr == true )
      {
        //std::cout << "not flipping W" << std::endl;
        Unew_el = GetLink(gfUnew, Inum, true);
        Uold_el = GetLink(gfUold, Inum, false);
      }
      else
      {
        //std::cout << "flipping W" << std::endl;
        Unew_el = GetLink(gfUnew, Inum, false);
        Uold_el = GetLink(gfUold, Inum, true);
      }

      su2<double> rho_el = GetRho(gfrhoint, Elnum, true);
      //std::cout << "rho_el: " << std::endl << rho_el << std::endl;


      Vec<3,double> var_W{0., 0., 0.};
      /// take care of relative directions when varying a loop with respect to the reversly oriented parallel transporter
      if (ElLinkOr)
      {
        var_W = (SU2<double>(var_dir[var_i].GetMatrix()) * Unew_el * Uold_el).Trta();
      }
      else
      {
        var_W = -1.* ( Unew_el * (SU2<double>(var_dir[var_i].GetMatrix())) * Uold_el).Trta();
      }
      //std::cout << "var_W: " << var_W << std::endl;

      // calculate variation of gauss constraint violation
      Vec<3,double> gauss_diff = GaussElementVecError(meshvars, gfUold, gfUnew, gfrhoint, Elnum, eps);
      //std::cout << "gauss_diff: " << gauss_diff << std::endl;

      for(int c = 0; c < 3; c++)
      {
        var_err += var_W[c] * gauss_diff[c];
      }
    }
    var_err = 4*Ce *var_err;
    //std::cout << "var_err " << var_err << std::endl;

    // calculate update element
    su2<double> update_vec = -1. * eps * var_err * var_dir[var_i];
    //su2<double> update_vec = eps * var_err * var_dir[var_i];
    //std::cout << "update_vec: " << update_vec << std::endl;

    /// constructor uses exponential map
    SU2<double> update(update_vec);

    //std::cout << "update: " << std::endl << update << std::endl;
    // reconstruct and save new gauge link
    //std::cout << "U before update: " << std::endl << Uupdated << std::endl;
    Uupdated = update * Uupdated;
    //std::cout << "U after update: " << std::endl << Uupdated << std::endl;
  }

  SetLink(gfUnew, Inum, Uupdated, true);
};


/// global Gauss update
auto GaussUpdate = [&GaussLinkUpdate](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfrhoint, const double& eps)
{
  /// iterate over all interfaces
  /// edges in 2D
  /// and perform the Gauss update

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  //std::cout << "Number of edges: " << ma->GetNEdges() << std::endl;

  /// can be parallelized
  for(int Inum = 0; Inum < ma->GetNEdges(); Inum++)
  {
    //std::cout << "edge: " << Inum << std::endl;
    if (meshvars.FacetVorb(Inum) == true)
    {
      //std::cout << "is in VOL" << std::endl;
      GaussLinkUpdate(meshvars, gfUold, gfUnew, gfrhoint, Inum, eps);  
    }
    else
    {
      //std::cout << "is in BND" << std::endl;
    }
    
  }

};





auto Calc_ln_normg = [](SU2<double> W)
{

  SU2<double> Id{1., 0., 0., 0.};
  //std::cout << Id << std::endl << std::flush;
  
  Matrix<MyComplex<double>> lnWmat{2,2};
  lnWmat = Id.GetMatrix() - W.GetMatrix();

  /// the matrices below are strictly speaking not SU2 ...
  /// ... but quaternions
  /// alternatively use a logarithm?
  SU2<double> lnW(lnWmat);
  //std::cout << "lnW" << lnW << std::endl << std::flush;
  //std::cout << "lnWconj"<< lnW.Conj() << std::endl << std::flush;
  //std::cout << "lnW2" << lnW * (lnW.Conj()) << std::endl << std::flush;
  //std::cout << "Ct"  << Ct << std::endl << std::flush;
  //std::cout << "S for vnum=" << Vnum << ":" << Ct * (lnW * (lnW.Conj()) ).ReTr() << std::endl << std::flush;
  return (lnW * (lnW.Conj()) ).ReTr();

  /// TODO: is an action with a log an option?
  /// variation and all theory needs to be revisited
  //su2<double> A_W = Log(A_W);
  //return gInnerProduct(A_W, A_W)
};


/// evaluate representative of the gauge dependent magnetic field
auto CalcgfBW = [&GetWilsonLoop, &SetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, std::shared_ptr<GridFunction> gfBW)
{
  /// iterate over all temporal bones
  /// vertices in 2D
  /// and calculate curvature/Wilson loops

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  for(int Vnum = 0; Vnum < ma->GetNV(); Vnum++)
  {

    SU2<double> BW = GetWilsonLoop(meshvars, gfU, Vnum, -1, true);
    /// TODO: get pure Hodge weights?
    SetLink(gfBW, Vnum, BW, true);
  }

};

/// evaluate representative of the gauge dependent electric field
auto CalcgfEW = [&GettWilsonLoop, &SetLink](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfEW)
{
  /// iterate over all temporal bones
  /// vertices in 2D
  /// and calculate curvature/Wilson loops

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();

  for(int Enum = 0; Enum < ma->GetNEdges(); Enum++)
  {

    SU2<double> EW = GettWilsonLoop(meshvars, gfUold, gfUnew, Enum, true);
    /// TODO: get pure Hodge weight?
    SetLink(gfEW, Enum, EW, true);
  }

};


auto CalcHBV = [&GetWilsonLoop, &Calc_ln_normg](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, int & Vnum)
{
  double Ct = meshvars.GettBoneWeight(size_t(Vnum));
  SU2<double> W = GetWilsonLoop(meshvars, gfU, Vnum, -1, true);
  return Ct * Calc_ln_normg(W);

};

auto CalcHEE = [&GettWilsonLoop, Calc_ln_normg](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, int & Enum)
{
  double Ce= meshvars.GetBoneWeight(size_t(Enum));
  SU2<double> W = GettWilsonLoop(meshvars, gfUold, gfUnew, Enum, true);
  return Ce * Calc_ln_normg(W);

};

auto CalcgfHB = [&GetWilsonLoop, &CalcHBV](MeshVars& meshvars, std::shared_ptr<GridFunction> gfU, std::shared_ptr<GridFunction> gfHB)
{
  /// iterate over all temporal bones
  /// vertices in 2D
  /// and calculate curvature/Wilson loops

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();
  FlatVector<> gfHBvec = gfHB->GetVector().FV<double>();

  for(int Vnum = 0; Vnum < ma->GetNV(); Vnum++)
  {

    double HBV = CalcHBV(meshvars, gfU, Vnum);
    //std::cout << "S for vnum=" << Vnum << ":" << HBV << std::endl << std::flush;

    gfHBvec(Vnum) = HBV;
  }

};


auto CalcgfHE = [&GetWilsonLoop, &CalcHEE](MeshVars& meshvars, std::shared_ptr<GridFunction> gfUold, std::shared_ptr<GridFunction> gfUnew, std::shared_ptr<GridFunction> gfHE)
{
  /// iterate over all  bones
  /// edges in 2D
  /// and calculate curvature/Wilson loops

  std::shared_ptr<MeshAccess> ma = meshvars.GetMeshAccess();
  FlatVector<> gfHEvec = gfHE->GetVector().FV<double>();

  for(int Enum = 0; Enum < ma->GetNEdges(); Enum++)
  {
    double HEE = CalcHEE(meshvars, gfUold, gfUnew, Enum);
    //std::cout << "S for vnum=" << Vnum << ":" << HEE << std::endl << std::flush;

    gfHEvec(Enum) = HEE;
  }

};


m.def("ArrangeOrientedEdgesOfVertex", ArrangeOrientedEdgesOfVertex);
m.def("GetWilsonLoop", GetWilsonLoop);
m.def("EOMLinkUpdate", EOMLinkUpdate);
m.def("GaussLinkUpdate", GaussLinkUpdate);
m.def("GaussElementCheck", GaussElementCheck);
m.def("EOMUpdate", EOMUpdate);
m.def("GaussUpdate", GaussUpdate);
m.def("GaussCheck", GaussCheck);
m.def("CalcgfHB", CalcgfHB);
m.def("CalcgfHE", CalcgfHE);
m.def("CalcHBV", CalcHBV);
m.def("CalcHEE", CalcHEE);
m.def("CalcgfBW", CalcgfBW);
m.def("CalcgfEW", CalcgfEW);
m.def("Calc_ln_normg", Calc_ln_normg);
m.def("GetLink", GetLink);
m.def("SetLink", SetLink);
m.def("IntegratedualA", IntegratedualA);
m.def("Evaluatef", Evaluatef);
m.def("InterpolateAtoU", InterpolateAtoU);

//m.attr("tbone_weights") = &tbone_weights;

// END USER DEFINED CODE

res = m;    

}
