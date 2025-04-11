#include <python_comp.hpp>


#include "quaternion.hpp"
#include "liegroups.hpp"
#include "liegroupfe.hpp"

#include "SO3coefficient.hpp"

using namespace ngsolve;

PYBIND11_MODULE(liegroupfem,m) {
  auto ngs = py::module::import("ngsolve");

/*
  py::class_<Vec<3,double>> (m, "qVec")
    .def(py::init<Vec<3,double>>())
    ;
*/
  py::class_<Quaternion<>>
     (m, "Quaternion")
    .def(py::init([](double s, Vec<3,double> v) { return Quaternion<>(s, v); }))
    .def(py::init([](double s, std::array<double,3> v) { return Quaternion<>(s, {v[0], v[1], v[2]}); }))
    
    .def_property_readonly ("scal", &Quaternion<double>::getScal)
    .def_property_readonly ("vec", &Quaternion<double>::getVec)
    
    .def("__repr__", [] (Quaternion<> q) { return ToString(q); })
    .def("__add__", [](Quaternion<> q1, Quaternion<> q2) { return q1+q2; })
    .def("__mul__", [](Quaternion<> q1, Quaternion<> q2) { return q1*q2; })
    .def("conjugate", [](Quaternion<> q) { return q.conjugate(); })
    ;

  py::class_<UnitQuaternion<>>
    (m, "UnitQuaternion")
    .def(py::init([](double s, Vec<3,double> v) { return UnitQuaternion<>(s, v); }))
    .def(py::init([](double s, std::array<double,3> v) { return UnitQuaternion<>(s, {v[0], v[1], v[2]}); }))

    .def_property_readonly ("scal", &UnitQuaternion<double>::getScal)
    .def_property_readonly ("vec", &UnitQuaternion<double>::getVec)

    .def("__repr__", [] (UnitQuaternion<> q) { return ToString(q); })
    .def("__add__", [](UnitQuaternion<> q1, UnitQuaternion<> q2) { return q1+q2; })
    .def("__mul__", [](UnitQuaternion<> q1, UnitQuaternion<> q2) { return q1*q2; })
    .def("conjugate", [](UnitQuaternion<> q) { return q.conjugate(); })
;

m.def("GetRotQuaternion", [](double angle, std::array<double,3> axis) { return ngsolve::GetRotQuaternion(angle, {axis[0], axis[1], axis[2]}); });

  py::class_<SO3<>> (m, "SO3")
    .def(py::init<UnitQuaternion<>>())
    .def("__repr__", [] (SO3<> rot) { return ToString(rot); })
    .def("__mul__", [](SO3<> rot, Vec<3,double> v) { return rot.ApplyMatrix(v); })
    .def("getq", [](SO3<> rot) { return rot.getq(); })
    .def("GetMatrix", [](SO3<> rot) { return rot.GetMatrix(); })
    .def("Inverse", [](SO3<> rot) { return rot.Inverse(); })
    ;

  py::class_<so3<>> (m, "so3")
    .def(py::init<Vec<3,double>>())
    .def(py::init<double, double, double>())
    .def("__repr__", [] (so3<> sigma) { return ToString(sigma); })
    .def("__mul__", [](so3<> sigma, Vec<3,double> v) { return sigma.ApplyMatrix(v); })
    .def("getvec", [](so3<> sigma) { return sigma.getvec(); })
    .def("GetMatrix", [](so3<> sigma) { return sigma.GetMatrix(); })
  ;
  
  // py::class_<SO3CoefficientFunction> (m, "SO3CoefficientFunction")
  //   .def(py::init<>())
  // ;

m.def("Log", [](SO3<> rot) { return Log(rot); });
m.def("Exp", [](so3<> sigma) { return Exp( sigma ); });
m.def("Exp", [](std::array<double,3> v) { return Exp( so3(Vec<3,double>{v[0], v[1], v[2]}) ); });

m.def("Interpolate", [](double t, SO3<> a, SO3<> b) { return Interpolate(t, a, b); });
m.def("SO3CF", [] (shared_ptr<CoefficientFunction> cf) { return SO3CF(cf); });
m.def("VecCrossMatCF", [] (shared_ptr<CoefficientFunction> cf) { return VecCrossMatCF(cf); });

}