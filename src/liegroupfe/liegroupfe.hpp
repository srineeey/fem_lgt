#ifndef FILE_LIEGROUPFE
#define FILE_LIEGROUPFE

#include <fem.hpp>
#include "liegroups.hpp"


namespace ngsolve
{
  using namespace ngfem;


  template <typename G = SO3<double>>
  class LieGroupSegm : public FiniteElement
  {
  public:
    LieGroupSegm () : FiniteElement(2, 1) { }

    G Evaluate (IntegrationPoint & ip, std::vector<G> values)
    {
      double x = ip(0);
      G g0 = values[0];
      G g1 = values[1];
      /*
        for x = 0 --> g1
        for x = 1 --> g0
       */
      return Interpolate (x, g1, g0); 
    }
  };

  
}


#endif
