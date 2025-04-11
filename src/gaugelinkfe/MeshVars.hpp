#include <comp.hpp>
#include "comp.hpp"

#include <vector>


using namespace ngcomp;



/// Define BilinearFormIntegrator for H1 mass matrix - NOT NEEDED RIGHT NOW
// class H1HatMassIntegrator : public BilinearFormIntegrator
// {
//     public:
//       H1HatMassIntegrator (){ ; }
//       string Name () const override { return "MyH1HatMass"; }
//       int DimElement () const override { return 2; }
//       int DimSpace () const override { return 2; }
//       xbool IsSymmetric () const override { return true; }
//       VorB VB() const override { return VOL; }
//       // Calculates the element matrix
//       // ngbla::FlatMatrix<double, ngbla::RowMajor> ??
//       void CalcElementMatrix (const FiniteElement & base_fel,
//                               const ElementTransformation & eltrans, 
//                               FlatMatrix<double> elmat,
//                               LocalHeap & lh) const override
//       {
//         const ELEMENT_TYPE ETType = base_fel.ElementType();
//         /*
//           tell the compiler that we are expecting one of our elements.
//           if not, an exception will be raised
//         */
//         auto & fel = dynamic_cast<H1HighOrderFE<ETType> &> (base_fel);
//         //Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);
//         // number of element basis functions:
//         int ndof = fel.GetNDof();
//         elmat = 0;
//         Vector<> shape_ref(ndof); // basis on reference element
//         Vector<> shape(ndof);     // basis on mapped element
//         /*
//           get integration rule for element geometry, 
//           integration order is 2 times element order
//         */
//         IntegrationRule ir(fel.ElementType(), 2*fel.Order());
//         // loop over integration points
//         for (int i = 0 ; i < ir.GetNIP(); i++)
//           {
//             // calculate Jacobi matrix in the integration point
//             MappedIntegrationPoint<2,2> mip(ir[i], eltrans);
//             fel.CalcShape (ir[i], shape_ref);         
//             // transform it for the mapped element
//             shape = shape_ref;          
//             // integration weight and Jacobi determinant
//             double fac = mip.IP().Weight() * mip.GetMeasure();
//             // elmat_{i,j} += (fac*lam) * InnerProduct (grad shape_i, grad shape_j)
//             elmat += (fac) * shape * Trans(shape);
//           }     
//       }
// };





/// Various mesh variables and functions for lattice gauge theoretic calculations
class MeshVars
{
  std::shared_ptr<MeshAccess> ma = nullptr;

  /// keeps track of whether facets are inside volume or boundary
  std::unordered_map<size_t, bool> facet_vorb;

  /// keep track of the timestamp when facet_vorb was last evaluated
  /// TODO: currently not used - need to implement this
  std::unordered_map<size_t, size_t> facet_vorb_timestamp;

  /// keeps track of whether facets are inside volume or boundary
  std::unordered_map<size_t, bool> vertex_vorb;

  /// keep track of the timestamp when facet_vorb was last evaluated
  /// TODO: currently not used - need to implement this
  std::unordered_map<size_t, size_t> vertex_vorb_timestamp;

  /// save the weights for each bone in a list
  /// caution: this only takes care of the spatial weight
  /// temporal weights must also be added
  /// right now only spatially static meshes work
  std::unordered_map<size_t, double> tbone_weights;
  std::unordered_map<size_t, double> bone_weights;

  std::unordered_map<size_t, double> tbonemass;
  std::unordered_map<size_t, double> bonemass;

  double dt = 0.;

  /// keep track fo the timestamp when the weights was last evaluated
  std::unordered_map<size_t, size_t> tbone_weights_timestamp;
  std::unordered_map<size_t, size_t> bone_weights_timestamp;

  /// save the edges hanging at a vertex in an oriented fashion
  /// required for evaluating Wilson loops - multiplication order matters here!
  std::unordered_map<size_t, std::vector<std::tuple<int, bool>>> oriented_edges_of_vertex;

  /// keep track of the timestamp when the edge vertex map was last evaluated
  std::unordered_map<size_t, size_t> oriented_edges_of_vertex_timestamp;


  /// TODO: save global orientation gridfunction here
  //std::shared_ptr<GridFunction> gfor;

  public:

  MeshVars(std::shared_ptr<MeshAccess> _ma, double _dt) : ma(_ma), dt(_dt)
  {
    /// Refresh members
    EvaluateFacetVorb();
    EvaluateVertexVorb();

    //CalcBoneMass();
    //CalctBoneMass();

    RefreshBoneWeights();
    RefreshtBoneWeights();
    RefreshOrientedEdgesOfVertex();
  };

  MeshVars(std::shared_ptr<MeshAccess> _ma) : ma(_ma), dt(1.)
  {
    /// Refresh members
    EvaluateFacetVorb();
    EvaluateVertexVorb();

    RefreshBoneWeights();
    RefreshtBoneWeights();
    RefreshOrientedEdgesOfVertex();
  };

  std::shared_ptr<MeshAccess> GetMeshAccess()
  {
    return ma;
  }

  double get_dt()
  {
    return this->dt;
  }

  void set_dt(double _dt)
  {
    this->dt = _dt;
  }


  /// refresh facet_vorb variable by iterating through each edge
  void EvaluateFacetVorb()
  {
    auto ma = this->GetMeshAccess(); 
    for(int Inum = 0; Inum < ma->GetNEdges(); Inum++)
    {
      Array<int> EElnums;
      ma->GetEdgeElements(Inum, EElnums);
      if(EElnums.Size() > 1)
      { 
        //std::cout << Inum << " " << "VOL" << std::endl;
        this->facet_vorb[Inum] = true;
      }
      else
      {
        //std::cout << Inum << " " << "BND" << std::endl;
        this->facet_vorb[Inum] = false;
      }
      //std::cout << Inum << " " << facet_vorb[Inum] << std::endl;
    }
  }

  bool FacetVorb(size_t Inum)
  {
    return this->facet_vorb[Inum];
  }

  std::unordered_map<size_t, bool> GetFacetVorb(size_t Inum)
  {
    return this->facet_vorb;
  }



  /// refresh facet_vorb variable by iterating through each edge
  void EvaluateVertexVorb()
  {
    auto ma = this->GetMeshAccess(); 
    for(int Vnum = 0; Vnum < ma->GetNEdges(); Vnum++)
    {
      //Array<int> VElnums; 
      //ma->GetVertexElements(Vnum, VElnums);
      auto VElnums( ma->GetVertexElements(Vnum) );
      //VElnums[0] = _VElnums[0];
      //VElnums[1] = _VElnums[1];

      //Array<int> VElnums = ma->GetVertexElements(Vnum);
      if(VElnums.Size() > 2)
      { 
        //std::cout << Inum << " " << "VOL" << std::endl;
        this->vertex_vorb[Vnum] = true;
      }
      else
      {
        //std::cout << Inum << " " << "BND" << std::endl;
        this->vertex_vorb[Vnum] = false;
      }
      //std::cout << Inum << " " << facet_vorb[Inum] << std::endl;
    }
  }

  bool VertexVorb(size_t Vnum)
  {
    return this->vertex_vorb[Vnum];
  }

  std::unordered_map<size_t, bool> GetVertexVorb()
  {
    return this->vertex_vorb;
  }


  double CalcMassDiagEntry(const int & tBnum)
  {
    /// TODO: get Entry with mat(i,j) or mat[i,j]
    //return massbfmat(tBnum, tBnum);
    return 0;
  }

  /// TODO: calculate mass matrix entry for weights
  void CalctBoneMass()
  {
    auto ma = this->GetMeshAccess(); 
    /// TODO: get H1space of order 1
    //shared_ptr<H1HighOrderFESpace> fesh1(1);
    Flags flags;
    flags.SetFlag ("order", int(1) );
    
    //shared_ptr<H1HighOrderFESpace> fesh1_ptr = make_shared<H1HighOrderFESpace> (fesh1);
    auto fesh1 = make_shared<H1HighOrderFESpace>(ma, flags);
    
    auto u = fesh1->GetTrialFunction();
    auto up = fesh1->GetTestFunction();

    /// TODO: define bilinear form for mass matrix
    Flags flags_massbf;
    auto massbf = make_shared<T_BilinearFormSymmetric<double>> (fesh1, "massh1", flags_massbf);

    /// TODO: add mass term (BilinearFormIntegrator)
    massbf->AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(u * up, VOL, VOL));


    /// TODO: (Re)Assemble()
    size_t heap_size = sizeof(double)*ma->GetNV()*ma->GetNV()*5;
    LocalHeap lh(heap_size, "massbf", true);
    massbf->Assemble(lh);


    /// TODO: GetMatrix() - fix compiler error!
    SparseMatrix<double> tbonemassmat = massbf->GetMatrix();
    //auto& tbonemassmat = massbf->GetMatrix();
    
    /// tBnum = row number
    for(size_t tBnum = 0; tBnum < ma->GetNV(); tBnum++)
    {
        //std::cout <<  "tbone " << tBnum << std::endl << std::flush;
        /// TODO: find way to index matrix object
        /// j = number of non-zero entry in row
        // for(size_t j = tbonemassmat.GetRowIndices[tBnum]; j < tbonemassmat.GetRowIndices[tBnum+1]; j++)
        // {
        //   if (tbonemassmat.col[j] == tBnum)
        //   {
        //     this->tbonemass[tBnum] = tbonemassmat.GetRowValues(j);
        //   }
        // }
        

        this->tbonemass[tBnum] = tbonemassmat(tBnum, tBnum);
        //this->tbonemass[tBnum] = tbonemassmat(tBnum, tBnum);
    }
  }

  /// TODO: fix imbalance between volume and boundary weights in the case of general geometry
  double CalctBoneWeight(const int & tBnum)
  {
    /// calculate weight of the temporal bone (inverse area of barycentric dual vertex patch)
    double weight = 0.;

    //Array<int> Elnums;
    //ma->GetVertexElements(tBnum, Elnums);
    Array<int> Elnums( ma->GetVertexElements(tBnum) );
    
    Vec<2, double> b = ma->GetPoint<2>(size_t(tBnum));

    double dual_area = 0;
    for(int Elnum : Elnums)
    {
      //Array<int> ElVnums;
      //ma->GetElVertices(Elnum, ElVnums);
      //ma->GetElVertices(ElementId(VOL, Elnum), ElVnums);
      //Array<int> ElVnums( (ma->GetElVertices(ElementId(VOL, Elnum))) );
      //Array<int> ElVnums = (ma->GetElVertices(ElementId(VOL, Elnum)));

      /// get barycenter of element 
      //Vec<2, double> elc{0., 0.};
      Vec<2, double> elc = ElBarycenter(Elnum);

      /// vector from bone (vertex) to element center
      Vec<2, double> l_b_elc = elc - b;

      /// get lengths from bone to edge barycenters
      std::vector<Vec<2, double>> l_b_ecs;
      for (int Enum : ma->GetElEdges(Elnum))
      {
        INT<2> EVnums = ma->GetEdgePNums(Enum);

        if(tBnum == EVnums[0] || tBnum == EVnums[1])
        {
          Vec<2, double> l_b_ec = 0.5*( ma->GetPoint<2>(size_t(EVnums[0])) - ma->GetPoint<2>(size_t(EVnums[1])) );
          l_b_ecs.push_back(l_b_ec);
        }
      }
      if (l_b_ecs.size() != 2)
      {
        std::cout << l_b_ecs.size() <<  " edges associated to vertex " <<  tBnum <<" in element" << Elnum << std::endl << std::flush;
      }

      /// 0.5 * (c - b) x 0.5 * (v - b)
      for (Vec<2, double> l_b_ec: l_b_ecs)
      {
        dual_area += 0.5 * fabs( l_b_elc[0]*l_b_ec[1] - l_b_elc[1]*l_b_ec[0] ) ;
      }

    }

    //std::cout << "tBnum: " << tBnum << std::endl << std::flush;
    //std::cout << "dual_area: " << dual_area << std::endl << std::flush;

    /// The vertices at the border have distorted weights (less dual_area in the mesh)
    /// fix this for each border vertex
    if(!VertexVorb(tBnum))
    {
      if(Elnums.Size() == 2)
      {
        dual_area = 2.*dual_area;
      }
      else if(Elnums.Size() == 1)
      {
        dual_area = 4.*dual_area;
      }
    }

    if (dual_area*dual_area > 0.000000001)
    {
      //weight = dual_area;
      //weight = 1./dual_area;  
      //weight = 1./(dual_area*dual_area);  
      weight = dt*dt/(dual_area*dual_area);
      /// TODO: 3D dual cell factor? = mass matrix factor?
      //weight *= 1./3.;  
    }
    else
    {
      /// handle exceptional case?
      weight = 0;
      std::cout << "dividing by " << dual_area <<" dual_area, setting weight="<< weight << " for bone number" << tBnum << std::endl << std::flush;
    }
    //weight = 0;
    
    

    return weight;
  }



  /// obtain temporal bone weight (vertex in 2D)
  /// either by calculating it
  /// or retrieving from saved list 
  /// depending on timestamp
  double GettBoneWeight(const int & tBnum)
  {
    if (tbone_weights_timestamp[tBnum] < ma->GetTimeStamp() || tbone_weights_timestamp.find(tBnum) == tbone_weights_timestamp.end() )
    {
      double weight = CalctBoneWeight(tBnum);
      
      tbone_weights[tBnum] = weight;
      tbone_weights_timestamp[tBnum] = ma->GetTimeStamp();
    }
    else
    {
      //std::cout << "retrieving saved tbone weights" << std::endl << std::flush;
    }

    return tbone_weights[tBnum];

  };


  void RefreshtBoneWeights()
  {
    /// iterate trough all temporal bones
    /// vertices in 2D

    std::cout << "Starting tbone refresh" << std::endl << std::flush;
    std::cout << "Number of tbones (vertices): " << ma->GetNV() << std::endl << std::flush;
    //std::cout << ma->GetNV() << std::endl << std::flush;
    for(size_t tBnum = 0; tBnum < ma->GetNV(); tBnum++)
    {
        //std::cout <<  "tbone " << tBnum << std::endl << std::flush;
        double boneweight = GettBoneWeight(tBnum);
    }

  };


  double CalcBoneWeight(const int & Bnum)
  {
    /// calculate weight of the spatial bone (inverse length of barycentric dual edge)
    double weight = 0.;

    /// calculate center of edge
    Array<int> EVnums{size_t(0), size_t(0)};
    //ma->GetEdgePNums(Bnum, EVnums);
    INT<2> _EVnums = ma->GetEdgePNums(Bnum);
    EVnums[0] = _EVnums[0];
    EVnums[1] = _EVnums[1];
    //Array<int> EVnums ( ma->GetEdgePNums(Bnum) );
    
    Vec<2, double> Bc = 0.5*( ma->GetPoint<2>(size_t(EVnums[0])) + ma->GetPoint<2>(size_t(EVnums[1])) );
    Vec<2, double> vecl = ma->GetPoint<2>(size_t(EVnums[1])) - ma->GetPoint<2>(size_t(EVnums[0])) ;
    //Vec<2, double> Bc = EBarycenter(Enum);


    /// obtain elements hanging at edge
    Array<int> Elnums;
    ma->GetEdgeElements(Bnum, Elnums);

    /// sum up the distance from edge barycenter to element barycenters
    double dual_length = 0;
    for(int Elnum : Elnums)
    {
      //Array<int> ElVnums;
      //ma->GetElVertices(Elnum, ElVnums);
      Array<int> ElVnums (ma->GetElVertices(ElementId(VOL, Elnum)));
      //Array<int> ElVnums = (ma->GetElVertices(ElementId(VOL, Elnum)));

      /// get barycenter of element 
      // Vec<2, double> elc{0., 0.};
      // for(int i = 0; i < ElVnums.Size() ; i++)
      // {
      //   elc = elc + ma->GetPoint<2>(size_t(ElVnums[i]));
      // }
      // elc = (1./double(ElVnums.Size()))*elc;

      Vec<2, double> elc = ElBarycenter(Elnum);

      Vec<2,double> l_elc_Bc = elc - Bc;

      dual_length += sqrt(l_elc_Bc[0]*l_elc_Bc[0] + l_elc_Bc[1]*l_elc_Bc[1]);
    }

    //std::cout << "Bnum: " << Bnum << std::endl << std::flush;
    //std::cout << "length: " << length << std::endl << std::flush;

    /// The facets at the border have distorted weights (smaller dual edge in the mesh)
    /// fix this for each border facet
    if(!FacetVorb(Bnum))
    {
      if(Elnums.Size() == 1)
      {
        dual_length = 2.*dual_length;
      }
    }

    double length = sqrt(vecl[0]*vecl[0] + vecl[1]*vecl[1]);


    if (dual_length*dual_length > 0.00000001)
    {
      //weight = dual_length;
      //weight = 1./dual_length;  
      //weight = 1./(dual_length*dual_length);  
      weight = length*length/(dt*dt*dual_length*dual_length); 
      /// TODO: 3D dual cell factor? = mass matrix factor?
      ///weight *= 1./3.; 
    }
    else
    {
      /// handle exceptional case?
      weight = 0;
      std::cout << "dividing by " << dual_length <<" area, setting weight="<< weight << " for temporal bone number" << Bnum << std::endl << std::flush;
    }
    //weight = 0;

    /// careful: space-time metric with negative signature in time!
    return weight;
    //return -1.*weight;

  }

  /// obtain the bone weight (edges in 2D)
  /// either by calculating it
  /// or retrieving from saved list 
  /// depending on timestamp
  double GetBoneWeight(const int & Bnum)
  {
    if (bone_weights_timestamp[Bnum] < ma->GetTimeStamp() || bone_weights_timestamp.find(Bnum) == bone_weights_timestamp.end() )
    {
      double weight = CalcBoneWeight(Bnum);

      bone_weights[Bnum] = weight;
      bone_weights_timestamp[Bnum] = ma->GetTimeStamp();
    }
    else
    {
      //std::cout << "retrieving saved bone weights" << std::endl << std::flush;
    }

    return bone_weights[Bnum];

  };

  void RefreshBoneWeights()
  {
    /// iterate trough all bones
    /// edges in 2D

    std::cout << "Starting bone refresh" << std::endl << std::flush;
    std::cout << "Number of bones (edges): " << ma->GetNEdges() << std::endl << std::flush;
    for(int Bnum = 0; Bnum < ma->GetNEdges(); Bnum++)
    {
        //std::cout <<  "bone " << Bnum << std::endl << std::flush;
        double boneweight = GetBoneWeight(Bnum);
    }

  };



  ///  Manually pick the edges connected to a vertex in an oriented fashion (based on mesh geometry)
  /// return vector contains : int (edge number), bool (orientation: CCW = true ,CW = false)
  std::vector<std::tuple<int, bool>> CalcVertexOrientedEnums(const size_t& Vnum)
  {
    
    struct compareE
    {
      bool operator()(const std::tuple<int, bool, double>& lhs, const std::tuple<int, bool, double>& rhs) const
      {
        return std::get<0>(lhs) < std::get<0>(rhs);
      }
    };

    /// unordered set of edge numbers coupled with respective orientation and angle with respect to x-axis
    std::set<std::tuple<int, bool, double>, compareE> E_or_angles;


    /// go through all Elements hanging at the Vertex
    //Array<int> VElnums; 
    //ma->GetVertexElements(Vnum, VElnums);
    auto VElnums( ma->GetVertexElements(Vnum) );
    //VElnums[0] = _VElnums[0];
    //VElnums[1] = _VElnums[1];

    for (int Elnum : VElnums)
    {
      /// go through all Edges for a fixed Element
      for (int Enum : ma->GetElEdges(Elnum))
      {
        INT<2> EVnums = ma->GetEdgePNums(Enum);

        /// CAUTION: EVnums[0] < EVnums[1] ?
        /// from observation: seems so ...
        /// just in case: flip both
        if (EVnums[0] > EVnums[1])
        {
          std::cout << "EVnums[0] > EVnums[1] :" << std::endl;
          std::cout << EVnums[0] << " > " << EVnums[1] << std::endl;
          std::cout << "switching..." << std::endl;
          auto num = EVnums[0];
          EVnums[0] = EVnums[1];
          EVnums[1] = num;
        }

        /// required for globally setting a consistent link direction?
        /// here done: oriented edges that go from lower to higher vertex number are 'left tilted' (CCW tilted when viewed from "above")
        //std::cout << "Edge " << Enum << " has vertices " << EVnums[0] << " " << EVnums[1] << std::endl;

        /// does this edge hang at the desired Vertex?
        if(Vnum == EVnums[0] || Vnum == EVnums[1])
        {
          /// yes it does - add it!
          //std::cout << "Appending Edge " << Enum << std::endl;

          /// edge orientation - does the edge start at the desired vertex?
          /// p0i has the meaning of an index first, later as a bool
          bool p0i =  false;
          if(Vnum == EVnums[0])
          {
            /// this should be interpreted as 0 (index)
            p0i = false;
          }
          else if(Vnum == EVnums[1])
          {
            /// this should be interpreted as 1 (index)
            p0i = true;
          }

          /// calculate the edge angle with respect to the x axis (for sorting)
          
          /// starting point of the edge
          Vec<2, double> p0 = ma->GetPoint<2>(size_t(EVnums[int(p0i)]));
          //std::cout << "p0: " << p0[0] << " " << p0[1] << std::endl;
          /// ending point of the edge
          Vec<2, double> p1 = ma->GetPoint<2>(size_t(EVnums[int(!p0i)]));
          //std::cout << "p1: " << p1[0] << " " << p1[1] << std::endl;

          /// calculate unit vector starting at desired vertex
          double invnormvec01 = 1./sqrt( (p1 - p0)[0]*(p1 - p0)[0] + (p1 - p0)[1]*(p1 - p0)[1] ); 
          Vec<2, double> vec01 = invnormvec01*(p1 - p0);

          //std::cout << "vec01: " << vec01[0] << " " << vec01[1] << std::endl;

          /// calculate angle with respect to x axis
          double angle = atan2(vec01[1], vec01[0]);
          //double angle = double(Enum);

          /// bool(!p0i) value reflects, whether the edge starts at desired vertex
          /// this at the same time reflects the direction of the gauge link (forward = true, backward = false)
          /// geometrically this can be understood via 'left tilting'
          /// consistent with orientation based on global normal vector?
          /// i.e. Hdiv GridFunction with vec.data[:]=1.?
          /// this has been tested in python with 2D meshes

          E_or_angles.insert(std::make_tuple(Enum, !p0i, angle));
        }


      } 
    }

    /// Once all Edges hanging at the Vertex have been gathered
    /// Sort them with respect to the angle
    
    std::vector<std::tuple<int, bool, double>> E_or_angles_sorted(E_or_angles.begin(), E_or_angles.end());

    //std::cout << "Before sorting: " << std::endl;
    for (auto _E_or_angle : E_or_angles)
    {
      //std::cout << " " << std::get<0>(_E_or_angle);
    }
    //std::cout << std::endl;

    /// the compare operatore that orders the tuple elements based on second entry
    struct
    {
      bool operator()(const std::tuple<int, bool, double>& lhs, const std::tuple<int, bool, double>& rhs) const
      {
        return std::get<2>(lhs) < std::get<2>(rhs);
      }
    }
    compareangle;

    std::sort(E_or_angles_sorted.begin(), E_or_angles_sorted.end(), compareangle);

    //std::cout << "After sorting: " << std::endl;
    //for (auto _E_or_angle : E_or_angles_sorted)
    //{
    //  std::cout << " " << std::get<0>(_E_or_angle);
    //}
    //std::cout << std::endl;

    /// remove angles - they are unnecessary for return value
    std::vector<std::tuple<int, bool>> E_ors_sorted;

    for (auto _E_or_angle : E_or_angles_sorted)
    {
      E_ors_sorted.push_back( std::make_tuple(std::get<0>(_E_or_angle), std::get<1>(_E_or_angle)) );
    }

    return E_ors_sorted;
  };



  /// obtain the edges connected to a vertex
  /// either by calculation
  /// or retrieving from saved list 
  /// depending on timestamp
  std::vector<std::tuple<int, bool>> GetOrientedEdgesOfVertex(const int & Vnum)
  {
    //if (int(oriented_edges_of_vertex_timestamp[Vnum]) >= ma->GetTimeStamp())
    if (oriented_edges_of_vertex_timestamp[Vnum] < ma->GetTimeStamp() || oriented_edges_of_vertex_timestamp.find(Vnum) == oriented_edges_of_vertex_timestamp.end() )
    {  
        //std::cout << "calculating sorted oriented edges" << std::endl;
    
        std::vector<std::tuple<int, bool>> E_ors = CalcVertexOrientedEnums(Vnum);

        oriented_edges_of_vertex[Vnum] = E_ors;
        //std::cout << "old timestamp: " << oriented_edges_of_vertex_timestamp[Vnum] << std::endl;
        oriented_edges_of_vertex_timestamp[Vnum] = ma->GetTimeStamp();
        //std::cout << "new timestamp: " << oriented_edges_of_vertex_timestamp[Vnum] << std::endl;
    }
    else
    {
        //std::cout << "retrieving sorted oriented edges" << std::endl;
    }    
    return oriented_edges_of_vertex[Vnum];

  };


  void RefreshOrientedEdgesOfVertex()
  {
    /// iterate trough all vertices
    /// edges in 2D

    std::cout << "Starting vertex edge map refresh" << std::endl << std::flush;
    std::cout << "Number of vertices: " << ma->GetNV() << std::endl << std::flush;
    for(int Vnum = 0; Vnum < ma->GetNV(); Vnum++)
    {
        //std::cout <<  "vertex " << Vnum << std::endl << std::flush;
        GetOrientedEdgesOfVertex(Vnum);
    }

  };


  /// determine the relative orientation of gauge link at an edge
  /// returns true if link points out of the element
  /// returns false if link points into the element
  bool ElLinkOrientation(const int& Inum, const int& Elnum)
  {
    std::shared_ptr<MeshAccess> ma = this->GetMeshAccess();

    /// get vertices at the end of edges
    /// make sure that EVnums[0] < EVnums[1]
    /// which seems to be the default ...
    INT<2> EVnums = ma->GetEdgePNums(Inum);
    if (EVnums[0] > EVnums[1])
      {
        std::cout << "EVnums[0] > EVnums[1] :" << std::endl;
        std::cout << EVnums[0] << " > " << EVnums[1] << std::endl;
        std::cout << "switching..." << std::endl;
        auto num = EVnums[0];
        EVnums[0] = EVnums[1];
        EVnums[1] = num;
      }

    /// determine barycenter of element
    //Vec<2, double> ElC{0., 0.};
    Vec<2, double> ElC = ElBarycenter(Elnum);
    //std::cout << "Element " << Elnum << " has barycenter " << ElC << std::endl;

    /// left tilt the edge vector reflects global orientation
    Vec<2, double>EV0 = ma->GetPoint<2>(size_t(EVnums[0]));
    Vec<2, double>EV1 = ma->GetPoint<2>(size_t(EVnums[1]));

    // std::cout << "Edge " << Inum << " has points " << std::endl;
    // std::cout << "EV0 " << EV0 << std::endl;
    // std::cout << "EV1 " << EV1 << std::endl;

    /// this is the direction of the link
    Vec<2, double> tiltedgevec{ -1.*(EV1 - EV0)[1], (EV1 - EV0)[0] };
    Vec<2, double> dirvec = ( 0.5*(EV0 + EV1) ) - ElC;

    // std::cout << "tiltedgevec " << tiltedgevec << std::endl;
    // std::cout << "dirvec " << dirvec << std::endl;

    /// determine whether the element center is below or above the edge (plane)
    /// by projecting vectors (inner product)
    double s =  tiltedgevec[0]*dirvec[0] + tiltedgevec[1]*dirvec[1];
    //std::cout << "s " << s << std::endl;

    //std::cout << "Element " << Elnum << " Inum " << Inum << ": " << bool(s > 0) << std::endl;
    return bool(s > 0);
  };

  Vec<2, double> ElBarycenter(const int & Elnum) const
  {
    Vec<2, double> ElC{0., 0.};
    //Array<int> ElVnums;
    //ma->GetElVertices(Elnum, ElVnums);
    Array<int> ElVnums (ma->GetElVertices(ngcomp::ElementId(VOL, Elnum)));
    //Array<int> ElVnums = (ma->GetElVertices(ElementId(VOL, Elnum)));
    for(int i = 0; i < ElVnums.Size() ; i++)
    {
      ElC = ElC + ma->GetPoint<2>(size_t(ElVnums[i]));
    }
    ElC = (1./double(ElVnums.Size()))*ElC;

    return ElC;
  };

  Vec<2, double> EBarycenter(const int & Enum) const
  {
    Vec<2, double> EC{0., 0.};
    INT<2> EVnums = ma->GetEdgePNums(Enum);
    //Array<int> EVnums (ma->GetEdgePNums(ngcomp::ElementId(EDGE, Enum)));
    for(int i = 0; i < EVnums.Size() ; i++)
    {
      EC = EC + ma->GetPoint<2>(size_t(EVnums[i]));
    }
    EC = (1./double(EVnums.Size()))*EC;

    return EC;
  };

};

