namespace gas {

  struct Material{
      ftype a;
      Material (ftype _a=1.): a{_a} {};
  };
  Material get_material(int im) {
      return Material(1);
  }
  struct Gas{
    static const int NQ{ 2 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};
    int material_index {};

    ftype& operator[] (int i){
      return u.at(i);
    }

    void fromInit(ftype x, ftype Lx, ftype t =0){
      ftype rho = 1; 
      ftype small_u = 10.1;
      u =  {rho, rho*small_u};
    }

    auto Flux(){
      Gas w {};
      Material m = get_material(material_index);
      w[0] = u[1];
      w[1] = u[1] * u[1] / u[0] + m.a * m.a * u[0];
      return w;
    }


    inline auto Eigenvector(int i){
      auto m = get_material(material_index);
      SolutionVector K {};

      ftype small_u = u[1]/u[0];
      if (i==0) {
        K[0] = 1; 
        K[1] = small_u - m.a;
      } else if (i==1) {
        K[0] = 1; 
        K[1] = small_u + m.a;
      } 
      fmt::print("{} ||| {}  - {} ={}\n", u, small_u, m.a, K);
      return K;
    }


    auto Amatrix() {
      arma::mat Aflux(NQ, NQ);
      auto m = get_material(material_index);
      Aflux(0,0) = 0;
      Aflux(0,1) = 1;
      ftype s = u[1]/u[0];
      Aflux(1,0) = - s * s + m.a * m.a; 
      Aflux(1,1) = 2*s;
      return Aflux;
    }

    auto Source(ftype t, ftype x){
      Gas w {Eigenvector(0)};
      w = {0,0};
      return w;
    }


    auto FluxMinus(){
      Gas w {};
      return w;
    }

    auto FluxPlus(){
      Gas w {};
      return w;
    }
  };
}
