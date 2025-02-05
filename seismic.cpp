namespace seismic {
  struct Material{
      ftype rho, cp;
      ftype drho; 
      ftype rhocp2;
      Material (ftype _rho=1., ftype _cp=1.): rho{_rho}, cp{_cp}, 
            drho {1./_rho},
            rhocp2 {_rho*_cp*_cp} {};
  };

  Material get_material(int im) {
    if (im==0) { 
      return Material(1,1);
    } else {
      return Material(1,1);
    }
  }

  ftype Ricker(ftype t, ftype a) {
    ftype delay = 5*a;
    ftype delayed_t = t-delay; 
    
    ftype R =  delayed_t * exp( -delayed_t*delayed_t/a/a)*sqrt(2)*exp(.5)/a; 
    //ftype R =  exp( -delayed_t*delayed_t/a/a) ; 
    fmt::print("{} {}\n", t, R);
    return R;
  }

  struct Seismic{
    static const int NQ{ 2 };
    static constexpr ftype Tmax{ .5 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};
    int material_index {};

    ftype& operator[] (int i){
      return u.at(i);
    }

    void fromInit(ftype x, ftype Lx, ftype t =0){

      Material m = get_material(material_index);
      ftype k = 2*M_PI/Lx;
      ftype omega = k*m.cp;
      ftype wave = sin(k*x - omega*t);

      u =  {wave, - wave/m.rho/m.cp};
    }

    auto Flux(){
      Seismic w {};
      Material m = get_material(material_index);
      w[0] = -m.rhocp2*u[1];
      w[1] = -m.drho*u[0];
      return w;
    }


    inline auto Eigenvector(int i){
      auto m = get_material(material_index);
      SolutionVector u{};
      if (i==0) {
        u[0] = m.rho; 
        u[1] = -1./m.cp;
      } else if (i==1) {
        u[0] = m.rho; 
        u[1] = 1./m.cp;
      } 
      return u;
    }

    auto Amatrix() {
      arma::mat Aflux(NQ, NQ);
      auto m = get_material(material_index);
      Aflux(0,0) = 0;
      Aflux(0,1) = -m.rhocp2;
      Aflux(1,0) = -m.drho; 
      Aflux(1,1) = 0;
      return Aflux;
    }

    auto Source(ftype t, ftype x, ftype dlt = 1){
      Seismic w {Eigenvector(0)};
      ftype s = dlt*Ricker(t,0.02);
      w[0] *= s;
      w[1] *= s;
      return w;
    }

    auto FluxMinus(){
      Seismic w {};
      Material m = get_material(material_index);
      w[0] = 0.5*( - m.rhocp2 * u[1] - m.cp * u[0] ) ;
      w[1] = 0.5*( - m.cp * u[1] - m.drho * u[0] ) ;
      return w;
    }

    auto FluxPlus(){
      Seismic w {};
      Material m = get_material(material_index);
      w[0] = 0.5*( - m.rhocp2 * u[1] + m.cp * u[0] ) ;
      w[1] = 0.5*( m.cp * u[1] - m.drho * u[0] ) ;
      return w;
    }
  };
}
