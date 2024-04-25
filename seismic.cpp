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

  ftype Ricker(ftype t, ftype delay, ftype a) {
    ftype delayed_t = t-delay; 
    return delayed_t * exp( -delayed_t*delayed_t/a/a); 
  }

  struct Seismic{
    static const int NQ{ 2 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};
    int material_index {};

    ftype& operator[] (int i){
      return u.at(i);
    }

    void fromInit(ftype x, ftype Lx){
      //u =  {sin(2*M_PI*x/DEFINED_N), sin(2*M_PI*x/DEFINED_N)/sqrt(model.a)};
      Material m = get_material(material_index);
      ftype wave = sin(2*M_PI*x/Lx);

      u =  {wave, - wave/m.rho/m.cp};
      u =  {0, 0};
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

    auto Source(ftype t, ftype x){
      Seismic w {Eigenvector(0)};
      ftype s = Ricker(t,5,1);
      //fmt::print("{:16} {:16}\n",t,s);
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
