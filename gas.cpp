namespace gas {


  struct Gas{
    static const int NQ{ 3 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};

    static constexpr ftype Tmax{ .2 };
    static constexpr ftype gamma = 1.4;

    ftype& operator[] (int i){
      return u.at(i);
    }

    auto from_rhoup(ftype _rho, ftype _u, ftype _p){
      SolutionVector w {};
      w[0] = _rho; 
      w[1] = w[0]*_u;
      w[2] = w[0] * ( 0.5 * w[1] * w[1] / w[0] / w[0] + _p / w[0] / (gamma - 1) );
      return w;
    }

    void fromInit(ftype x, ftype Lx, ftype t =0){
      ftype rho=1, small_u = 0;
      //if (x < Lx/2) { u = from_rhoup(1., 0., 1.); } else { u = from_rhoup(0.1, 0., 1.); };   //RP1
      if (x < Lx/2) { u = from_rhoup(1., 0.75, 1.); } else { u = from_rhoup(0.125, 0., 0.1); };   //RP2
    }

    auto Flux(){ //Toro 3.1 page 89
      Gas w {};
      w[0] = u[1];
      w[1] = 0.5 * (3-gamma) * u[1] * u[1] / u[0] + (gamma - 1) * u[2];
      w[2] = gamma * u[1] * u[2] / u[0] - 0.5 * (gamma - 1) * u[1] * u[1] * u[1] / u[0] / u[0];
      return w;
    }


    inline auto Eigenvector(int i){
      SolutionVector K {};
      ftype p = (gamma - 1) * (u[2] - 0.5 * u[1] * u[1] / u[0]  );
      return K;
    }



    auto Amatrix() {
      arma::mat Aflux(NQ, NQ);

      Aflux(0,0) = 0.;
      Aflux(0,1) = 1.;
      Aflux(0,2) = 0.;

      ftype s = u[1] / u[0];

      Aflux(1,0) = -0.5 * (gamma - 3.) * s * s; 
      Aflux(1,1) = (3. - gamma) * s;
      Aflux(1,2) = gamma - 1.;

      ftype t = u[2] / u[0];

      Aflux(2,0) = - gamma * s * t + (gamma-1.) *s*s*s; 
      Aflux(2,1) = gamma * t - 1.5 * (gamma-1.) *s*s;
      Aflux(2,2) = gamma * s;

      return Aflux;
    }

    auto Source(ftype t, ftype x, int dlt){
      Gas w {};
      return w;
    }

  }; 
}
