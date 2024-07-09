namespace eqs4testing {
  struct Model{
      ftype Au, Av, U0, V0, nu, k, omega;
    //} model = {.a=1, .b=1};
      void set(ftype Lx){
        Au = 0.1;
        Av = 0.3;
        U0 = 4;
        V0 = 6;
        nu = 10;
        k = 2*M_PI/Lx;
        omega = 2*M_PI/Lx;
      }
  }; 
  Model model;// = {.Au = 0.1, .Av = 0.3, .U0 = 4, .V0 = 6, .nu=0};
  struct Eqs4testing{
    static const int NQ{ 2 };

    static constexpr ftype Tmax{ .5 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};

    ftype& operator[] (int i){
      return u.at(i);
    }

    auto Solution(ftype t, ftype x){
      ftype phase = model.k * x - model.omega * t;
      Eqs4testing w {
        model.U0 + model.Au * sin(phase), 
        model.V0 + model.Av * cos(phase)
      };
      return w;
    }
    auto Solution_diff_t(ftype t, ftype x){
      ftype phase = model.k * x - model.omega * t;
      Eqs4testing w {
        - model.omega * model.Au * cos(phase), 
         model.omega * model.Av * sin(phase)
      };
      return w;
    }

    auto Solution_F_diff_x(ftype t, ftype x){
      ftype phase = model.k * x - model.omega * t;
      Eqs4testing w {
       -( model.Av*cos(phase) + model.V0 ) * model.Av * model.k * sin(phase),
        ( model.Au*sin(phase) + model.U0 ) * model.Au * model.k * cos(phase) 
      };
      return w;
    }

    void fromInit(ftype x, ftype Lx, ftype t = 0){
      //ftype step = (x>Lx*0.5)? 1: 2;
      //u = {step,step};
      u = Solution(t,x).u;
      //SolutionVector w { 1,1 };
      //u = w;
      //fmt::print("u {}\n", u);
    }

    auto Flux(){
      Eqs4testing w {};
      w[0] = 0.5*u[1]*u[1];
      w[1] = 0.5*u[0]*u[0];
      return w;
    }

    auto Amatrix() {
      arma::mat Aflux(NQ, NQ);
      Aflux(0,0) = 0;
      Aflux(0,1) = u[1];
      Aflux(1,0) = u[0]; 
      Aflux(1,1) = 0;
      return Aflux;
    }

    auto Source(ftype t, ftype x, ftype dlt=1){
      Eqs4testing w {};
      Eqs4testing ue {Solution(t,x)};
      Eqs4testing duedt {Solution_diff_t(t,x)} ;
      Eqs4testing dFedx {Solution_F_diff_x(t,x)} ;
      w[0] = - model.nu * ( u[0] - ue[0] ) + duedt[0] + dFedx[0];
      w[1] = - model.nu * ( u[1] - ue[1] ) + duedt[1] + dFedx[1];
      return w;
    }

  };
}
