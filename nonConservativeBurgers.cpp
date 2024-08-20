namespace nonConservativeBurgers {
  struct Model{
      ftype e1 {}, e2 {};
      ftype uminus {}, vminus {};
      ftype a {}, b {};
  }; 
  Model model = {.e1 = 1.e-4, .e2 = 1.e-4, 
                 .uminus=7.99, .vminus=11.01,
                 .a=1, .b=.1
  };

  struct Burgers{
    static const int NQ{ 2 };

    static constexpr ftype Tmax{ .1 };

    using SolutionVector = std::array<ftype,NQ> ;
    SolutionVector u {};

    ftype& operator[] (int i){
      return u.at(i);
    }

    auto Solution(ftype t, ftype x, ftype Lx){
      ftype sigma = 10;
      ftype k = 2*M_PI/Lx;
      Burgers w {
        x - model.a*t,
        x - model.b*t
      };
      return w;
    }

    void fromInit(ftype x, ftype Lx, ftype t = 0){
      u = Solution(t,x,Lx).u;
    }

    auto Flux(){
      Burgers w {};
      //Burgers w {
      //  model.a*(u[0]+0*u[1]),
      //  model.b*(0*u[0]+u[1]),
      //};
      return w;
    }

    auto Bmatrix() {
      arma::mat B(NQ, NQ, arma::fill::zeros);
      B(0,0) = model.a;//u[0];
      B(0,1) = 0*model.a;//u[0];
      B(1,0) = 0*model.b;//u[1]; 
      B(1,1) = model.b;//u[1];
      return B;
    }

    auto Amatrix() {
      //auto B = Bmatrix();
      arma::mat B(NQ, NQ, arma::fill::zeros);
      //B(0,0) = model.a;//u[0];
      //B(0,1) = 0*model.a;//u[0];
      //B(1,0) = 0*model.b;//u[1]; 
      //B(1,1) = model.b;//u[1];
      return B;
    }

    auto Source(ftype t, ftype x, ftype dlt=1){
      Burgers w {};
      return w;
    }

  };
}
