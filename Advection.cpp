struct Advection{
  static const int NQ{ 2 };

  using SolutionVector = std::array<ftype,NQ> ;
  SolutionVector u {};

  static constexpr struct {
    ftype a, b;
  //} model = {.a=1, .b=1};
  } model = {.a=0.5*0.5, .b=1};

  ftype& operator[] (int i){
    return u.at(i);
  }

  void fromInit(ftype x, ftype Lx){
    ftype step = (x>DEFINED_N/2)? 1:0;
    step = 2*step-1;
    //u =  {step, step/sqrt(model.a)};
    u =  {sin(2*M_PI*x/Lx), sin(2*M_PI*x/Lx)/sqrt(model.a)};
    //u =  {x, 2*x};
  }

  auto Flux(){
    Advection w {};
    w[0] = model.a*u[1];
    w[1] = model.b*u[0];
    return w;
  }

  auto Source(ftype t, ftype x, ftype dlt = 1){
    Advection w {};
    return w;
  }

  auto FluxMinus(){
    Advection w {};
    w[0] = 0 - 0.5*sqrt(model.a*model.b) * u[0] + 0.5*model.a * u[1];
    w[1] = 0 + 0.5*model.b * u[0] - 0.5*sqrt(model.a*model.b) * u[1];
    return w;
  }

  auto FluxPlus(){
    Advection w {u};
    w[0] = 0 + 0.5*sqrt(model.a*model.b) * u[0] + 0.5*model.a * u[1];
    w[1] = 0 + 0.5*model.b * u[0] + 0.5*sqrt(model.a*model.b) * u[1];
    return w;
  }
};
