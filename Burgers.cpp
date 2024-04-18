struct Burgers{
  static const int NQ{ 1 };

  using SolutionVector = std::array<ftype,NQ> ;
  SolutionVector u {};

  ftype& operator[] (int i){
    return u.at(i);
  }

  void fromInit(ftype x, ftype Lx){
    ftype step = (x>Lx*0.5)? 1: 0;
    u = {step};
  }

  auto Flux(){
    Burgers w {};
    w[0] = 0.5*u[0]*u[0];
    return w;
  }

};
