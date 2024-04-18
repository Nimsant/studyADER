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

  auto FluxMinus(){
    Burgers w {};
    w[0] = 0;
    return w;
  }

  auto FluxPlus(){
    Burgers w {u};
    w[0] = u[0];
    return w;
  }

};
auto LaxFlux(Burgers uL, Burgers uR, ftype dx, ftype dt){
  Burgers w {};
  w[0] = 0.5*(uL.Flux()[0] + uR.Flux()[0]) - 0.5 * (dx/dt) *(uR[0] - uL[0]);
  return w;
}
auto RusanovFlux(Burgers uL, Burgers uR, ftype dx, ftype dt){
  Burgers w {};
  auto fL = uL.Flux();
  auto fR = uR.Flux();
  Burgers f{};
  for (int iq=0; iq<Burgers::NQ; iq++){
    f [iq]= (fR[iq]>fL[iq])?fR[iq]:fL[iq];
  }
  w[0] = 0.5*(uL.Flux()[0] + uR.Flux()[0]) - 0.5 * f[0] *(uR[0] - uL[0]);
  return w;
}
