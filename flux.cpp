template<typename T>
auto CIRFlux(T uL, T uR, ftype dx, ftype dt){
  T w {};
  for (int iq=0; iq< T::NQ; iq++) {
    w[iq] = uL.FluxPlus()[iq] + uR.FluxMinus()[iq];
  }
  return w;
}
template<typename T>
auto LaxFlux(T uL, T uR, ftype dx, ftype dt){
  T w {};
  for (int iq=0; iq< T::NQ; iq++) {
    w[iq] = 0.5*(uL.Flux()[iq] + uR.Flux()[iq]) - 0.5 * (dx/dt) *(uR[iq] - uL[iq]);
  }
  return w;
}
template<typename T>
auto RusanovFlux(T uL, T uR, ftype dx, ftype dt){
  T w {};
  auto fL = uL.Flux();
  auto fR = uR.Flux();
  T f{};
  for (int iq=0; iq< T::NQ; iq++){
    f [iq]= (fR[iq]>fL[iq])?fR[iq]:fL[iq];
  }
  for (int iq=0; iq< T::NQ; iq++) {
    w[iq] = 0.5*(uL.Flux()[iq] + uR.Flux()[iq]) - 0.5 * f[iq] *(uR[iq] - uL[iq]);
  }
  return w;
}
