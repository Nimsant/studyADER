template<typename T>
auto CIRFlux(T uL, T uR, ftype dx, ftype dt){
  T w {};
  for (int iq=0; iq< T::NQ; iq++) {
    w[iq] = uL.FluxPlus()[iq] + uR.FluxMinus()[iq];
  }
  return w;
}

template<typename T>
auto Amod(T u) {
  auto Aflux = u.Amatrix();

  arma::cx_vec eigval;
  arma::cx_mat cx_eigvec;
  eig_gen(eigval,cx_eigvec, Aflux);
  arma::mat eigvec = arma::real(cx_eigvec);

  arma::mat Lmod (T::NQ, T::NQ, arma::fill::zeros);
  for (int iq = 0; iq<T::NQ; iq++) {
    ftype l = eigval(iq).real();
    Lmod(iq,iq) = (l>0)? l: -l;
  }
  arma::mat _Amod = eigvec * Lmod *eigvec.i();
  return _Amod;
}

template<typename T>
auto SolomonOsherFlux(T uL, T uR, ftype dx, ftype dt){
  
  auto fL = uL.Flux();
  auto fR = uR.Flux();

  arma::mat Aint_dQ (T::NQ, T::NQ);

  const int Oorder = 1;
  for (int ik=0; ik<Oorder+1; ik++){
    T wk {};
    ftype s = GAUSS_ROOTS[Oorder][ik];
    for (int iq=0; iq< T::NQ; iq++) {
      wk[iq] = uL[iq]*s + (1-s)*uR[iq];
    }
    arma::mat Amodk = Amod(wk);
    Aint_dQ += GAUSS_WEIGHTS[Oorder][ik] * Amodk;
  }

  T w {};

  for (int iq = 0; iq < T::NQ; iq++) {
    T AmodU {};
    for (int ip = 0; ip < T::NQ; ip++) {
      AmodU[iq] += Aint_dQ(iq,ip) * (uR[ip] - uL[ip]);
    }

    w[iq] = 0.5 * (fL[iq] + fR[iq]) - 0.5 * AmodU[iq];
  }
  
  return w;

}
