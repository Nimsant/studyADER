#include <cstdio>
#include <math.h>
#include <armadillo>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/ranges.h>

using ftype=double;

#define MAX_MORDER 10
#define MAX_ITERS 1

template <int M> ftype polynominal(int l, ftype x){
#include "gen_polynomials_l.cpp"
  return 1;
}
template <int M> ftype deriv_polynomial(int l, ftype x){
#include "gen_derivatives_l.cpp"
  return 1;
}

#include "gen_polynom_at_zero.cpp"
#include "gen_polynom_at_one.cpp"
#include "gen_weights.cpp"
#include "gen_roots.cpp"

template<int Mx, int Mt, typename T, int N, int ddx>
struct DumbserMethod {

  static constexpr ftype Lx {1.*N*ddx};
  static constexpr ftype dx {1./ddx};
  static constexpr ftype dt {.1}; ///???????????????
  static constexpr int NBx {Mx+1};
  static constexpr int NBt {Mt+1};
  static constexpr int NQ {T::NQ};

  arma::mat K1inv;
  arma::mat K2;
  arma::mat I4volflux; //matrix for volume flux

  using xDGdecomposition = std::array<T,NBx>; // U[ibx][iq] 
  using tDGdecomposition = std::array<T,NBt>; // U[ibt][iq] 
  using xtDGdecomposition = std::array<T,NBt*NBx>; // q[ibt*NBx+ibx][iq] 
  using Mesh = std::array<xDGdecomposition,N>; // cells[ix][ibx][iq]

  xtDGdecomposition left_cell_q; 
  T left_cell_left_flux; 

  Mesh cells {};

  DumbserMethod () :
    left_cell_q(), 
    left_cell_left_flux() {
      // Init K1inv, K2, I
      ftype K1[NBt*NBx][NBt*NBx];
      ftype Kxi[NBt*NBx][NBt*NBx];
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          //int ik {ikt*NBx + ikx}; 
          int ik {ikx*NBt + ikt}; 
          for (int ilx=0; ilx<NBx; ilx++) {
            for (int ilt=0; ilt<NBt; ilt++) {
              //int il {ilt*NBx + ilx}; 
              int il {ilx*NBt + ilt}; 
              ftype delta_kxlx = (ikx==ilx)? 1:0; 
              ftype delta_ktlt = (ikt==ilt)? 1:0; 
              K1[il][ik] = delta_kxlx * GAUSS_WEIGHTS[Mx][ikx] * ( 
                     POLYNOM_AT_ONE[Mt][ilt] * POLYNOM_AT_ONE[Mt][ikt]  
                   - GAUSS_WEIGHTS[Mt][ilt] * deriv_polynomial<Mt>(ikt, GAUSS_ROOTS[Mt][ilt]) 
                  );
              Kxi[il][ik] = GAUSS_WEIGHTS[Mx][ikx] * deriv_polynomial<Mx>(ilx, GAUSS_ROOTS[Mt][ikx]) * GAUSS_WEIGHTS[Mt][ilt] * delta_ktlt;
            }
          }
        }
      }
      arma::mat K1mat(&K1[0][0], NBt*NBx, NBt*NBx);
      arma::mat Kximat(&Kxi[0][0], NBt*NBx, NBt*NBx);
      K1inv = K1mat.i();
      K2 = K1mat.i()*Kximat;
      I4volflux.set_size(NBx,NBx);
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ilx=0; ilx<NBx; ilx++) {
          I4volflux(ilx,ikx) = GAUSS_WEIGHTS[Mx][ilx] * deriv_polynomial<Mx>(ikx, GAUSS_ROOTS[Mx][ilx]); 
        }
      }
  };

  // For Flux
  //
  T boundary_0_project_at_ti(xtDGdecomposition q, int ikt) {
    T qR; 
    for (int iq=0; iq<NQ; iq++){
      for (int ikx=0; ikx<NBx; ikx++) {
        int ik {ikx*NBt + ikt}; 
        qR[iq] += POLYNOM_AT_ZERO[Mx][ikx] * q[ik][iq];
      }
    }
    return qR;
  }
  T boundary_1_project_at_ti(xtDGdecomposition q, int ikt) {
    T qL; 
    for (int iq=0; iq<NQ; iq++){
      for (int ikx=0; ikx<NBx; ikx++) {
        int ik {ikx*NBt + ikt}; 
        qL[iq] += POLYNOM_AT_ONE[Mx][ikx] * q[ik][iq];
      }
    }
    return qL;
  }

  xDGdecomposition TakeRootValue (int ix, ftype dx){
    xDGdecomposition U {}; 
    for (int ibx=0; ibx<NBx; ibx++) {
      ftype xikx = dx * (ix + GAUSS_ROOTS[Mx][ibx]);
      T udata; udata.fromInit(xikx);
      for (int iq=0; iq<NQ; iq++) {
        U[ibx][iq] = udata[iq];
      }
    }
    return U; 
  };

  void init() {
    for (int ix; ix<N; ix++){
      cells[ix] = TakeRootValue(ix, dx);
    }
  };

  xtDGdecomposition ADER (xDGdecomposition u){
    xtDGdecomposition q {}; 
    xtDGdecomposition W {}; 
    for (int ikx=0; ikx<NBx; ikx++) {
      for (int ikt=0; ikt<NBt; ikt++) {
        int ik {ikx*NBt + ikt}; 
        for (int iq=0; iq<NQ; iq++) {
          q[ik][iq] = u[ikx][iq];
          W[ik][iq] = GAUSS_WEIGHTS[Mx][ikx] * u[ikx][iq] * POLYNOM_AT_ZERO[Mt][ikt];
        }
      }
    }
    for (int iter=0; iter<MAX_ITERS; iter++) {
      xtDGdecomposition F {};
      for (int il=0; il<NBt*NBx; il++) {
        F[il] = q[il].Flux();
      }
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          int ik {ikx*NBt + ikt}; 
          for (int iq=0; iq<NQ; iq++) {
            q[ik][iq] = 0; 
            for (int il=0; il<NBt*NBx; il++) {
              q[ik][iq] += K1inv(ik,il)*W[il][iq] - K2(ik,il)* (dt/dx)* F[il][iq];
            }
          }
        }
      }
    }

    return q;
  }


  void ADER_update () {

    for (int ix=0; ix<N; ix++){
      xtDGdecomposition q {ADER(cells[(ix+N)%N])}; 

      //Temporary: collect back to u
      for (int iq=0; iq<NQ; iq++) {
        for (int ikx=0; ikx<NBx; ikx++) {
          cells[ix][ikx][iq] = 0; 
          for (int ikt=0; ikt<NBt; ikt++) {
            int ik {ikx*NBt + ikt}; 
            cells[ix][ikx][iq] += q[ik][iq] * POLYNOM_AT_ONE[Mt][ikt]; 
          }
        }
      }


    } //end ix loop 

  };

  void update () {
    for (int ix=0; ix<N+2; ix++){ // Two extra updates to prepare periodic boundary data; if(ix> 0 (or 1)) is below to account for it 

      xtDGdecomposition q {ADER(cells[(ix+N)%N])}; 

      // dg-update. compute flux
      
      T left_cell_right_flux {};

      if (ix > 0) {
        T left_cell_Fx1_integrated_dt {};
        T Fx0_integrated_dt {};
        for (int ikt=0; ikt<NBt; ikt++) {
          T qL {boundary_1_project_at_ti(left_cell_q, ikt)};
          T qR {boundary_0_project_at_ti(q, ikt)};
          T left_cell_Fx1_at_ti {qL.FluxPlus()}; 
          T Fx0_at_ti {qR.FluxMinus()}; 
          for (int iq=0; iq<NQ; iq++) {
            Fx0_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * Fx0_at_ti[iq];
            left_cell_Fx1_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * left_cell_Fx1_at_ti[iq];
          }
        }
        for (int iq=0; iq<NQ; iq++) {
          left_cell_right_flux[iq] = - Fx0_integrated_dt[iq] + left_cell_Fx1_integrated_dt[iq];
        }
      }


      // dg-update. update cell
      if (ix > 1) {
        for (int iq=0; iq<NQ; iq++) {
          for (int ikx=0; ikx<NBx; ikx++) {
            T addfluxp {};
            for (int ilx=0; ilx<NBx; ilx++) {
              for (int ilt=0; ilt<NBt; ilt++) {
                int il {ilx*NBt + ilt}; 
                T F = left_cell_q[il].Flux();
                for (int iq=0; iq<NQ; iq++){
                  addfluxp[iq] += F[iq] * GAUSS_WEIGHTS[Mt][ilt] * I4volflux(ilx,ikx);
                }
              }
            }
            cells[(N+ix-1)%N][ikx][iq] -= (1/(dx*GAUSS_WEIGHTS[Mx][ikx])) * ( 
                                                  left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx] - 
                                                  left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx] -
                                                  dt * addfluxp[iq]
                                                  );
          }
        }
      }

      left_cell_q = q;
      left_cell_left_flux = left_cell_right_flux;
      

    } //end ix loop 

  };

  void print_all (int istep) {
    std::string funcfilename = fmt::format("drop_{}.dat",istep);
    std::FILE* file = std::fopen( funcfilename.c_str(), "w");
    //fmt::print(file,"# {{ {:?}: {}, {:?}: {}, {:?}: {} }}\n", "dx", dx,"NQ",NQ, "NB", NB);
    fmt::print(file,"# {{");
    fmt::print(file,"  {:?}: {},", "dx", dx);
    fmt::print(file,"  {:?}: {},", "dt", dt);
    fmt::print(file,"  {:?}: {},", "T", istep*dt);
    fmt::print(file,"  {:?}: {},", "N", N);
    fmt::print(file,"  {:?}: {},", "NQ", NQ);
    fmt::print(file,"  {:?}: {},", "NBx", NBx);
    fmt::print(file,"  {:?}: {}", "NBt", NBt);
    fmt::print(file," }}\n");
    for (int ix=0; ix<N; ix++){
      for (int iq=0; iq<NQ; iq++) {
        for (int ibx=0; ibx<NBx; ibx++) {
          fmt::print(file,"{:20.10g}", cells[ix][ibx][iq]);
        }
      }
      fmt::print(file,"\n");
    }
    std::fclose(file);
  };
};


#include "seismic.cpp"

int main() {
   DumbserMethod<DEFINED_MX, DEFINED_MT, seismic::Seismic, DEFINED_N, 1> mesh_calc;
   mesh_calc.init();
   int istep = 0;
   mesh_calc.print_all(istep);
   for (; istep < DEFINED_NT; istep++) {
     mesh_calc.update();
   }
   //mesh_calc.ADER_update();
   //istep++;
   mesh_calc.print_all(istep);
   fmt::print("Finished\n");
   return 0;
}



