#include <cstdio>
#include <math.h>
#include <armadillo>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/ranges.h>

using ftype=double;

#define MAX_MORDER 10
#define MAX_ITERS 10

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

#include "nonConservativeBurgers.cpp"
#include "flux.cpp"

template<int Mx, int Mt, typename T, int N>
struct DumbserMethod {

  ftype Lx {12};
  ftype dx {Lx/N};
  //ftype dt {dx/.4/10};
  ftype dt {dx/10};
  ftype Tmax {1};
  static constexpr int NBx {Mx+1};
  static constexpr int NBt {Mt+1};
  static constexpr int NQ {T::NQ};

  bool is_source_cell(int ix){ 
    //return true; 
    //if (ix==N/2)  return true; 
    return false;
  };

  arma::mat K1inv;
  arma::mat K2;
  arma::mat I4volflux; //matrix for volume flux

  using xDGdecomposition = std::array<T,NBx>; // U[ibx][iq] 
  using tDGdecomposition = std::array<T,NBt>; // U[ibt][iq] 
  using xtDGdecomposition = std::array<T,NBt*NBx>; // q[ibt*NBx+ibx][iq] 
  using Mesh = std::array<xDGdecomposition,N>; // cells[ix][ibx][iq]

  xtDGdecomposition left_cell_q; 
  xtDGdecomposition saved_q; 
  T left_cell_left_D; 
  T left_cell_left_flux; 

  Mesh cells {};

  DumbserMethod () :
    left_cell_q(), 
    saved_q(),
    left_cell_left_D() {
      // Init K1inv, K2, I
      ftype K1[NBt*NBx][NBt*NBx];
      ftype Kxi[NBt*NBx][NBt*NBx];
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          int ik {ikx*NBt + ikt}; 
          for (int ilx=0; ilx<NBx; ilx++) {
            for (int ilt=0; ilt<NBt; ilt++) {
              int il {ilx*NBt + ilt}; 
              ftype delta_kxlx = (ikx==ilx)? 1:0; 
              ftype delta_ktlt = (ikt==ilt)? 1:0; 
              K1[il][ik] = delta_kxlx * GAUSS_WEIGHTS[Mx][ikx] * ( 
                     POLYNOM_AT_ONE[Mt][ilt] * POLYNOM_AT_ONE[Mt][ikt]  
                   - GAUSS_WEIGHTS[Mt][ilt] * deriv_polynomial<Mt>(ikt, GAUSS_ROOTS[Mt][ilt]) 
                  );
              Kxi[il][ik] = GAUSS_WEIGHTS[Mx][ikx] * deriv_polynomial<Mx>(ilx, GAUSS_ROOTS[Mx][ikx]) * GAUSS_WEIGHTS[Mt][ilt] * delta_ktlt;
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

  void set_Lx_Courant(ftype _Lx, ftype courant){
    Lx = _Lx; 
    dx = Lx/N;
    dt = courant*dx;
    Tmax = T::Tmax;
    //fmt::print("#Lx = {}; Nx = {}; dx = {}, dt = {}\n", Lx, N, dx, dt);
  }

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

  xDGdecomposition TakeRootValue (int ix, ftype t = 0){
    xDGdecomposition U {}; 
    for (int ibx=0; ibx<NBx; ibx++) {
      ftype xikx = dx * (ix + GAUSS_ROOTS[Mx][ibx]);
      T udata; udata.fromInit(xikx, Lx, t);
      for (int iq=0; iq<NQ; iq++) {
        U[ibx][iq] = udata[iq];
      }
    }
    return U; 
  };

  void init(ftype t = 0) {
    for (int ix=0; ix<N; ix++){
      cells[ix] = TakeRootValue(ix, t);
    }
  };

  xtDGdecomposition ADER (xDGdecomposition u, int istep, int ix){
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
      xtDGdecomposition S {};
      xtDGdecomposition Q {};
      xtDGdecomposition B {};

      for (int ikx=0; ikx<NBx; ikx++) {  // for each
        for (int ikt=0; ikt<NBt; ikt++) { // for each
          int ik {ikx*NBt + ikt}; 
          auto Bm = q[ik].Bmatrix(); // eval the whole iq*ip matrix
          for (int ilx=0; ilx<NBx; ilx++) { // sum in l and in m
            auto Flxkt = q[ilx * NBt + ikt].Flux(); // eval the whole iq vector
            for (int iq=0; iq<NQ; iq++) { // for each
              F[ik][iq] += Flxkt[iq] * deriv_polynomial<Mx>(ilx, GAUSS_ROOTS[Mx][ikx]);  
              Q[ik][iq] += q[ilx * NBt + ikt][iq] * deriv_polynomial<Mx>(ilx, GAUSS_ROOTS[Mx][ikx]);  
            }
          }
          for (int iq=0; iq<NQ; iq++) {
            F[ik][iq] *= GAUSS_WEIGHTS[Mx][ikx] * GAUSS_WEIGHTS[Mt][ikt];  // simplify by dividing everything by w_ikx?
            for (int ip=0; ip<NQ; ip++) {
              B[ik][iq] += Bm(iq,ip) * Q[ik][ip];
            }
            B[ik][iq] *= GAUSS_WEIGHTS[Mx][ikx] * GAUSS_WEIGHTS[Mt][ikt];  
          }
          if (is_source_cell(ix)) {
            S[ik] = q[ik].Source( dt*(istep + GAUSS_ROOTS[Mt][ikt]), dx*(ix + GAUSS_ROOTS[Mx][ikx]) , ((ikt)==0)?1:0 );
            for (int iq=0; iq<NQ; iq++) {
              S[ik][iq] *=  GAUSS_WEIGHTS[Mx][ikx] * GAUSS_WEIGHTS[Mt][ikt];
              //fmt::print("{} {}     ||S||    {:20.18}\n",il,iq,S[il][iq]);
            }
          }
        }
      }

      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          int ik {ikx*NBt + ikt}; 
          for (int iq=0; iq<NQ; iq++) {
            q[ik][iq] = 0; 
            for (int il=0; il<NBt*NBx; il++) { //int il {ilx*NBt + ilt}; 
              q[ik][iq] += K1inv(ik,il) * (W[il][iq] - (dt/dx) * F[il][iq] -(dt/dx) * B[il][iq]) ;
              if (is_source_cell(ix)) {
                q[ik][iq] += K1inv(ik,il) * dt *S[il][iq];
              }
            }
          }
        }
      }
    }

    return q;
  }


  void ADER_update (int istep) {

    for (int ix=0; ix<N; ix++){
      xtDGdecomposition q {ADER(cells[(ix+N)%N], istep, ix)}; 

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

  void update (int istep) {
    for (int ix=0; ix<N+2; ix++){ // Two extra updates to prepare periodic boundary data; if(ix> 0 (or 1)) is below to account for it 

      xtDGdecomposition q {ADER(cells[(ix+N)%N], istep, ix)}; 
      if (ix==N+1) {
        q = saved_q;
      } else if (ix==1) {
        saved_q = q;
      }

      // dg-update. compute flux
      
      T left_cell_right_D {};
      T left_cell_right_flux {};

      if (ix > 0) {

        T D_integrated_dt {};
        T flux_integrated_dt {};
        for (int ikt=0; ikt<NBt; ikt++) {
          T qL {boundary_1_project_at_ti(left_cell_q, ikt)};
          T qR {boundary_0_project_at_ti(q, ikt)};
          auto nCD = nonConservativeD(qL,qR,dx,dt);
          auto SOF = SolomonOsherFlux(qL,qR,dx,dt);
          for (int iq=0; iq<NQ; iq++) {
            D_integrated_dt[iq] += GAUSS_WEIGHTS[Mt][ikt] * nCD[iq]; 
            flux_integrated_dt[iq] += GAUSS_WEIGHTS[Mt][ikt] * SOF[iq]; 
          }
        }
        for (int iq=0; iq<NQ; iq++) {
          left_cell_right_D[iq] = D_integrated_dt[iq];
          left_cell_right_flux[iq] = flux_integrated_dt[iq];
        }
      }

      // dg-update. update cell
      if (ix > 1) {
        for (int iq=0; iq<NQ; iq++) {
          for (int ikx=0; ikx<NBx; ikx++) {
            T volfluxp {};
            T addsourcep {};
            T smoothBq {};
            for (int ilx=0; ilx<NBx; ilx++) {
              for (int ilt=0; ilt<NBt; ilt++) {
                int il {ilx*NBt + ilt}; 
                auto Bm = q[ikx*NBt+ilt].Bmatrix(); // eval the whole iq*ip matrix
                T F = left_cell_q[il].Flux();
                for (int iq=0; iq<NQ; iq++){
                  volfluxp[iq] += F[iq] * GAUSS_WEIGHTS[Mt][ilt] * I4volflux(ilx,ikx);
                  for (int ip=0; ip<NQ; ip++){
                    smoothBq[iq] += Bm(iq,ip) * left_cell_q[il][ip] * GAUSS_WEIGHTS[Mt][ilt] * I4volflux(ikx,ilx);
                  }
                }
              }
            }
            if (is_source_cell((ix+N-1)%N)){ 
              for (int ilt=0; ilt<NBt; ilt++) {
                T S = left_cell_q[ikx*NBt + ilt].Source( dt*(istep + GAUSS_ROOTS[Mt][ilt]), dx*((ix+N-1)%N + GAUSS_ROOTS[Mx][ikx]) , (ikx==0)?1:0 );
                for (int iq=0; iq<NQ; iq++){
                  addsourcep[iq] += S[iq] * GAUSS_WEIGHTS[Mt][ilt] * GAUSS_WEIGHTS[Mx][ikx];
                }
              }
            }
            //cells[(N+ix-1)%N][ikx][iq] -= (1/(dx*GAUSS_WEIGHTS[Mx][ikx])) * ( 
                                                  //left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx] - 
                                                  //left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx] -
                                                  //dt * addfluxp[iq]
                                                  //- dt * dx * addsourcep[iq]
                                                  //);
            cells[(N+ix-1)%N][ikx][iq] += (1/GAUSS_WEIGHTS[Mx][ikx]) * ( 
                                                   + (dt/dx) * volfluxp[iq]
                                                   - (dt/dx) * smoothBq[iq]
                                                   - (dt/dx) * left_cell_right_D[iq] * POLYNOM_AT_ONE[Mx][ikx]
                                                   - (dt/dx) * left_cell_left_D[iq] * POLYNOM_AT_ZERO[Mx][ikx]
                                                   - (dt/dx) * left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx]
                                                   + (dt/dx) * left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx]
                                                  - dt * addsourcep[iq]
                );
            
            /*
            fmt::print(" + {} - {}\n - {} - {} - {} + {}\n", 
                (dt/dx) * volfluxp[iq],
                (dt/dx) * smoothBq[iq],
                (dt/dx) * left_cell_right_D[iq] * POLYNOM_AT_ONE[Mx][ikx],
                (dt/dx) * left_cell_left_D[iq] * POLYNOM_AT_ZERO[Mx][ikx],
                (dt/dx) * left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx],
                (dt/dx) * left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx]);
                */
          }
        }
      }

      left_cell_q = q;
      left_cell_left_D = left_cell_right_D;
      left_cell_left_flux = left_cell_right_flux;

    } //end ix loop 

  };

  void print_all (int istep) {
    std::string funcfilename = fmt::format("drop_{:09}.dat",istep);
    std::FILE* file = std::fopen( funcfilename.c_str(), "w");
    //fmt::print(file,"# {{ {:?}: {}, {:?}: {}, {:?}: {} }}\n", "dx", dx,"NQ",NQ, "NB", NB);
    fmt::print(file,"# {{");
    fmt::print(file,"  {:?}: {},", "dx", dx);
    fmt::print(file,"  {:?}: {},", "dt", dt);
    fmt::print(file,"  {:?}: {},", "NT", istep);
    fmt::print(file,"  {:?}: {},", "N", N);
    fmt::print(file,"  {:?}: {},", "NQ", NQ);
    fmt::print(file,"  {:?}: {},", "NBx", NBx);
    fmt::print(file,"  {:?}: {}", "NBt", NBt);
    fmt::print(file," }}\n");
    for (int ix=0; ix<N; ix++){
      for (int iq=0; iq<NQ; iq++) {
        for (int ibx=0; ibx<NBx; ibx++) {
          fmt::print(file,"{:30.18g}", cells[ix][ibx][iq]);
        }
      }
      fmt::print(file,"\n");
    }
    std::fclose(file);
  };
};

int main() {
  fmt::print(" {:>8} {:>8} {:>8} {:>6} {:>4}\n","dx", "dt", "T", "N", "NQ");
  {

      DumbserMethod<0, 0, nonConservativeBurgers::Burgers, 1000> mesh_calc;
      mesh_calc.init();
      int istep = 0;
  fmt::print(" {:>8} {:>8} \n",mesh_calc.dx, mesh_calc.dt);
      mesh_calc.print_all(istep);

      //for (; istep < floor(mesh_calc.Tmax/mesh_calc.dt+0.5); istep++) {
      for (; istep < 1000; istep++) {
     //   fmt::print("istep {}\n",istep);
        mesh_calc.update(istep);
        if (istep%100==0) mesh_calc.print_all(istep);
      }
      //mesh_calc.ADER_update(istep);
      //istep++;
      mesh_calc.print_all(istep);
      mesh_calc.init(istep*mesh_calc.dt);
      mesh_calc.print_all(-istep);
  }

   fmt::print("#Finished\n");
   return 0;
}



