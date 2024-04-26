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

#include "seismic.cpp"
#include "Advection.cpp"
#include "Burgers.cpp"
#include "flux.cpp"
#include "eqs4testing.cpp"

template<int Mx, int Mt, typename T, int N>
struct DumbserMethod {

  ftype Lx {1.*N};
  ftype dx {1.};
  ftype dt {.1}; ///???????????????
  static constexpr int NBx {Mx+1};
  static constexpr int NBt {Mt+1};
  static constexpr int NQ {T::NQ};

  bool is_source_cell(int ix){ 
    return true ; 
    if (ix==N/2)  return true; 
    return false ; 
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
      T udata; udata.fromInit(xikx, Lx);
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
      for (int il=0; il<NBt*NBx; il++) {
        F[il] = q[il].Flux();
        if (is_source_cell(ix)) {
          S[il] = q[il].Source( dt*(istep + GAUSS_ROOTS[Mt][il%NBt]), dx*(ix + GAUSS_ROOTS[Mx][il/NBt]) );
          for (int iq=0; iq<NQ; iq++) {
            fmt::print("q[{}][{}] = {:16}\n", il, iq, q[il][iq]);
            fmt::print("S[{}][{}] = {:16}\n", il, iq, S[il][iq]);
            fmt::print("F[{}][{}] = {:16}\n", il, iq, F[il][iq]);
            fmt::print("(dt/dx)*F[{}][{}] = {:16}\n", il, iq, (dt/dx)*F[il][iq]);
            S[il][iq] *= dt * GAUSS_WEIGHTS[Mx][il/NBt] * GAUSS_WEIGHTS[Mt][il%NBt];
            fmt::print("dt*w*w*S[{}][{}] = {:16}\n", il, iq, S[il][iq]);
          }
        }
      }

      if (ix==0){
        K1inv.print("K1inv");
        I4volflux.print("I");
      }
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          int ik {ikx*NBt + ikt}; 
          for (int iq=0; iq<NQ; iq++) {
            q[ik][iq] = 0; 
            for (int il=0; il<NBt*NBx; il++) { //int il {ilx*NBt + ilt}; 
              q[ik][iq] += K1inv(ik,il)*W[il][iq] - K2(ik,il) * (dt/dx) * F[il][iq];
              if (is_source_cell(ix)) {
                q[ik][iq] += K1inv(ik,il) * S[il][iq];
              }
              if (ix==0){
                fmt::print("iter={} ik={} iq={} il={}| {} = {} ++++ {}*{}\n", iter, ik, iq, il, 
                    q[ik][iq],  
                    K2(ik,il) * (dt/dx)* F[il][iq],
                    S[il][iq], 
                    K1inv(il,ik) );
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
      if (true) {
        for (int ikx=0; ikx<NBx; ikx++) {
          for (int ikt=0; ikt<NBt; ikt++) {
            int ik {ikx*NBt + ikt}; 
            for (int iq=0; iq<NQ; iq++) {
              fmt::print("ikx={} ikt={} || q[{}][{}]={}\n",ik/NBt, ik%NBt, ik,iq,q[ik][iq]);
            }
          }
        }
      }

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
      } else if (ix==0) {
        saved_q = q;
      }

      // dg-update. compute flux
      
      T left_cell_right_flux {};

      if (ix > 0) {

        T Flux_integrated_dt {};
        for (int ikt=0; ikt<NBt; ikt++) {
          T qL {boundary_1_project_at_ti(left_cell_q, ikt)};
          T qR {boundary_0_project_at_ti(q, ikt)};
          for (int iq=0; iq<NQ; iq++) {
            //Flux_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * CIRFlux(qL,qR,dx,dt)[iq];  // for linear systems 
            //Flux_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * RusanovFlux(qL,qR,dx,dt)[iq]; // for all systems
            Flux_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * LaxFlux(qL,qR,dx,dt)[iq];  // doesn't look well
          }
            if (ix==10 || ix==9 || ix==11) {
              fmt::print(" {} {} {} Lax {} Lax\n", ix, qL.u, qR.u, LaxFlux(qL,qR,dx,dt).u);
            }
        }
        for (int iq=0; iq<NQ; iq++) {
          left_cell_right_flux[iq] = Flux_integrated_dt[iq];
        }
      }

      // dg-update. update cell
      if (ix > 1) {
        for (int iq=0; iq<NQ; iq++) {
          for (int ikx=0; ikx<NBx; ikx++) {
            T addfluxp {};
            T addsourcep {};
            for (int ilx=0; ilx<NBx; ilx++) {
              for (int ilt=0; ilt<NBt; ilt++) {
                int il {ilx*NBt + ilt}; 
                T F = left_cell_q[il].Flux();
                for (int iq=0; iq<NQ; iq++){
                  addfluxp[iq] += F[iq] * GAUSS_WEIGHTS[Mt][ilt] * I4volflux(ilx,ikx);
                }
                if (is_source_cell((ix+N-1)%N)){ 
                  T S = left_cell_q[ikx*NBt + ilt].Source( dt*(istep + GAUSS_ROOTS[Mt][ilt]), dx*((ix+N-1)%N + GAUSS_ROOTS[Mx][ikx]) );
                  for (int iq=0; iq<NQ; iq++){
                    addsourcep[iq] += S[iq] * GAUSS_WEIGHTS[Mt][ilt] * GAUSS_WEIGHTS[Mx][ikx];
                  }
                }
              }
            }
            if (ix==10) {
              fmt::print("Fr={:16.4} Fl={:16.4} vf={:16.4} s={:16.4} u={}\n",left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx], left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx], dt * addfluxp[iq],  dt * dx * addsourcep[iq], cells[(N+ix-1)%N][ikx][iq]);
            }
            cells[(N+ix-1)%N][ikx][iq] -= (1/(dx*GAUSS_WEIGHTS[Mx][ikx])) * ( 
                                                  left_cell_right_flux[iq] * POLYNOM_AT_ONE[Mx][ikx] - 
                                                  left_cell_left_flux[iq] * POLYNOM_AT_ZERO[Mx][ikx] -
                                                  dt * addfluxp[iq]
                                                  - dt * dx * addsourcep[iq]
                                                  );
          }
        }
      }

      left_cell_q = q;
      left_cell_left_flux = left_cell_right_flux;

    } //end ix loop 

  };


  void print_error (int istep) {
    //std::string funcfilename = fmt::format("error_{}.dat",istep);
    //std::FILE* file = std::fopen( funcfilename.c_str(), "w");
    //std::fclose(file);
    
    ftype L2 {0};
    ftype Linf {0};
    for (int ix=0; ix<N; ix++){
      for (int ibx=0; ibx<NBx; ibx++) {
        ftype xikx = dx * (ix + GAUSS_ROOTS[Mx][ibx]) - istep * dt ;
        T udata; udata.fromInit(xikx, Lx);
        for (int iq=0; iq<NQ; iq++) {
          ftype t = udata[iq] - cells[ix][ibx][iq];
          L2 += t*t;
          ftype diff {fabs(udata[iq] - cells[ix][ibx][iq])};
          if (diff > Linf) { Linf = diff; }
        }
      }
    }
    L2 = sqrt(L2)/N;

    fmt::print(" {:8}",  dx);
    fmt::print(" {:8}",  dt);
    fmt::print(" {:8}",  istep*dt);
    fmt::print(" {:6}",  N);
    fmt::print(" {:4}",  NQ);
    fmt::print(" {:4}",  NBx);
    fmt::print(" {:4}",  NBt);
    fmt::print(" {:16.9}",  L2);
    fmt::print(" {:16.9}",  Linf);
    fmt::print(" \n");
  }
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



template<int MX, int MT, int N>
void one_full_calc(){
    DumbserMethod<MX, MT, seismic::Seismic, N> mesh_calc;
    mesh_calc.init();
    ftype courant = mesh_calc.dt * seismic::get_material(0).cp / mesh_calc.dx;
    int istep = 0;
    for (; istep < N/courant; istep++) {
      mesh_calc.update();
    }
    mesh_calc.print_error(istep);
}

int main() {
  //fmt::print(" {:>8} {:>8} {:>8} {:>6} {:>4} {:>4} {:>4} {:>16} {:>16}\n","dx", "dt", "T", "N", "NQ", "NBx", "NBt", "L2", "Linf");
  {
    DumbserMethod<2, 2, Eqs4testing, 1> mesh_calc;
    mesh_calc.init();
    int istep = 0;
    mesh_calc.print_all(istep);
    for (; istep < 0; istep++) {
      mesh_calc.update(istep);
    }
    mesh_calc.ADER_update(istep);
    istep++;
    mesh_calc.print_all(istep);
    mesh_calc.print_error(istep);
  }


   fmt::print("#Finished\n");
   return 0;
}



