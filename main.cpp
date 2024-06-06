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

#include "seismic.cpp"
#include "gas.cpp"
#include "Advection.cpp"
#include "Burgers.cpp"
#include "flux.cpp"
#include "eqs4testing.cpp"



template<int Mx, int Mt, typename T, int N>
struct DumbserMethod {

  ftype Lx {1.*N};
  ftype dx {1.};
  ftype dt {.01};
  ftype Tmax {1};
  static constexpr int NBx {Mx+1};
  static constexpr int NBt {Mt+1};
  static constexpr int NQ {T::NQ};

  bool is_source_cell(int ix){ 
    return true; 
    if (ix==N/2)  return true; 
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
    saved_q(),
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
      double checksum = 0;
      for (int ikx=0; ikx<NBx; ikx++) {
        checksum += GAUSS_WEIGHTS[Mx][ikx];
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
      for (int il=0; il<NBt*NBx; il++) {
        F[il] = q[il].Flux();
        if (is_source_cell(ix)) {
          S[il] = q[il].Source( dt*(istep + GAUSS_ROOTS[Mt][il%NBt]), dx*(ix + GAUSS_ROOTS[Mx][il/NBt]) );
          for (int iq=0; iq<NQ; iq++) {
            S[il][iq] *= dt * GAUSS_WEIGHTS[Mx][il/NBt] * GAUSS_WEIGHTS[Mt][il%NBt];
            //fmt::print("{} {}     ||S||    {:20.18}\n",il,iq,S[il][iq]);
          }
        }
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
                //fmt::print("{} {}     ||q||    {:20.18}\n",ik,iq,q[ik][iq]);
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
      
      T left_cell_right_flux {};

      if (ix > 0) {

        T Flux_integrated_dt {};
        for (int ikt=0; ikt<NBt; ikt++) {
          T qL {boundary_1_project_at_ti(left_cell_q, ikt)};
          T qR {boundary_0_project_at_ti(q, ikt)};
          for (int iq=0; iq<NQ; iq++) {
            //Flux_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * CIRFlux(qL,qR,dx,dt)[iq];  // for linear systems 
            Flux_integrated_dt[iq] += dt * GAUSS_WEIGHTS[Mt][ikt] * SolomonOsherFlux(qL,qR,dx,dt)[iq]; // for all systems
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
              }
            }
            if (is_source_cell((ix+N-1)%N)){ 
              for (int ilt=0; ilt<NBt; ilt++) {
                T S = left_cell_q[ikx*NBt + ilt].Source( dt*(istep + GAUSS_ROOTS[Mt][ilt]), dx*((ix+N-1)%N + GAUSS_ROOTS[Mx][ikx]) );
                for (int iq=0; iq<NQ; iq++){
                  addsourcep[iq] += S[iq] * GAUSS_WEIGHTS[Mt][ilt] * GAUSS_WEIGHTS[Mx][ikx];
                }
              }
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
    /*
    for (int ibx=0; ibx<NBx; ibx++) {
      ftype xikx = dx * (ix + GAUSS_ROOTS[Mx][ibx]);
      T udata; udata.fromInit(xikx, Lx, t);
      for (int iq=0; iq<NQ; iq++) {
        U[ibx][iq] = udata[iq];
      }
    }
    */
    
    ftype L2 {0};
    ftype Linf {0};
    for (int ix=0; ix<N; ix++){
      for (int ibx=0; ibx<NBx; ibx++) {
        ftype xikx = dx * (ix + GAUSS_ROOTS[Mx][ibx]) - 0 * istep * dt ;
        
        T udata; udata.fromInit(xikx, Lx, istep*dt);
        for (int iq=0; iq<NQ; iq++) {
          ftype t = udata[iq] - cells[ix][ibx][iq];
          L2 += t*t;
          ftype diff {fabs(udata[iq] - cells[ix][ibx][iq])};
          if (diff > Linf) { Linf = diff; }
        }
      }
    }
    L2 = sqrt(L2)/N;

    fmt::print(" {:8.5}",  dx);
    fmt::print(" {:8.5}",  dt);
    fmt::print(" {:8.5}",  istep*dt);
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



template<int MX, int MT, int N>
void one_full_calc(ftype courant){
    DumbserMethod<MX, MT, eqs4testing::Eqs4testing, N> mesh_calc;
    mesh_calc.set_Lx_Courant(1,courant);
    eqs4testing::model.set(mesh_calc.Lx);// >>>>>>>>>>>>>>> ? <<<<<<<<<<<<<<<<
    mesh_calc.init();
    int istep = 0;
    //int Nperiod = floor(mesh_calc.Lx/1/mesh_calc.dt + .5);
    int Nperiod = floor(mesh_calc.Tmax / mesh_calc.dt + .5);
    mesh_calc.print_all(istep);
    for (; istep < Nperiod; istep++) {
      //fmt::print("istep = {}\n",istep);
      mesh_calc.update(istep);
    }
    mesh_calc.print_error(istep);
}

template<int MX, int MT>
void several_full_calc(ftype courant){
  //one_full_calc<MX,MT,4>(courant);
  one_full_calc<MX,MT,8>(courant);
  one_full_calc<MX,MT,10>(courant);
  one_full_calc<MX,MT,16>(courant);
  one_full_calc<MX,MT,32>(courant);
}

int main() {
  fmt::print(" {:>8} {:>8} {:>8} {:>6} {:>4} {:>4} {:>4} {:>16} {:>16}\n","dx", "dt", "T", "N", "NQ", "NBx", "NBt", "L2", "Linf");
    //for (const auto courant_i : {1., .5, .1, .05, .01, 0.005, .001}) {
    for (const auto courant_i : {0.005, .001}) {
  //{ ftype courant_i = .001;
      several_full_calc<0,0>(courant_i);

      several_full_calc<1,1>(courant_i);
      several_full_calc<1,0>(courant_i);

      several_full_calc<2,0>(courant_i);
      several_full_calc<2,1>(courant_i);
      several_full_calc<2,2>(courant_i);

      several_full_calc<3,3>(courant_i);
      several_full_calc<4,4>(courant_i);
      several_full_calc<5,5>(courant_i);
      several_full_calc<6,6>(courant_i);
    }
  /*
  {

      DumbserMethod<6, 6, eqs4testing::Eqs4testing, 10> mesh_calc;
      mesh_calc.set_Lx_Courant(1,.001);

      eqs4testing::model.set(mesh_calc.Lx);// >>>>>>>>>>>>>>> ? <<<<<<<<<<<<<<<<

      mesh_calc.init();
      int istep = 0;
      mesh_calc.print_all(istep);
      
      for (; istep < .5/mesh_calc.dt; istep++) {
        mesh_calc.update(istep);
      }
      //mesh_calc.ADER_update(istep);
      //istep++;
      mesh_calc.print_all(istep);
      //mesh_calc.init(istep*mesh_calc.dt);
      //mesh_calc.print_all(-istep);
      mesh_calc.print_error(istep);
  }
  */

  /*
   gas::Gas u; 
   u.fromInit(0,1,0);
   auto K0 = u.Eigenvector(0);
   auto K1 = u.Eigenvector(1);
   auto f = SolomonOsherFlux(u,u,1,1);
   fmt::print("{} {}\n", u.u, K0, K1);
   */

   fmt::print("#Finished\n");
   return 0;
}



