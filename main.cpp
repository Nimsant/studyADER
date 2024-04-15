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


struct AdvectionRelaxation{
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

  void fromInit(ftype x){
    ftype step = (x>DEFINED_N/2)? 1:0;
    step = 2*step-1;
    //u =  {step, step/sqrt(model.a)};
    //u =  {sin(2*M_PI*x/DEFINED_N), sin(2*M_PI*x/DEFINED_N)/sqrt(model.a)};
    u =  {x, 2*x};
  }

  auto Flux(){
    AdvectionRelaxation w {};
    w[0] = model.a*u[1];
    w[1] = model.b*u[0];
    return w;
  }

  auto FluxMinus(){
    SolutionVector w {u};
    u[0] = 0 - 0.5*sqrt(model.a*model.b) * w[0] + 0.5*model.a * w[1];
    u[1] = 0 + 0.5*model.b * w[0] - 0.5*sqrt(model.a*model.b) * w[1];
    return *this;
  }

  auto FluxPlus(){
    SolutionVector w {u};
    u[0] = 0 + 0.5*sqrt(model.a*model.b) * w[0] + 0.5*model.a * w[1];
    u[1] = 0 + 0.5*model.b * w[0] + 0.5*sqrt(model.a*model.b) * w[1];
    return *this;
  }
};

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
  };

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
//        for (int iter=0; iter<MAX_ITERS; iter++) {
//        arma::mat Fp = flux(qp);
//        for (int iq=0; iq<NQ; iq++) {
//          qp.col(iq) = K1inv*Wp.col(iq) - K2 * (dt/dx) * Fp.col(iq);

    K1inv.print("K1inv");
    K2.print("K2");
    fmt::print("qp = \n");
    for (int il=0; il<NBt*NBx; il++) {
      for (int iq=0; iq<NQ; iq++) {
        fmt::print("{} ",q[il][iq]);
      }
      fmt::print("\n");
    }
    fmt::print("Wp = \n");
    for (int il=0; il<NBt*NBx; il++) {
      for (int iq=0; iq<NQ; iq++) {
        fmt::print("{} ",W[il][iq]);
      }
      fmt::print("\n");
    }
    for (int iter=0; iter<MAX_ITERS; iter++) {
      xtDGdecomposition F {};
      for (int il=0; il<NBt*NBx; il++) {
        F[il] = q[il].Flux();
      }
      fmt::print("Fp = \n");
      for (int il=0; il<NBt*NBx; il++) {
        for (int iq=0; iq<NQ; iq++) {
          fmt::print("{} ",F[il][iq]);
        }
        fmt::print("\n");
      }
      for (int ikx=0; ikx<NBx; ikx++) {
        for (int ikt=0; ikt<NBt; ikt++) {
          int ik {ikx*NBt + ikt}; 
          for (int iq=0; iq<NQ; iq++) {
            q[ik][iq] = 0; 
            for (int il=0; il<NBt*NBx; il++) {
              q[ik][iq] += K1inv(ik,il)*W[il][iq];
              q[ik][iq] += - K2(ik,il)* (dt/dx)* F[il][iq];
            }
          }
        }
      }
    }
    fmt::print("qp = \n");
    for (int il=0; il<NBt*NBx; il++) {
      for (int iq=0; iq<NQ; iq++) {
        fmt::print("{} ",q[il][iq]);
      }
      fmt::print("\n");
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



//template<int Mx, int Mt, typename T, int N, int ddx>
int main() {
   DumbserMethod<DEFINED_MX, DEFINED_MT, AdvectionRelaxation, DEFINED_N, 1> mesh_calc;
   mesh_calc.init();
   int istep = 0;
   mesh_calc.print_all(istep);
   mesh_calc.ADER_update();
   istep++;
   mesh_calc.print_all(istep);
   fmt::print("Finished\n");
   return 0;
}



