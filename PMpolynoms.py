import numpy as np
from icecream import ic
import scipy.integrate as integrate

MAX_MORDER = 8

def polynomial_basis(M):
    coef = (0,)*max(0,(M+1)) +(1,)
    domain = [0,1]
    PMp1 = np.polynomial.legendre.Legendre(coef,domain=domain)

    roots = PMp1.roots()

    data = {
            'polynomials': [],
            'weights': [],
            'roots': roots,
            'derivatives': []
            }
    for l in range(M+1):
        dlt = np.array([1 if i==l else 0 for i in range(M+1)])
        M2findpsil = np.array([[p**power for power in range(M+1)] for p in roots])
        X = np.linalg.solve(M2findpsil, dlt)
        psil = np.polynomial.polynomial.Polynomial(X)
        data['polynomials'].append(psil)
        data['derivatives'].append(psil.deriv())
        data['weights'].append(integrate.quad(psil.__call__,0,1)[0])
    return data 

def MakePolynomialCPP(what='polynomials'):
  s = '''if constexpr (M==0) {'''
  if what=='polynomials':
    s+= "  return 1;\n"
  if what=='derivatives':
    s+= "  return 0;\n"
  for Morder in range(1,MAX_MORDER):
    s += f'}} else if (M=={Morder}) {{\n'
    plmns = polynomial_basis(Morder)
    for l, plmn_l in enumerate(plmns[what]):
      s += f'  if (l=={l}) {{ \n'
      s += f'    const ftype x1 = x;\n'
      max_power = len(plmn_l.coef)
      for power in range(2,max_power):
          s += f'    const ftype x{power} = x*x{power-1};\n'
      ss = ''
      for power in range(max_power):
          ss+= f' {plmn_l.coef[power]:+}*x{power} '.replace('*x0','')
      s += f'    return {ss}; }}\n'
  s += f'}}\n\n'
  with open(f'gen_{what}_l.cpp','w') as f:
    f.write(s)

def MakePolynomialRoots():
  s = f' const std::array<std::array<ftype,{MAX_MORDER+1}>,{MAX_MORDER}> GAUSS_ROOTS {{{{'
  for Morder in range(MAX_MORDER):
    pl = polynomial_basis(Morder)
    output = ', '.join(
            [f"{i:.17g}" for i in pl['roots']] +\
            ['0' for i in range(MAX_MORDER+1 - len(pl['roots']))]
            )
    s += f"    {{ {output} }},\n"
  s += '}};\n\n'
  with open('gen_roots.cpp','w') as f:
    f.write(s)

def MakePolynomialWeights():
  s = f' const std::array<std::array<ftype,{MAX_MORDER+1}>,{MAX_MORDER}> GAUSS_WEIGHTS {{{{'
  for Morder in range(MAX_MORDER):
    pl = polynomial_basis(Morder)
    output = ', '.join(
            [f"{i:.17g}" for i in pl['weights']] +\
            ['0' for i in range(MAX_MORDER+1 - len(pl['weights']))]
            )
    s += f"    {{ {output} }},\n"
  s += '}};\n\n'
  with open('gen_weights.cpp','w') as f:
    f.write(s)

def MakePolynomAt(where):
  s = f' const std::array<std::array<ftype,{MAX_MORDER+1}>,{MAX_MORDER}> POLYNOM_AT_{["ZERO","ONE"][where]} {{{{'
  for Morder in range(MAX_MORDER):
    pl = polynomial_basis(Morder)
    output = ', '.join(
            [f"{i(where):.17g}" for i in pl['polynomials']] +\
            ['0' for i in range(MAX_MORDER+1 - len(pl['weights']))]
            )
    s += f"    {{ {output} }},\n"
  s += '}};\n\n'
  with open('gen_polynom_at_'+['zero','one'][where] + '.cpp','w') as f:
    f.write(s)
  s += f'}}\n\n'

def MakeI():
  for Morder in range(MAX_MORDER):
    pl = polynomial_basis(Morder)
    print(f"M={Morder}")
    for l in range(Morder):
        print(pl['polynomials'][l])
        print(pl['derivatives'][l])
        print(pl['derivatives'][l])
    print(pl['weights'])
    print(pl['roots'])
      
   
if __name__ == "__main__":
    MakePolynomialCPP('polynomials')
    MakePolynomialCPP('derivatives')
    MakePolynomialRoots()
    MakePolynomialWeights()
    MakePolynomAt(0)
    MakePolynomAt(1)
    MakeI()

