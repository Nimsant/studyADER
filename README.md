Introduction
============

This is am educational project which serves to understand the ADER-DG method present in the many works of Michael Dumbser. 

The source paper with the mathematical formulation which is closest to the code in this project is 

 [1] arXiv:1304.5408v1 High Order Space-Time Adaptive WENO Finite Volume Schemes for Non-Conservative Hyperbolic Systems, Michael Dumbser, Arturo Hidalgo, Olindo Zanotti


The following papers are also used:

 [2] https://doi.org/10.1016/j.jcp.2008.05.025
 
 [3] https://doi.org/10.1016/j.jcp.2007.12.005

The present code is to be used to 

 - understand the ADER-DG mathematical expressions for the mesh updates
 - copy-paste expressions into more complex projects

The code features are limited:
 
 - the code is 1D
 - the elements are rectangular, and the mesh is uniform
 - no limiter is used
 - PnPm WENO reconstruction is not applied
 - several phyisical models are included, but only Seismic and eqs4testing are working in the latest version (the others may miss a function parameter or a source method, which you can add yousefl)

You need
========

 * C++ compiler (g++)
   * the Armadillo library is used for matrix inversion 
 * python3 
   * packages used: numpy, matplotlib, icecream

How to use
==========

1. Compile the LaTeX files with the theoretical background

 > make theory.pdf

2. Read the theory.pdf file which appeared. 

3. (Important!) Ask questions whenever any detail is not clear.

4. Go to main.cpp file and find the 'main' function. Fill it with the initializion of the mesh and data output. For example, 
      
```
  DumbserMethod<5, 5, seismic::Seismic, 16> mesh_calc;
  mesh_calc.init();
  mesh_calc.print_all(0);
```

The parameters are explained in theory.pdf. 

5. Compile with `make`, 
  
   or compile and run with `make run`, 

   or compile, run and plot the output with `make show`

6. In the plot, 

   6.1 The first component of the solution vector is plotted vs the x axis.
   
   6.2 If several outputs are present (i.e, the `print_all` method is called several times with different argument)

```
 mesh_calc.print_all(1);
 // ... update mesh
 mesh_calc.print_all(2);
```
   you can use interactive controls by pressing 'left' and 'right' arrow keys to switch between files. 

7. Try to understand what is being plotted and how to change this inital condition.

8. (Important!) Ask questions whenever anything is not working.

9. Vary the mesh template parameters, plot the output and see how the parameter choice affects the representation of the initial function on the mesh. 

10. Add one `ADER_update` method, or several `update` and `print_all` methods to see that the method works.

11. Try to include the `print_error` function to find the order of approximation of the method. Use the `several_full_calc` method for convenience

12. Look in the `ADER_update` and `update` expressions, understand them  and replicate them in your code (which is probably better than this one and is targeted to solve real problems)





