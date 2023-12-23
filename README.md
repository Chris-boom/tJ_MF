# Slave Boson Mean Field Calculation For t-J Model

This project is a simple slave boson mean field theory of square lattice. It define the class tJSBMF. You need to specify the parametes of t-J model and self-consistent calculation parameters first, then use self_consistent() member function to calculate.

## parameters

We have define t=1. J is antiferromagnetic coupling, x is hole dopping, N*N is the lattice size, rtol and atol are self-consistent calculation converge condition, T is temperature, wave is a character to specify the symmetry('s' or 'd').

## output

Delta is the RVB parameter; B is the hopping parameter defined as $$B = \frac{1}{2}\sum_{\sigma}\langle f_{i,\sigma}^{\dagger}f_{i+x,\sigma} \rangle$$; h is the boson hopping parameter defined as $$h = \langle b_{i}^{\dagger}b_{i+x} \rangle$$; DeltaSC is the BCS superconductor parameter.

## others

The tJSBMF_Plot.py and main.cpp files are not important. You can rewrite it as you wish.
CMakeList.txt should be rewrite to fit your environment.
