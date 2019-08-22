# VMD_python

Variational Mode Decomposition for Python

This is python realization for Variatioanl Mode Decomposition 

Authors: Konstantin Dragomiretskiy and Dominique Zosso

## Input and Parameters:

---------------------
signal  - the time domain signal (1D) to be decomposed

alpha   - the balancing parameter of the data-fidelity constraint

tau     - time-step of the dual ascent ( pick 0 for noise-slack )

K       - the number of modes to be recovered

DC      - true if the first mode is put and kept at DC (0-freq)

init    - 0 = all omegas start at 0
        - 1 = all omegas start uniformly distributed          
        - 2 = all omegas initialized randomly
                    
tol     - tolerance of convergence criterion; typically around 1e-6


## Output:
-------
u       - the collection of decomposed modes

u_hat   - spectra of the modes

omega   - estimated mode center-frequencies

## When using this code, please do cite the paper:
-----------------------------------------------

K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans. on Signal Processing (in press)
please check here for update reference: 
http://dx.doi.org/10.1109/TSP.2013.2288675
