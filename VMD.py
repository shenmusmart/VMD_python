def VMD(signal, alpha, tau, K, DC, init, tol):
    # ---------------------
    #  signal  - the time domain signal (1D) to be decomposed
    #  alpha   - the balancing parameter of the data-fidelity constraint
    #  tau     - time-step of the dual ascent ( pick 0 for noise-slack )
    #  K       - the number of modes to be recovered
    #  DC      - true if the first mode is put and kept at DC (0-freq)
    #  init    - 0 = all omegas start at 0
    #                     1 = all omegas start uniformly distributed
    #                     2 = all omegas initialized randomly
    #  tol     - tolerance of convergence criterion; typically around 1e-6
    #
    #  Output:
    #  -------
    #  u       - the collection of decomposed modes
    #  u_hat   - spectra of the modes
    #  omega   - estimated mode center-frequencies
    #

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    # Period and sampling frequency of input signal
    save_T=len(signal)
    fs=1/float(save_T)

    # extend the signal by mirroring
    T=save_T
    # print(T)
    f_mirror=np.zeros(2*T)
    print(f_mirror)
    f_mirror[0:T//2]=signal[T//2-1::-1]
    # print(f_mirror)
    f_mirror[T//2:3*T//2]= signal
    # print(f_mirror)
    f_mirror[3*T//2:2*T]=signal[-1:-T//2-1:-1]
    # print(f_mirror)
    f=f_mirror
    # print('f_mirror')
    # print(f_mirror)
    print('-------')

    # Time Domain 0 to T (of mirrored signal)
    T=float(len(f))
    # print(T)
    t=np.linspace(1/float(T),1,int(T),endpoint=True)
    # print(t)

    # Spectral Domain discretization
    freqs=t-0.5-1/T
    # print(freqs)
    # print('-----')
    # Maximum number of iterations (if not converged yet, then it won't anyway)
    N=500

    # For future generalizations: individual alpha for each mode
    Alpha=alpha*np.ones(K,dtype=complex)
    # print(Alpha.shape)
    # print(Alpha)
    # print('-----')

    # Construct and center f_hat
    f_hat=np.fft.fftshift(np.fft.fft(f))
    # print('f_hat')
    # print(f_hat.shape)
    # print(f_hat)
    # print('-----')
    f_hat_plus=f_hat
    f_hat_plus[0:int(int(T)/2)]=0
    # print('f_hat_plus')
    # print(f_hat_plus.shape)
    # print(f_hat_plus)
    # print('-----')
    # matrix keeping track of every iterant // could be discarded for mem
    u_hat_plus=np.zeros((N,len(freqs),K),dtype=complex)
    # print('u_hat_plus')
    # print(u_hat_plus.shape)
    # print(u_hat_plus)
    # print('-----')


    # Initialization of omega_k
    omega_plus=np.zeros((N,K),dtype=complex)
    # print('omega_plus')
    # print(omega_plus.shape)
    # print(omega_plus)
                        
    if (init==1):
        for i in range(1,K+1):
            omega_plus[0,i-1]=(0.5/K)*(i-1)
    elif (init==2):
        omega_plus[0,:]=np.sort(math.exp(math.log(fs))+(math.log(0.5)-math.log(fs))*np.random.rand(1,K))
    else:
        omega_plus[0,:]=0

    if (DC):
        omega_plus[0,0]=0

    # print('omega_plus')
    # print(omega_plus.shape)
    # print(omega_plus)

    # start with empty dual variables
    lamda_hat=np.zeros((N,len(freqs)),dtype=complex)

    # other inits
    uDiff=tol+2.2204e-16 #updata step
    # print('uDiff')
    # print(uDiff)
    # print('----')
    n=1 #loop counter
    sum_uk=0 #accumulator

    T=int(T)


    # ----------- Main loop for iterative updates

    while uDiff > tol and n<N:
        # update first mode accumulator
        k=1
        sum_uk = u_hat_plus[n-1,:,K-1]+sum_uk-u_hat_plus[n-1,:,0]
    #     print('sum_uk')
    #     print(sum_uk)
        #update spectrum of first mode through Wiener filter of residuals
        u_hat_plus[n,:,k-1]=(f_hat_plus-sum_uk-lamda_hat[n-1,:]/2)/(1+Alpha[k-1]*np.square(freqs-omega_plus[n-1,k-1]))
    #     print('u_hat_plus')
    #     print(u_hat_plus.shape)
    #     print(u_hat_plus[n,:,k-1])
    #     print('-----')
        
        

        #update first omega if not held at 0
        if DC==False:
            omega_plus[n,k-1]=np.dot(freqs[T//2:T],np.square(np.abs(u_hat_plus[n,T//2:T,k-1])).T)/np.sum(np.square(np.abs(u_hat_plus[n,T//2:T,k-1])))


        for k in range(2,K+1):

            #accumulator
            sum_uk=u_hat_plus[n,:,k-2]+sum_uk-u_hat_plus[n-1,:,k-1]
    #         print('sum_uk'+str(k))
    #         print(sum_uk)


            #mode spectrum
            u_hat_plus[n,:,k-1]=(f_hat_plus-sum_uk-lamda_hat[n-1,:]/2)/(1+Alpha[k-1]*np.square(freqs-omega_plus[n-1,k-1]))
    #         print('u_hat_plus'+str(k))
    #         print(u_hat_plus[n,:,k-1])
            
            #center frequencies
            omega_plus[n,k-1]=np.dot(freqs[T//2:T],np.square(np.abs(u_hat_plus[n,T//2:T,k-1])).T)/np.sum(np.square(np.abs(u_hat_plus[n,T//2:T:,k-1])))
    #         print('omega_plus'+str(k))
    #         print(omega_plus[n,k-1])
        #Dual ascent
    #     print(u_hat_plus.shape)
        lamda_hat[n,:]=lamda_hat[n-1,:]+tau*(np.sum(u_hat_plus[n,:,:],axis=1)-f_hat_plus)
    #     print('lamda_hat'+str(n))
    #     print(lamda_hat[n,:])

        #loop counter
        n=n+1

        #converged yet?
        uDiff=2.2204e-16

        for i in range(1,K+1):
            uDiff=uDiff+1/float(T)*np.dot(u_hat_plus[n-1,:,i-1]-u_hat_plus[n-2,:,i-1],(np.conj(u_hat_plus[n-1,:,i-1]-u_hat_plus[n-2,:,i-1])).conj().T)

            
        
        uDiff=np.abs(uDiff)
        # print('uDiff')
        # print(uDiff)
        
    # print('f_hat_plus')
    # print(f_hat_plus.shape)
    # print(f_hat_plus)
    # print('-----')   
    # print('u_hat_plus')
    # print(u_hat_plus.shape)
    # print(u_hat_plus)
    # print('-----')
    # print('sum_uk')
    # print(sum_uk)
    # print('-----')
        
    # ------ Postprocessing and cleanup

    # discard empty space if converged early

    N=np.minimum(N,n)
    omega = omega_plus[0:N,:]

    # Signal reconstruction
    u_hat = np.zeros((T,K),dtype=complex)
    u_hat[T//2:T,:]= np.squeeze(u_hat_plus[N-1,T//2:T,:])
    # print('u_hat')
    # print(u_hat.shape)
    # print(u_hat)
    u_hat[T//2:0:-1,:]=np.squeeze(np.conj(u_hat_plus[N-1,T//2:T,:]))
    u_hat[0,:]=np.conj(u_hat[-1,:])
    # print('u_hat')
    # print(u_hat)
    u=np.zeros((K,len(t)),dtype=complex)

    for k in range(1,K+1):
        u[k-1,:]= np.real(np.fft.ifft(np.fft.ifftshift(u_hat[:,k-1])))


    # remove mirror part 
    u=u[:,T//4:3*T//4]

    # print(u_hat.shape)
    #recompute spectrum
    u_hat = np.zeros((T//2,K),dtype=complex)

    for k in range(1,K+1):
        u_hat[:,k-1]=np.fft.fftshift(np.fft.fft(u[k-1,:])).conj().T
        
    # print('-----')
    # print('-----')
    # print(u)
    # print('-----')
    # print(u_hat)
    # print('-----')
    # print(omega)
    plt.plot(signal)
    plt.show()
    for i in range(1,K+1):
         plt.plot(u[i-1,:])
         plt.show()

    return (u,u_hat,omega)

            



















