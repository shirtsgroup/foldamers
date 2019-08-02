import numpy as np

kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)

def calculate_heat_capacity(E_expect,E2_expect,Temp_k):
    """
    """

    #------------------------------------------------------------------------
    # Compute Cv for NVT simulations as <E^2> - <E>^2 / (RT^2)
    #------------------------------------------------------------------------

    n=0

    if (n==0):
        print("")
        print("Computing Heat Capacity as ( <E^2> - <E>^2 ) / ( R*T^2 ) and as d<E>/dT")

    allCv_expect = np.zeros([K,ntypes,nBoots_work], np.float64)
    dCv_expect = np.zeros([K,ntypes],np.float64)

    allCv_expect[:,0,n] = (E2_expect - (E_expect*E_expect)) / ( kB * Temp_k**2)


    ####################################
    # C_v by fluctuation formula
    ####################################

    #Cv = (A - B^2) / (kT^2)
    # d2(Cv) = [1/(kT^2)]^2 [(dCv/dA)^2*d2A + 2*dCv*(dCv/dA)*(dCv/dB)*dAdB + (dCv/dB)^2*d2B]
    # = [1/(kT^2)]^2 [d2A - 4*B*dAdB + 4*B^2*d2B]
    # But this formula is not working for uncertainies!

    if (n==0):
        N_eff = (E2_expect - (E_expect*E_expect))/dE_expect**2  # sigma^2 / (sigma^2/n) = effective number of samples
        dCv_expect[:,0] = allCv_expect[:,0,n]*np.sqrt(2/N_eff)

    # only loop over the points that will be plotted, not the ones that

    ####################################
    # C_v by fluctuation formula
    ####################################

    #Cv = (A - B^2) / (kT^2)
    # d2(Cv) = [1/(kT^2)]^2 [(dCv/dA)^2*d2A + 2*dCv*(dCv/dA)*(dCv/dB)*dAdB + (dCv/dB)^2*d2B]
    # = [1/(kT^2)]^2 [d2A - 4*B*dAdB + 4*B^2*d2B]
    # But this formula is not working for uncertainies!

    if (n==0):
        N_eff = (E2_expect - (E_expect*E_expect))/dE_expect**2  # sigma^2 / (sigma^2/n) = effective number of samples
        dCv_expect[:,0] = allCv_expect[:,0,n]*np.sqrt(2/N_eff)

    # only loop over the points that will be plotted, not the ones that
    for i in range(originalK, K):

        # Now, calculae heat capacity by T-differences
        im = i-1
        ip = i+1
        if (i==originalK):
            im = originalK
        if (i==K-1):
            ip = i
        ####################################
        # C_v by first derivative of energy
        ####################################

        if (dertype == 'temperature'):  # temperature derivative
            # C_v = d<E>/dT
            allCv_expect[i,1,n] = (DeltaE_expect[im,ip])/(Temp_k[ip]-Temp_k[im])
            if (n==0):
                dCv_expect[i,1] = (dDeltaE_expect[im,ip])/(Temp_k[ip]-Temp_k[im])
        elif (dertype == 'beta'):  # beta derivative
            #Cv = d<E>/dT = dbeta/dT d<E>/beta = - kB*T(-2) d<E>/dbeta  = - kB beta^2 d<E>/dbeta
            allCv_expect[i,1,n] = kB * beta_k[i]**2 * (DeltaE_expect[ip,im])/(beta_k[ip]-beta_k[im])
            if (n==0):
                dCv_expect[i,1] = -kB * beta_k[i]**2 *(dDeltaE_expect[ip,im])/(beta_k[ip]-beta_k[im])

        ####################################
        # C_v by second derivative of free energy
        ####################################

        if (dertype == 'temperature'):
            # C_v = d<E>/dT = d/dT k_B T^2 df/dT = 2*T*df/dT + T^2*d^2f/dT^2

            if (i==originalK) or (i==K-1):
                 # We can't calculate this, set a number that will be printed as NAN
                 allCv_expect[i,2,n] = -10000000.0
            else:
                allCv_expect[i,2,n] = kB*Temp_k[i]*(2*df_ij[ip,im]/(Temp_k[ip]-Temp_k[im]) +
                                                    Temp_k[i]*(df_ij[ip,i]-df_ij[i,im])/
                                                    ((Temp_k[ip]-Temp_k[im])/(ip-im))**2)

            if (n==0):
                # Previous work to calculate the uncertainty commented out, should be cleaned up eventually
                # all_Cv_expect[i,2,n] = kB*Temp_k[i]*(2*df_ij[ip,i]+df_ij[i,im]/(Temp_k[ip]-Temp_k[im]) + Temp_k[i]*(df_ij[ip,i]-df_ij[i,im])/(Temp_k[ip]-Temp_k[i])**2)
                #all_Cv_expect[i,2,n] = kB*([2*Temp_k[i]/(Temp_k[ip]-Temp_k[im]) + Temp_k[i]**2/(Temp_k[ip]-Temp_k[i])**2]*df_ij[ip,i] + [2*Temp_k[i]/(Temp_k[ip]-Temp_k[im]) - Temp_k[i]**2/(Temp_k[ip]-Temp_k[i])**2]) df_ij[i,im]
                #all_Cv_expect[i,2,n] = kB*(A df_ij[ip,i] + B df_ij[i,im]
                A = 2*Temp_k[i]/(Temp_k[ip]-Temp_k[im]) + 4*Temp_k[i]**2/(Temp_k[ip]-Temp_k[im])**2
                B = 2*Temp_k[i]/(Temp_k[ip]-Temp_k[im]) + 4*Temp_k[i]**2/(Temp_k[ip]-Temp_k[im])**2
                #dCv_expect[i,2,n] = kB* [(A ddf_ij[ip,i])**2 + (B sdf_ij[i,im])**2 + 2*A*B*cov(df_ij[ip,i],df_ij[i,im])
                # This isn't it either: need to figure out that last term.
                dCv_expect[i,2] = kB*((A*ddf_ij[ip,i])**2 + (B*ddf_ij[i,im])**2)
                # Would need to add function computing covariance of DDG, (A-B)-(C-D)

        elif (dertype == 'beta'):
            # if beta is evenly spaced, rather than t, we can do 2nd derivative in beta
            # C_v = d<E>/dT = d/dT (df/dbeta) = dbeta/dT d/dbeta (df/dbeta) = -k_b beta^2 df^2/d^2beta
            if (i==originalK) or (i==K-1):
                #Flag as N/A -- we don't try to compute at the endpoints for now
                allCv_expect[i,2,n] = -10000000.0
            else:
                allCv_expect[i,2,n] = kB * beta_k[i]**2 *(df_ij[ip,i]-df_ij[i,im])/((beta_k[ip]-beta_k[im])/(ip-im))**2
            if (n==0):
                dCv_expect[i,2] = kB*(beta_k[i])**2 * (ddf_ij[ip,i]-ddf_ij[i,im])/((beta_k[ip]-beta_k[im])/(ip-im))**2
                # also wrong, need to be fixed.

    return()
