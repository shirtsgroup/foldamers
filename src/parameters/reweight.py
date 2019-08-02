#!/usr/bin/python

import numpy as np
# OpenMM utilities
import mdtraj as md
from simtk import unit
import pymbar

def get_free_energy_differences(mbar):
        """
        """

        results = mbar.getFreeEnergyDifferences()
        results = {'Delta_f': results[0], 'dDelta_f': results[1]}
        df_ij = results['Delta_f']
        ddf_ij = results['dDelta_f']
        return(df_ij,ddf_ij)

def get_mbar_expectation(E_kln,temperature_list,NumIntermediates,NumIterations=1000,nBoots=0,dertype='temperature',output=None,state_dependent=False):
        """
        """

        NumTemps = len(temperature_list) # Last TEMP # + 1 (start counting at 1)

        if (dertype == 'temperature'): # if the temperatures are equally spaced
          types = ['var','dT','ddT']
        elif (dertype == 'beta'): # if the inverse temperatures are equally spaced.
          types = ['var','dbeta','ddbeta']
        else:
          print('type of finite difference not recognized must be \'beta\' or \'temperature\'')
          quit()
        ntypes = len(types)

        kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)
        T_from_file = np.array([temperature._value for temperature in temperature_list])
        E_from_file = E_kln
        K = len(T_from_file)
        N_k = np.zeros(K,np.int32)
        g = np.zeros(K,np.float64)

        #for k in range(K):  # subsample the energies
#          g[k] = pymbar.timeseries.statisticalInefficiency(E_from_file[k][k])
#          indices = np.array(pymbar.timeseries.subsampleCorrelatedData(E_from_file[k][k],g=g[k])) # indices of uncorrelated samples
          #N_k[k] = len(indices) # number of uncorrelated samples
          #E_from_file[k,0:N_k[k]] = E_from_file[k,indices]

        #------------------------------------------------------------------------
        # Insert Intermediate T's and corresponding blank U's and E's
        #------------------------------------------------------------------------
        Temp_k = T_from_file
        minT = T_from_file[0]
        maxT = T_from_file[len(T_from_file) - 1]
        #beta = 1/(k*BT)
        #T = 1/(kB*beta)
        if dertype == 'temperature':
          minv = minT
          maxv = maxT
        elif dertype == 'beta':   # actually going in the opposite direction as beta for logistical reasons
          minv = 1/(kB*minT)
          maxv = 1/(kB*maxT)
        delta = (maxv-minv)/(NumIntermediates-1)
        originalK = len(Temp_k)

        val_k = []
        currentv = minv
        if dertype == 'temperature':
           # Loop, inserting equally spaced T's at which we are interested in the properties
           while (currentv <= maxv):
              val_k = np.append(val_k, currentv)
              currentv = currentv + delta
           Temp_k = np.concatenate((Temp_k,np.array(val_k)))
        elif dertype == 'beta':
        # Loop, inserting equally spaced T's at which we are interested in the properties
           while (currentv >= maxv):
              val_k = np.append(val_k, currentv)
              currentv = currentv + delta
           Temp_k = np.concatenate((Temp_k,(1/(kB*np.array(val_k)))))

        # Update number of states
        K = len(Temp_k)
        # Loop, inserting E's into blank matrix (leaving blanks only where new Ts are inserted)

        Nall_k = np.zeros([K], np.int32) # Number of samples (n) for each state (k) = number of iterations/energies
        E_kn = np.zeros([K, NumIterations], np.float64)

        for k in range(originalK):
            E_kn[k,0:N_k[k]] = E_from_file[k,k,0:N_k[k]]
            Nall_k[k] = N_k[k]

        #------------------------------------------------------------------------
        # Compute inverse temperatures
        #------------------------------------------------------------------------
        beta_k = 1 / (kB * Temp_k)

        #------------------------------------------------------------------------
        # Compute reduced potential energies
        #------------------------------------------------------------------------

        u_kln = np.zeros([K,K,NumIterations], np.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
        E_kn_samp = np.zeros([K,NumIterations], np.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l

        nBoots_work = nBoots + 1 # we add +1 to the bootstrap number, as the zeroth bootstrap sample is the original

        allE_expect = np.zeros([K,nBoots_work], np.float64)
        allE2_expect = np.zeros([K,nBoots_work], np.float64)
        dE_expect = np.zeros([K],np.float64)


        for n in range(nBoots_work):
           if (n > 0):
              print("Bootstrap: %d/%d" % (n,nBoots))
           for k in range(K):
           # resample the results:
               if Nall_k[k] > 0:
                  if (n == 0):  # don't randomize the first one
                     booti = np.array(range(N_k[k]))
                  else:
                     booti=np.random.randint(Nall_k[k],size=Nall_k[k])
                  E_kn_samp[k,0:Nall_k[k]] = E_kn[k,booti]

           for k in range(K):
               for l in range(K):
                   u_kln[k,l,0:Nall_k[k]] = beta_k[l] * E_kn_samp[k,0:Nall_k[k]]

        #------------------------------------------------------------------------
        # Initialize MBAR
        #------------------------------------------------------------------------

           if (n==0):  # only print this information the first time
              print("")
              print("Initializing MBAR:")
              print("--K = number of Temperatures with data = %d" % (originalK))
              print("--L = number of total Temperatures = %d" % (K))
              print("--N = number of Energies per Temperature = %d" % (np.max(Nall_k)))

           if (n==0):
              initial_f_k = None # start from zero 
           else:
              initial_f_k = mbar.f_k # start from the previous final free energies to speed convergence

           mbar = pymbar.MBAR(u_kln, Nall_k, verbose=False, relative_tolerance=1e-12)

           E_kln = u_kln  # not a copy, we are going to write over it, but we don't need it any more.
           for k in range(K):
               E_kln[:,k,:]*=beta_k[k]**(-1)  # get the 'unreduced' potential -- we can't take differences of reduced potentials because the beta is different; math is much more confusing with derivatives of the reduced potentials.
           if output != None:
               results = mbar.computeExpectations(E_kln,output='differences', state_dependent = True)
               results = {'mu': results[0], 'sigma': results[1]}
               result = results['mu']
               dresult = results['sigma']
           else:
               results = mbar.computeExpectations(E_kln, state_dependent = True)
               results = {'mu': results[0], 'sigma': results[1]}
               result = results['mu']
               dresult = results['sigma']


        return(mbar,E_kln,result,dresult)

def get_expectation_value(data_set):
        """
        Evaluate the dimensionless 

        Parameters
        ----------

        cgmodel: CGModel() class object.

        Returns
        -------

        cgmodel: CGModel() class object.

        """

        return(cgmodel)
