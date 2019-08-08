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

def get_intermediate_temperatures(T_from_file,NumIntermediates,dertype):
        """
        """
        #------------------------------------------------------------------------
        # Insert Intermediate T's and corresponding blank U's and E's
        #------------------------------------------------------------------------
        kB = unit.Quantity(0.008314462,unit.kilojoule_per_mole)  #Boltzmann constant (Gas constant) in kJ/(mol*K)
        minT = T_from_file[0]
        maxT = T_from_file[len(T_from_file) - 1]
        #beta = 1/(k*BT)
        #T = 1/(kB*beta)
        if dertype == 'temperature':
           minv = minT
           maxv = maxT
        elif dertype == 'beta':   # actually going in the opposite direction as beta for logistical reasons
           minv = 1/(kB._value*minT)
           maxv = 1/(kB._value*maxT)
        deltas = []
        for i in range(1,len(T_from_file)):
         deltas.append((T_from_file[i]._value-T_from_file[i-1]._value)/(NumIntermediates+1))
         deltas.append((T_from_file[i]._value-T_from_file[i-1]._value)/(NumIntermediates+1))
        originalK = len(T_from_file)

        Temp_k = []
        val_k = []
        currentv = minv._value
#        print(deltas)
        if dertype == 'temperature':
           # Loop, inserting equally spaced T's at which we are interested in the properties
           delta_index = 0
           while (currentv <= maxv._value):
              print(delta_index)
              print(currentv)
              delta = deltas[delta_index]
              val_k = np.append(val_k, currentv)
              currentv = currentv + delta
              Temp_k = np.concatenate((Temp_k,np.array(val_k)))
              delta_index = delta_index + 1
        elif dertype == 'beta':
        # Loop, inserting equally spaced T's at which we are interested in the properties
           while (currentv >= maxv):
              val_k = np.append(val_k, currentv)
              currentv = currentv + delta
           Temp_k = np.concatenate((Temp_k,(1/(kB._value*np.array(val_k)))))


        return(Temp_k)

def get_mbar_expectation(E_kln,temperature_list,NumIntermediates,dertype='temperature',output=None,mbar=None):
        """
        """

        if mbar == None:
         NumTemps = len(temperature_list) # Last TEMP # + 1 (start counting at 1)

         if (dertype == 'temperature'): # if the temperatures are equally spaced
          types = ['var','dT','ddT']
         elif (dertype == 'beta'): # if the inverse temperatures are equally spaced.
          types = ['var','dbeta','ddbeta']
         else:
          print('type of finite difference not recognized must be \'beta\' or \'temperature\'')
          quit()
         ntypes = len(types)

         kB = unit.Quantity(0.008314462,unit.kilojoule_per_mole)  #Boltzmann constant (Gas constant) in kJ/(mol*K)
         T_from_file = np.array([temperature._value for temperature in temperature_list])
         E_from_file = E_kln
         originalK = len(T_from_file)
         N_k = np.zeros(originalK,np.int32)

         #g = np.zeros(K,np.float64)
         for k in range(originalK):  # subsample the energies
#          g[k] = pymbar.timeseries.statisticalInefficiency(E_from_file[k][k])
#          indices = np.array(pymbar.timeseries.subsampleCorrelatedData(E_from_file[k][k],g=g[k])) # indices of uncorrelated samples
          try:
            N_k[k] = len(E_from_file[k][k]) # number of uncorrelated samples
          except:
            N_k[k] = len(E_from_file[k])
          #E_from_file[k,0:N_k[k]] = E_from_file[k,indices]

         Temp_k = get_intermediate_temperatures(temperature_list,NumIntermediates,dertype)
         print(Temp_k)

         # Update number of states
         K = len(Temp_k)
         # Loop, inserting E's into blank matrix (leaving blanks only where new Ts are inserted)

         Nall_k = np.zeros([K], np.int32) # Number of samples (n) for each state (k) = number of iterations/energies

         try:
          E_kn = np.zeros([K,len(E_from_file[0][0])], np.float64)
          for k in range(originalK):
            #if k != 0: 
              E_kn[k+k*NumIntermediates,0:N_k[k]] = E_from_file[k,k,0:N_k[k]]
                
            #else:
              #E_kn[k,0:N_k[k]] = E_from_file[k,k,0:N_k[k]]
              Nall_k[k+k*NumIntermediates] = N_k[k]
         except:
          E_kn = np.zeros([K,len(E_from_file[0])], np.float64)
          for k in range(originalK):
              E_kn[k+k*NumIntermediates,0:N_k[k]] = E_from_file[k,0:N_k[k]]
              Nall_k[k+k*NumIntermediates] = N_k[k]


         #------------------------------------------------------------------------
         # Compute inverse temperatures
         #------------------------------------------------------------------------
         beta_k = 1 / (kB._value * Temp_k)

         #------------------------------------------------------------------------
         # Compute reduced potential energies
         #------------------------------------------------------------------------

         allE_expect = np.zeros([K], np.float64)
         allE2_expect = np.zeros([K], np.float64)
         dE_expect = np.zeros([K],np.float64)
         u_kn = np.zeros([K,sum(Nall_k)], np.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
         #index = 0
         for k in range(K):
           index = 0
           for l in range(K):
              u_kn[k,index:index+Nall_k[l]] = beta_k[k] * E_kn[l,0:Nall_k[l]]
              index = index + Nall_k[l]

         #------------------------------------------------------------------------
         # Initialize MBAR
         #------------------------------------------------------------------------

         print("")
         print("Initializing MBAR:")
         print("--K = number of Temperatures with data = %d" % (originalK))
         print("--L = number of total Temperatures = %d" % (K))
         print("--N = number of Energies per Temperature = %d" % (np.max(Nall_k)))

         mbar = pymbar.MBAR(u_kn, Nall_k, verbose=False, relative_tolerance=1e-12,initial_f_k = None)

         E_kn = u_kn  # not a copy, we are going to write over it, but we don't need it any more.
         for k in range(K):
               E_kn[k,:]*=beta_k[k]**(-1)  # get the 'unreduced' potential -- we can't take differences of reduced potentials because the beta is different; math is much more confusing with derivatives of the reduced potentials.

        else:

          E_kn = E_kln          
          Temp_k = get_intermediate_temperatures(temperature_list,NumIntermediates,dertype) 

        if output != None:
               results = mbar.computeExpectations(E_kn,output='differences', state_dependent = True)
               results = {'mu': results[0], 'sigma': results[1]}
               result = results['mu']
               dresult = results['sigma']
        else:
               results = mbar.computeExpectations(E_kn, state_dependent = True)
               results = {'mu': results[0], 'sigma': results[1]}
               result = results['mu']
               dresult = results['sigma']

        return(mbar,E_kn,result,dresult,Temp_k)

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
