# Decompose the energy
 state = simulation.context.getState(getEnergy=True,getForces=True,getParameters=True)
 print("The nonbonded energy is: "+str(calculate_nonbonded_energy(cgmodel).in_units_of(unit.kilojoules_per_mole)))
 potential_energy = state.getPotentialEnergy()
 print("The potential energy is: "+str(potential_energy.in_units_of(unit.kilojoules_per_mole)))
 kinetic_energy = state.getKineticEnergy()
 print("The kinetic energy is: "+str(kinetic_energy.in_units_of(unit.kilojoules_per_mole)))
 forces = state.getForces()
 energy_parameters = state.getParameters()
 print("The energy parameters are: "+str(energy_parameters))
#parameter_derivatives = state.getEnergyParameterDerivatives()
#print("The energy parameter derivatives are: "+str(parameter_derivatives))
# Figure out how many terms we have in the energy function

# Verify the non-bonded energy
 sum_nonbonded_energies = unit.Quantity(0.0,cgmodel.epsilon.unit)
 for interaction in nonbonded_interactions:
  dist = distance(positions[interaction[0]],positions[interaction[1]])
  nonbonded_energy = calculate_nonbonded_energy(cgmodel,particle1=interaction[0],particle2=interaction[1])
  sum_nonbonded_energies = sum_nonbonded_energies.__add__(nonbonded_energy)

