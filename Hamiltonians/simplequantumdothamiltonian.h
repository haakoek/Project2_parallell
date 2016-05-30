#pragma once
#include "hamiltonian.h"
#include <vector>

class SimpleQuantumDotHamiltonian : public Hamiltonian
{
public:
    SimpleQuantumDotHamiltonian(class System *system,bool jastrow, bool interaction);
    double computeLocalEnergy(std::vector<Particle *> particles);
    double getKineticEnergy(std::vector<class Particle*> particles);
    double getPotentialEnergy(std::vector<class Particle*> particles);
    double computePotentialEnergy();
private:
    bool m_interaction = true;
    double m_potentialEnergy = 0.0;
    double m_kineticEnergy = 0.0;
};
