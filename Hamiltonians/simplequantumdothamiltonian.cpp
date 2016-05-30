#include "simplequantumdothamiltonian.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <cmath>
#include "particle.h"

using namespace std;

SimpleQuantumDotHamiltonian::SimpleQuantumDotHamiltonian(System* system,bool jastrow, bool interaction) :
    Hamiltonian(system, jastrow) {
        m_interaction = interaction;

}

double SimpleQuantumDotHamiltonian::computeLocalEnergy(std::vector<Particle *> particles) {

    m_kineticEnergy = computeKineticEnergy();
    m_potentialEnergy = computePotentialEnergy();

    return m_kineticEnergy + m_potentialEnergy;
}

double SimpleQuantumDotHamiltonian::computePotentialEnergy() {

    std::vector<Particle *> particles = m_system->getParticles();
    const double omega = m_system->getWaveFunction()->getParameters()[0];

    double potentialEnergy = 0.0;

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
        for(int j = 0; j < m_system->getNumberOfDimensions(); j++) {
            double xj = particles[i]->getPosition()[j];
            potentialEnergy += xj*xj;
        }
    }

    double interacting = 0;

    if(m_interaction) {

        for (int i=0; i < m_system->getNumberOfParticles(); i++) {

            // interacting potential energy
            for (int j=i+1; j < m_system->getNumberOfParticles(); j++) {
                double rij2 = 0;
                for (int dim=0; dim < m_system->getNumberOfDimensions(); dim++) {
                    rij2 += pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
                }
                interacting += 1.0 / sqrt(rij2);
            }
        }

    }


    return 0.5*omega*omega*potentialEnergy + interacting;
}

double SimpleQuantumDotHamiltonian::getKineticEnergy(std::vector<Particle *> particles) {
    return m_kineticEnergy; //*m_system->getWaveFunction()->evaluate(particles); //Multiply with psi_t
}

double SimpleQuantumDotHamiltonian::getPotentialEnergy(std::vector<Particle *> particles) {
    return m_potentialEnergy; //*m_system->getWaveFunction()->evaluate(particles); //Multiply with psi_t
}
