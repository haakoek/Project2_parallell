#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <iomanip>
#include <iostream>

using namespace std;

Hamiltonian::Hamiltonian(System* system, bool jastrow) {
    m_system = system;
    m_jastrow = jastrow;
}

double Hamiltonian::computeKineticEnergy() {

    double kineticEnergy = m_system->getWaveFunction()->computeLaplacian(0);

    if(m_jastrow) {
        return -0.5*(kineticEnergy+m_system->getWaveFunction()->compute_JastrowLaplacian());
    } else {
        return -0.5*kineticEnergy;
    }

}
