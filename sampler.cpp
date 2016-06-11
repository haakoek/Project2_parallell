#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <iostream>
#include <fstream>
#include <mpi/mpi.h>
#include <string>
#include <cstdarg>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

Sampler::Sampler(System* system, bool writeToFile) {
    m_system = system;
    m_stepNumber = 0;
    m_writeToFile = writeToFile;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::setWriteToFile(bool writeToFile) {
    m_writeToFile = writeToFile;
}

void Sampler::sample(bool acceptedStep) {

    // Make sure the sampling variable(s) are initialized at the first step.

    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0.0;
        m_accepted = 0;
        m_E2 = 0.0;

        m_dEdalpha = 0.0;
        m_dE1alpha = 0.0;
        m_dE2alpha = 0.0;

        m_dEbeta   = 0.0;
        m_dE1beta  = 0.0;
        m_dE2beta  = 0.0;

        m_variance = 0.0;
        m_meanKineticEnergy = 0.0;
        m_meanPotentialEnergy = 0.0;
        m_meanDistance = 0.0;

        for(int i = 0; i < 600; i++) {
            m_probR[i] = 0.0;
        }

        char filename1[50];
        char filename2[50];
        int n1,n2;

        if(m_system->getJastrow()) {
            n1 = sprintf(filename1,"energies%d_N=%d_w=%g_WithJastrow.txt",m_system->getRank()+1,m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
            n2 = sprintf(filename2,"positions%d_N=%d_w=%g_WithJastrow.txt",m_system->getRank()+1,m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
        } else {
            n1 = sprintf(filename1,"energies%d_N=%d_w=%g.txt",m_system->getRank()+1,m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
            n2 = sprintf(filename2,"positions%d_N=%d_w=%g.txt",m_system->getRank()+1,m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
        }

        if(m_writeToFile) {
            m_outfile.open(filename1);
            m_positionsfile.open(filename2);
        }

    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */

    double localEnergy   = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());

    m_cumulativeEnergy  += localEnergy;
    m_E2                += localEnergy*localEnergy;
    m_stepNumber++;



    m_meanKineticEnergy   += m_system->getHamiltonian()->getKineticEnergy(m_system->getParticles());
    m_meanPotentialEnergy += m_system->getHamiltonian()->getPotentialEnergy(m_system->getParticles());

    //Samples for steepest descent

    m_system->getWaveFunction()->computeParametersDerivative();


    double dpsi_dalpha = m_system->getWaveFunction()->getAlphaDerivative();
    double dpsi_dbeta  = m_system->getWaveFunction()->getBetaDerivative();

    m_dE1alpha += dpsi_dalpha*localEnergy;
    m_dE2alpha += dpsi_dalpha;

    m_dE1beta  += dpsi_dbeta*localEnergy;
    m_dE2beta  += dpsi_dbeta;

    if(m_writeToFile) {

        m_outfile << localEnergy << endl;

        //Write positions to file in order to compute one-body density.
        for(int p = 0; p < m_system->getNumberOfParticles(); p++) {

            double r = 0.0;

            for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
                double xk = m_system->getParticles()[p]->getPosition()[d];
                r+= xk*xk;
            }

            r = sqrt(r);

            double dr = 0;
            double step = 0.1;

            for(int i = 0; i < 600; i++) {
                if( r >= dr && r < dr+step ) {
                    m_probR[i] += 1;
                    break;
                } else {
                    dr += step;
                }
            }
        }

    }



    if (m_stepNumber == (1-m_system->getEquilibrationFraction())*m_system->getNumberOfMetropolisSteps()) {

        //Close open files.
        cout << "Process: " << m_system->getRank()+1 << endl;

        for(int i = 0; i < 600; i++) {
            m_positionsfile << m_probR[i] << endl;
        }

        m_positionsfile.close();
        m_outfile.close();
    }


    if(acceptedStep == true) {
        //Update number of accepted steps.
        m_accepted++;
    }
}

void Sampler::printOutputToTerminal() {

    if(m_system->getRank() == 0) {

        char filename3[50];
        int n1;

        if(m_system->getJastrow()) {
            n1 = sprintf(filename3,"Averages_N=%d_w=%g_WithJastrow.txt",m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
        } else {
            n1 = sprintf(filename3,"Averages_N=%d_w=%g.txt",m_system->getNumberOfParticles(),m_system->getWaveFunction()->getParameters()[0]);
        }

        m_averages.open(filename3);

        int     np = m_system->getNumberOfParticles();
        int     nd = m_system->getNumberOfDimensions();
        int     ms = m_system->getNumberOfMetropolisSteps();
        int     p  = m_system->getWaveFunction()->getNumberOfParameters();
        double  ef = m_system->getEquilibrationFraction();
        std::vector<double> pa = m_system->getWaveFunction()->getParameters();

        cout << endl;
        cout << "  -- System info -- " << endl;
        cout << " Number of particles  : " << np << endl;
        cout << " Number of dimensions : " << nd << endl;
        cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
        cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
        cout << endl;
        cout << "  -- Wave function parameters -- " << endl;
        cout << " Number of parameters : " << p << endl;

        m_averages << "  -- System info -- " << endl;
        m_averages << " Number of particles  : " << np << endl;
        m_averages << " Number of dimensions : " << nd << endl;
        m_averages << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
        m_averages << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
        m_averages << endl;
        m_averages << "  -- Wave function parameters -- " << endl;
        m_averages << " Number of parameters : " << p << endl;


        for (int i=0; i < p; i++) {
            cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
            m_averages << " Parameter " << i+1 << " : " << pa.at(i) << endl;
        }
        cout << endl;

        m_averages << " Number of processors: " << m_system->getSize() << endl;

        cout << "  -- Results -- " << endl;
        m_averages << "  -- Results -- " << endl;

        cout << "Energy: " << m_energy << endl;
        m_averages <<"Energy: " << m_energy << endl;

        cout << "<E^2>: " << m_E2 << endl;
        m_averages << "<E^2>: " << m_E2 << endl;

        cout << "Variance: " << m_variance << endl;
        m_averages << "Variance: " << m_variance << endl;

        cout << "Std: " << sqrt(m_variance) << endl;
        m_averages << "Std: " << sqrt(m_variance) << endl;

        cout << "Acceptance rate: " << m_acceptanceRate << endl;
        m_averages << "Acceptance rate: " << m_acceptanceRate << endl;

        cout << " dE/dalpha       : " << m_dEdalpha << endl;
        cout << " dE/dbeta        : " << m_dEbeta   << endl;

        cout << "<r12>: " << m_meanDistance << endl;
        m_averages << "<r12>: " << m_meanDistance << endl;

        cout << "<K>: " << m_meanKineticEnergy << endl;
        m_averages << "<K>: " << m_meanKineticEnergy << endl;

        cout << "<V>: " << m_meanPotentialEnergy;
        m_averages << "<V>: " << m_meanPotentialEnergy;

        cout << endl;

        m_averages.close();

    }
}

void Sampler::computeAverages() {

    m_energy   = m_cumulativeEnergy/(double)m_stepNumber;
    m_E2       = m_E2/(double)m_stepNumber;
    m_variance = (m_E2 - m_energy*m_energy)/(double)m_stepNumber;
    m_dEdalpha = 2.0*(m_dE1alpha/(double)m_stepNumber - (m_dE2alpha/(double)m_stepNumber)*m_energy);
    m_dEbeta   = 2.0*(m_dE1beta/(double)m_stepNumber  -  (m_dE2beta/(double)m_stepNumber)*m_energy);
    m_meanDistance = m_meanDistance/(double)m_stepNumber;
    m_meanKineticEnergy = m_meanKineticEnergy/(double)m_stepNumber;
    m_meanPotentialEnergy = m_meanPotentialEnergy/(double)m_stepNumber;
    m_acceptanceRate  = (double)m_accepted/(double)((1-m_system->getEquilibrationFraction())*m_numberOfMetropolisSteps);

    //Parallell communication

    double reducedEnergy  = 0;
    double reducedE2      = 0;
    double reducedKinetic = 0;
    double reducedPotential = 0;

    MPI_Reduce(&m_energy,&reducedEnergy,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_E2,&reducedE2,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_meanKineticEnergy,&reducedKinetic,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_meanPotentialEnergy,&reducedPotential,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (m_system->getRank()==0) {
        m_energy              = reducedEnergy/(double)m_system->getSize();
        m_meanKineticEnergy   = reducedKinetic/(double)m_system->getSize();
        m_meanPotentialEnergy = reducedPotential/(double)m_system->getSize();
        m_variance            = (reducedE2/(double)m_system->getSize() - m_energy*m_energy)/((double)m_stepNumber*m_system->getSize());
    }

}

void Sampler::reset() {

    m_cumulativeEnergy = 0.0;
    m_accepted = 0;
    m_E2 = 0.0;

    m_dEdalpha = 0.0;
    m_dE1alpha = 0.0;
    m_dE2alpha = 0.0;

    m_dEbeta   = 0.0;
    m_dE1beta  = 0.0;
    m_dE2beta  = 0.0;

    m_variance = 0.0;
    m_meanKineticEnergy = 0.0;
    m_meanPotentialEnergy = 0.0;
    m_meanDistance = 0.0;

    m_stepNumber = 0;

}


