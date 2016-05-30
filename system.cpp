#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;

System::System(bool importanceSampling, bool with_jastrow, bool writeToFile) {
    m_importanceSampling = importanceSampling;
    m_with_jastrow = with_jastrow;
    m_writeToFile  = writeToFile;
}

bool System::metropolisStepSlaterDet() {

    int random_particle  = Random::nextInt(m_numberOfParticles);  //Choose a random particle
    int random_dimension = Random::nextInt(m_numberOfDimensions); //Choose a random dimension

    double Rc_old = 1;

    if(m_with_jastrow) {
        Rc_old = m_waveFunction->compute_Jastrow(random_particle);
    }

    double change = Random::nextDouble()*2.0-1.0;
    m_particles[random_particle]->adjustPosition(change*m_stepLength,random_dimension);

    double Rc_new = 1;

    if(m_with_jastrow) {
        Rc_new = m_waveFunction->compute_Jastrow(random_particle);
    }

    double R_sd = m_waveFunction->computeRatio(random_particle);
    double s = Random::nextDouble();
    double w = R_sd*R_sd*(Rc_new*Rc_new)/(Rc_old*Rc_old);
    //double jastrow_ratio = Rc_new/Rc_old;

    if(w > 1) {
        m_waveFunction->updateInverse(random_particle);
        return true;
    } else if(w >= s) {
        m_waveFunction->updateInverse(random_particle);
        return true;
    } else {
        m_particles[random_particle]->adjustPosition(-change*m_stepLength,random_dimension);
        return false;
    }

    //return true;

}

bool System::metropolisStepSlaterDetImportance() {

    int random_particle  = Random::nextInt(m_numberOfParticles);  //Choose a random particle

    std::vector<double> F_old = quantumForce(random_particle);

    double Rc_old = 1;

    if(m_with_jastrow) {
        Rc_old = m_waveFunction->compute_Jastrow(random_particle);
    }

    std::vector<double> change(m_numberOfDimensions);

    for(int i = 0; i < m_numberOfDimensions; i++) {
        change[i] = 0.5*F_old[i]*m_timeStep + Random::nextGaussian(0.0,sqrt(m_timeStep));
        m_particles[random_particle]->adjustPosition(change[i],i);
    }

    double Rc_new = 1;

    if(m_with_jastrow) {
        Rc_new = m_waveFunction->compute_Jastrow(random_particle);
    }

    double R_sd = m_waveFunction->computeRatio(random_particle);
    double ratio = R_sd*R_sd*(Rc_new*Rc_new)/(Rc_old*Rc_old);

    std::vector<double> F_new = quantumForce(random_particle);



    //Compute greens function
    std::vector<double> oldPosition(m_numberOfDimensions);
    std::vector<double> newPosition(m_numberOfDimensions);



    for(int i = 0; i < m_numberOfDimensions; i++) {
        newPosition[i] = m_particles[random_particle]->getPosition()[i];
        oldPosition[i] = m_particles[random_particle]->getPosition()[i] - change[i];
    }

    double exponent = 0;


    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        double term1 = - (oldPosition[dim] - newPosition[dim]  -  0.5*m_timeStep*F_new[dim]) *
                         (oldPosition[dim] - newPosition[dim]  -  0.5*m_timeStep*F_new[dim]);
        double term2 =   (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*F_old[dim]) *
                         (-oldPosition[dim] + newPosition[dim] - 0.5*m_timeStep*F_old[dim]);
        exponent += term1 + term2;
    }

    double greensRatio = exp(exponent / 2*m_timeStep);

    //double greens_old = evaluateGreensFunction(newPosition,oldPosition,F_old);
    //double greens_new = evaluateGreensFunction(oldPosition,newPosition,F_new);

    ratio = ratio*ratio*greensRatio;

    if(ratio >= Random::nextDouble()) {
        m_waveFunction->updateInverse(random_particle);
        return true;
    } else {

        for(int i = 0; i < m_numberOfDimensions; i++) {
            m_particles[random_particle]->adjustPosition(-change[i],i);
        }

        return false;
    }
}

std::vector<double> System::quantumForce(int random_particle) {

    std::vector<double> qForce = m_waveFunction->computeGradientImportance(random_particle);

    for(int i = 0; i < m_numberOfDimensions; i++) {
        qForce[i] *= 2.0;
    }

    return qForce;

}

double System::evaluateGreensFunction(std::vector<double> newPosition, std::vector<double> oldPosition,
                                      std::vector<double> quantumForce) {

    std::vector<double> greensVector;

    // make vector that needs to be dotted
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        greensVector.push_back(newPosition[dim] - oldPosition[dim] - 0.5*m_timeStep*quantumForce[dim]);
    }

    double greensFunction = 0;
    // find length squared of vector
    for (int dim=0; dim < m_numberOfDimensions; dim++) {
        greensFunction += greensVector[dim]*greensVector[dim];
    }
    greensFunction /= 2*m_timeStep;
    return exp(-greensFunction);
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {


    m_sampler                   = new Sampler(this, m_writeToFile);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    bool acceptedStep = false;

    for (int i=0; i < numberOfMetropolisSteps; i++) {

        if(m_importanceSampling) {
           acceptedStep = metropolisStepSlaterDetImportance();
        } else {
           acceptedStep = metropolisStepSlaterDet();
        }

        if (i > m_equilibrationFraction * m_numberOfMetropolisSteps) {
            m_sampler->sample(acceptedStep);
        }

        if(m_rank == 0) {
            if (!(i%1000)) {
                cout << "Progress: " << i/((double) numberOfMetropolisSteps) * 100 << " % \r";
                fflush(stdout);
            }
        }

    }

    m_sampler->computeAverages();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setTimeStepImportanceSampling(double timeStep) {
    assert(timeStep >= 0);
    m_timeStep = timeStep;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
    m_particles    = m_initialState->getParticles();
}

double System::getEnergy() {
    return m_sampler->getEnergy();
}

double System::getAlphaDerivativeEnergy() {
    return m_sampler->getDEDalpha();
}

double System::getBetaDerivativeEnergy() {
    return m_sampler->getDEDbeta();
}

void System::printToTerminal() {
    m_sampler->printOutputToTerminal();
}

void System::setAlpha(double alpha) {
    m_waveFunction->setAlpha(alpha);
}

void System::setBeta(double beta) {
    m_waveFunction->setBeta(beta);
}

void System::printParticlePositions() {
    for(int i = 0; i < m_numberOfParticles; i++) {
        cout << m_particles[i]->getPosition()[0];
    }
}

void System::setWriteToFile(bool writeToFile) {
    m_writeToFile = writeToFile;
    m_sampler->setWriteToFile(writeToFile);
}


