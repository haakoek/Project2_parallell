#include "examples.h"
#include <iostream>
#include <iomanip>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <ctime>
#include "steepestdescent.h"
#include <cmath>

using namespace std;

SteepestDescent::SteepestDescent(System *system, int nSteps, int n_particles, int n_dim)
{
    m_system = system;
    m_nSteps = nSteps;
    m_iterations = 0;
    m_nparticles = n_particles;
    m_numberOfDimensions = n_dim;
}

void SteepestDescent::optimize(std::vector<double> parameters) {

    int maxNumberOfSteps = 100;
    int stepNumber = 0;
    int numberOfParameters = parameters.size();
    double tolerance = 1e-6;
    double m_stepLengthOptimize = 0.0001;
    double oldAbsoluteGradient = 1;

    while (oldAbsoluteGradient > tolerance && stepNumber < maxNumberOfSteps) {

        cout << "****** Iteration " << stepNumber << " ******" << endl;

        m_system->setInitialState(new RandomUniform(m_system, m_numberOfDimensions, m_nparticles));

        m_system->setAlpha(parameters[0]);
        m_system->setBeta(parameters[1]);

        m_system->runMetropolisSteps((int) m_nSteps);

        std::vector<double> localEnergyGradient(2);

        localEnergyGradient[0] = m_system->getAlphaDerivativeEnergy();
        localEnergyGradient[1] = m_system->getBetaDerivativeEnergy();

        double newAbsoluteGradient = 0;

        for (int j=0; j < 2; j++) {
            newAbsoluteGradient += localEnergyGradient[j]*localEnergyGradient[j];
        }

        newAbsoluteGradient = sqrt(newAbsoluteGradient);

        cout << "Gradient: " << newAbsoluteGradient << endl;

        // update alpha and beta if new gradient is less than old
        // if not, reduce step length
        if (newAbsoluteGradient < oldAbsoluteGradient) {
            for (int j=0; j < 2; j++) {
                parameters[j] -= m_stepLengthOptimize*localEnergyGradient[j];
            }
        }
        else {
            m_stepLengthOptimize /= 1.2;
            cout << "New step length: " << m_stepLengthOptimize << endl;
        }

        for (int j=0; j < 2; j++) {
            cout << " Parameter " << j+1 << " : " << parameters.at(j) << endl;
        }
        cout << endl;

        // before new iteration
        oldAbsoluteGradient = newAbsoluteGradient;
        stepNumber++;

    }

    for (int j=0; j < numberOfParameters; j++) {
        cout << " Optimal parameter " << j+1 << " : " << parameters.at(j) << endl;
    }

    cout << endl;

    m_optimalAlpha = parameters[0];
    m_optimalBeta  = parameters[1];

}

void SteepestDescent::altOptimize(std::vector<double> parameters) {

    double E_new;
    double alpha_new;
    double beta_new;

    double GAMMA = 0.001;    //step size (learning rate)
    int MAX_ITER = 100;  //maximum number of iterations
    //double FUNC_TOL = 0.1;  //termination tolerance for F(x)

    int iter = 0;
    int iters_sinceLastImprovement = 0;

    double alpha_old = parameters[0];
    double beta_old  = parameters[1];

    m_system->setInitialState(new RandomUniform(m_system, m_numberOfDimensions, m_nparticles));

    m_system->setAlpha(alpha_old);
    m_system->setBeta(beta_old);

    m_system->runMetropolisSteps((int) m_nSteps);

    std::vector<double> gradient_old(2);
    std::vector<double> gradient_new(2);

    gradient_old[0] = m_system->getAlphaDerivativeEnergy();
    gradient_old[1] = m_system->getBetaDerivativeEnergy();

    double E_old = m_system->getEnergy();

    cout << E_old << endl;

    while(iter < MAX_ITER) {

        iter+=1;

        alpha_new = alpha_old-GAMMA*gradient_old[0];
        beta_new  = beta_old-GAMMA*gradient_old[1];

        m_system->setInitialState(new RandomUniform(m_system, m_numberOfDimensions, m_nparticles));

        m_system->setAlpha(alpha_new);
        m_system->setBeta(beta_new);

        m_system->runMetropolisSteps((int) m_nSteps);

        E_new = m_system->getEnergy();

        gradient_new[0] = m_system->getAlphaDerivativeEnergy();
        gradient_new[1] = m_system->getBetaDerivativeEnergy();

        cout << "************" << endl;
        cout << "Iteration: " << iter  << endl;
        cout << "E_old:     " << E_old << endl;
        cout << "E_new:     " << E_new << endl;

        if(E_new < E_old) {

            iters_sinceLastImprovement = 0;

            cout << "alpha_new: " << alpha_new << endl;
            cout << "beta_new:  " << beta_new  << endl;
            cout << "************" << endl;

            E_old = E_new;
            alpha_old = alpha_new;
            beta_old  = beta_new;

            gradient_old[0] = gradient_new[0];
            gradient_old[1] = gradient_new[1];


        } else {


            if(iters_sinceLastImprovement == 20) {

                cout << "Reset gamma" << endl;

                GAMMA = 0.5;
                iters_sinceLastImprovement = 0;

            } else {

                GAMMA /= 1.2;
                iters_sinceLastImprovement += 1;

            }

            cout << "************" << endl;

        }


    }




    m_optimalAlpha = alpha_old;
    m_optimalBeta  = beta_old;

    cout << "Optimal alpha: " << m_optimalAlpha << endl;
    cout << "Optimal beta:  " << m_optimalBeta  << endl;


}

double SteepestDescent::getOptAlpha() {
    return m_optimalAlpha;
}

double SteepestDescent::getOptBeta() {
    return m_optimalBeta;
}

int SteepestDescent::getNrOfIters() {
    return m_iterations;
}

