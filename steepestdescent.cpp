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

void SteepestDescent::optimize(double alpha) {

    double E_new;
    double alpha_new;
    //double beta_new;

    double gamma_alpha = 0.3;    //step size (learning rate)
    //double gamma_beta  = 0.3;

    int MAX_ITER = 100;  //maximum number of iterations
    int iters_sinceLastImprovement = 0;

    int iter = 0;

    double alpha_old = alpha;
    //double beta_old  = parameters[1];

    m_system->setAlpha(alpha_old);
    //m_system->setBeta(beta_old);

    m_system->runMetropolisSteps(m_nSteps);

    m_system->printToTerminal();

    double gradient_old, gradient_new;

    gradient_old = m_system->getAlphaDerivativeEnergy();
    //gradient_old[1] = m_system->getBetaDerivativeEnergy();

    double E_old = m_system->getEnergy();

    cout << E_old << endl;

    m_system->getSampler()->reset();

    while(iter < MAX_ITER && iters_sinceLastImprovement < 20) {

        iter+=1;

        cout << "***********" << endl;
        cout << "Iteration: " << iter << endl;

        alpha_new = alpha_old - gamma_alpha*gradient_old;

        //m_system->setInitialState(new RandomUniform(m_system, m_numberOfDimensions, m_nparticles));

        m_system->setAlpha(alpha_new);

        m_system->setInitialState(new RandomUniform(m_system,m_numberOfDimensions,m_nparticles));
        m_system->runMetropolisSteps((int) m_nSteps);

        cout << "Acceptance_new: " << m_system->getSampler()->getAcceptanceRate() << endl;

        E_new = m_system->getEnergy();

        gradient_new = m_system->getAlphaDerivativeEnergy();

        if(gradient_new*gradient_new > gradient_old*gradient_old) {
            if(iters_sinceLastImprovement == 10) {
                gamma_alpha = 0.03;
            } else {
                gamma_alpha *= 0.9;
            }
            iters_sinceLastImprovement += 1;
        } else {

            iters_sinceLastImprovement = 0;

            alpha_old = alpha_new;
            gradient_old = gradient_new;

            cout << "alpha_old: " << alpha_old << endl;
            cout << "alpha_new: " << alpha_new << endl;
            cout << "E_old: "     << E_old << endl;
            cout << "E_new: " << E_new << endl;

            E_old = E_new;

        }
        cout << "***********" << endl;
        m_system->getSampler()->reset();


    }

    m_optimalAlpha = alpha_old;

    cout << "Optimal alpha: " << m_optimalAlpha << endl;

    //return m_optimalAlpha;
}

void SteepestDescent::altOptimize(std::vector<double> parameters) {

    double E_new;
    double alpha_new;
    double beta_new;

    double gamma_alpha = 0.3;    //step size (learning rate)
    double gamma_beta  = 0.3;

    int MAX_ITER = 30;  //maximum number of iterations
    int iters_sinceLastImprovement = 0;

    int iter = 0;

    double alpha_old = parameters[0];
    double beta_old  = parameters[1];

    m_system->setAlpha(alpha_old);
    m_system->setBeta(beta_old);

    m_system->runMetropolisSteps((int) m_nSteps);

    m_system->printToTerminal();

    std::vector<double> gradient_old(2);
    std::vector<double> gradient_new(2);

    gradient_old[0] = m_system->getAlphaDerivativeEnergy();
    gradient_old[1] = m_system->getBetaDerivativeEnergy();

    double E_old = m_system->getEnergy();

    cout << E_old << endl;

    m_system->getSampler()->reset();

    while(iter < MAX_ITER && iters_sinceLastImprovement < 5) {

        iter+=1;

        cout << "***********" << endl;
        cout << "Iteration: " << iter << endl;

        alpha_new = alpha_old - gamma_alpha*gradient_old[0];
        beta_new  = beta_old  - gamma_beta*gradient_old[1];

        //m_system->setInitialState(new RandomUniform(m_system, m_numberOfDimensions, m_nparticles));

        m_system->setAlpha(alpha_new);
        m_system->setBeta(beta_new);

        m_system->runMetropolisSteps((int) m_nSteps);

        E_new = m_system->getEnergy();

        gradient_new[0] = m_system->getAlphaDerivativeEnergy();
        gradient_new[1] = m_system->getBetaDerivativeEnergy();

        if(gradient_new[0]*gradient_new[0]+gradient_new[1]*gradient_new[1] > gradient_old[0]*gradient_old[0]+gradient_old[1]*gradient_old[1]) {
            gamma_alpha *= 0.5;
            gamma_beta  *= 0.5;
            iters_sinceLastImprovement += 1;
        } else {

            iters_sinceLastImprovement = 0;

            alpha_old = alpha_new;
            beta_old  = beta_new;
            gradient_old[0] = gradient_new[0];
            gradient_old[1] = gradient_new[1];


            cout << "alpha_new: " << alpha_new << endl;
            cout << "beta_new:  " << beta_new  << endl;
            cout << "E_old: "     << E_old << endl;
            cout << "E_new: " << E_new << endl;

            E_old = E_new;

        }
        cout << "***********" << endl;
        m_system->getSampler()->reset();


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

