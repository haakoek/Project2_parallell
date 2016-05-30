#include "slaterwavefunction.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "SingleParticleWaveFunctions/singleparticlewavefunctions.h"
#include <armadillo>

using namespace std;
using namespace arma;

SlaterWaveFunction::SlaterWaveFunction(System* system,double omega, double alpha, double beta, SingleParticleWaveFunctions* spwf) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta  >= 0);
    assert(omega >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(omega);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_spwf = spwf;
    this->intialize();
}

void SlaterWaveFunction::intialize() {

    m_quantumNumbers = zeros<mat>(10, 2);

    // quantum numbers for up to 20 particles
    m_quantumNumbers(0,0) = 0; m_quantumNumbers(0,1) = 0;
    m_quantumNumbers(1,0) = 1; m_quantumNumbers(1,1) = 0;
    m_quantumNumbers(2,0) = 0; m_quantumNumbers(2,1) = 1;
    m_quantumNumbers(3,0) = 2; m_quantumNumbers(3,1) = 0;
    m_quantumNumbers(4,0) = 1; m_quantumNumbers(4,1) = 1;
    m_quantumNumbers(5,0) = 0; m_quantumNumbers(5,1) = 2;
    m_quantumNumbers(6,0) = 3; m_quantumNumbers(6,1) = 0;
    m_quantumNumbers(7,0) = 2; m_quantumNumbers(7,1) = 1;
    m_quantumNumbers(8,0) = 1; m_quantumNumbers(8,1) = 2;
    m_quantumNumbers(9,0) = 0; m_quantumNumbers(9,1) = 3;


    int N = m_system->getNumberOfParticles();

    m_SlaterMatrixUp = zeros<mat>(N/2,N/2);
    m_SlaterMatrixDown = zeros<mat>(N/2,N/2);

    for(int i = 0; i < N/2; i++) {
        for(int j = 0; j < N/2; j++) {
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            Particle* xUp = m_system->getParticles()[i];
            Particle* xDown = m_system->getParticles()[i+N/2];
            m_SlaterMatrixUp(i,j)   = m_spwf->evaluateSp(nx,ny,xUp);
            m_SlaterMatrixDown(i,j) = m_spwf->evaluateSp(nx,ny,xDown);
        }
    }

    m_SlaterInverseUp   = m_SlaterMatrixUp.i();
    m_SlaterInverseDown = m_SlaterMatrixDown.i();

    m_a = zeros<mat>(N,N);

    //Jastrow part

    for (int i=0; i < N; i++) {
        for (int j=0; j < N; j++) {
            if (i < N/2) {
                if (j < N/2) {
                    m_a(i,j) = 1.0/3;
                }
                else {
                    m_a(i,j) = 1;
                }
            }
            else {
                if (j < N/2) {
                    m_a(i,j) = 1;
                }
                else {
                    m_a(i,j) = 1.0/3;
                }
            }
        }
    }
}

double SlaterWaveFunction::compute_Jastrow(int particle_nr) {

    double jastrow_factor = 0.0;
    double beta = m_parameters[2];

    double xi = m_system->getParticles()[particle_nr]->getPosition()[0];
    double yi = m_system->getParticles()[particle_nr]->getPosition()[1];

    for(int j = 0; j < m_system->getNumberOfParticles(); j++) {


        if(j != particle_nr) {

            const double xj = m_system->getParticles()[j]->getPosition()[0];
            const double yj = m_system->getParticles()[j]->getPosition()[1];
            double r_ij = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);
            r_ij = sqrt(r_ij);

            jastrow_factor += (m_a(particle_nr,j)*r_ij)/(1.0+beta*r_ij);

        }
    }

    return exp(jastrow_factor);
}

double SlaterWaveFunction::computeRatio(int particle_nr) {

    m_i = particle_nr;

    Particle* xi = m_system->getParticles()[particle_nr];
    double R = 0.0;

    int N = m_system->getNumberOfParticles();

    if(particle_nr < N/2) {
        //In the case of moving a particle "with spin up"
        for(int k = 0; k < N/2; k++) {
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            R += m_spwf->evaluateSp(nx,ny,xi)*m_SlaterInverseUp(k,particle_nr);
        }
    } else {
        //In the case of moving a particle "with spin down"
        for(int k = 0; k < N/2; k++) {
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            R += m_spwf->evaluateSp(nx,ny,xi)*m_SlaterInverseDown(k,particle_nr-N/2);
        }
    }

    m_Ratio = R;

    return R;

}

void SlaterWaveFunction::updateInverse(int particle_nr) {

    Particle* xi = m_system->getParticles()[particle_nr];
    int N = m_system->getNumberOfParticles();
    double R = m_Ratio; //*jastrow_ratio; //Jastrow must be included here

    if(particle_nr < N/2) {
        //Update SlaterInverseUp

        //Elements of i-th column
        mat m_SlaterInverseUpOld = m_SlaterInverseUp;
        for(int k = 0; k < N/2; k++) {
            m_SlaterInverseUp(k,particle_nr) = (1.0/R)*m_SlaterInverseUpOld(k,particle_nr);
        }

        //Elements of j-th column

        for(int j = 0; j < N/2; j++) {

            if(j != particle_nr) {
                double Sj = 0.0;

                for(int l = 0; l < N/2; l++) {
                    int nx = m_quantumNumbers(l,0);
                    int ny = m_quantumNumbers(l,1);
                    Sj += m_spwf->evaluateSp(nx,ny,xi)*m_SlaterInverseUpOld(l,j);
                }

                for(int k = 0; k < N/2; k++) {
                    m_SlaterInverseUp(k,j) = m_SlaterInverseUpOld(k,j) - (Sj/R)*m_SlaterInverseUpOld(k,particle_nr);
                }

            }
        }

    } else {
        //Update SlaterInverseDown

        //Elements of i-th column
        mat m_SlaterInverseDownOld = m_SlaterInverseDown;
        for(int k = 0; k < N/2; k++) {
            m_SlaterInverseDown(k,particle_nr-N/2) = (1.0/R)*m_SlaterInverseDownOld(k,particle_nr-N/2);
        }

        //Elements of j-th column

        for(int j = 0; j < N/2; j++) {

            if(j != (particle_nr-N/2)) {
                double Sj = 0.0;

                for(int l = 0; l < N/2; l++) {
                    int nx = m_quantumNumbers(l,0);
                    int ny = m_quantumNumbers(l,1);
                    Sj += m_spwf->evaluateSp(nx,ny,xi)*m_SlaterInverseDownOld(l,j);
                }

                for(int k = 0; k < N/2; k++) {
                    m_SlaterInverseDown(k,j) = m_SlaterInverseDownOld(k,j) - (Sj/R)*m_SlaterInverseDownOld(k,particle_nr-N/2);
                }

            }
        }

    }

    //Jastrow part?

}

std::vector<double> SlaterWaveFunction::computeGradient(int particle_nr) {

    std::vector<double> grad_ri(m_system->getNumberOfDimensions());
    Particle* xi = m_system->getParticles()[particle_nr];
    int N = m_system->getNumberOfParticles();

    if(particle_nr < N/2) {
        for(int k = 0; k < N/2; k++) {
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            std::vector<double> grad_phij_ri = m_spwf->computeGradient(nx,ny,xi);
            for(int i = 0; i < m_system->getNumberOfDimensions(); i++) {
                grad_ri[i] += grad_phij_ri[i]*m_SlaterInverseUp(k,particle_nr);
            }
        }
    } else {
        for(int k = 0; k < N/2; k++) {
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            std::vector<double> grad_phij_ri = m_spwf->computeGradient(nx,ny,xi);
            for(int i = 0; i < m_system->getNumberOfDimensions(); i++) {
                grad_ri[i] += grad_phij_ri[i]*m_SlaterInverseDown(k,particle_nr-N/2);
            }
        }
    }


    return grad_ri;
}

double SlaterWaveFunction::computeLaplacian(int particle_nr) {

    //Compute laplacian of particle k!

    int N = m_system->getNumberOfParticles();

    if(m_firstStepLaplacian) {
        for(int i = 0; i < N/2; i++) {
            for(int j = 0; j < N/2; j++) {
                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);
                Particle* xUp   = m_system->getParticles()[i];
                Particle* xDown = m_system->getParticles()[i+N/2];
                lapUp   += m_spwf->computeLaplacian(nx,ny,xUp)*m_SlaterInverseUp(j,i);
                lapDown += m_spwf->computeLaplacian(nx,ny,xDown)*m_SlaterInverseDown(j,i);
            }
        }
        m_firstStepLaplacian = false;
    } else {

        if(m_i < N/2) {
            lapUp = 0.0;
            for(int i = 0; i < N/2; i++) {
                for(int j = 0; j < N/2; j++) {
                    int nx = m_quantumNumbers(j,0);
                    int ny = m_quantumNumbers(j,1);
                    Particle* xUp   = m_system->getParticles()[i];
                    lapUp   += m_spwf->computeLaplacian(nx,ny,xUp)*m_SlaterInverseUp(j,i);
                }
            }
        } else {
            lapDown = 0.0;
            for(int i = 0; i < N/2; i++) {
                for(int j = 0; j < N/2; j++) {
                    int nx = m_quantumNumbers(j,0);
                    int ny = m_quantumNumbers(j,1);
                    Particle* xDown = m_system->getParticles()[i+N/2];
                    lapDown += m_spwf->computeLaplacian(nx,ny,xDown)*m_SlaterInverseDown(j,i);
                }
            }
        }
    }

    return lapUp+lapDown;

}

std::vector<double> SlaterWaveFunction::computeGradientImportance(int i) {


    std::vector<double> gradient(2);

    std::vector<double> gradSlater =  computeGradient(i);
    std::vector<double> gradJastrow = gradientJastrow(i);

    gradient[0] = (gradSlater[0] + gradJastrow[0])/m_Ratio; //gradSlater[0]/m_Ratio;
    gradient[1] = (gradSlater[1] + gradJastrow[1])/m_Ratio; //gradSlater[1]/m_Ratio;

    return gradient;
}


double SlaterWaveFunction::compute_JastrowLaplacian() {

    //Jastrow part?
    int N = m_system->getNumberOfParticles();
    double beta = m_parameters[2];
    double jastrow_laplacian = 0.0;

    for(int i = 0; i < N; i++) {

        std::vector<double> gradJastrow = gradientJastrow(i);
        jastrow_laplacian += gradJastrow[0]*gradJastrow[0] + gradJastrow[1]*gradJastrow[1];

            double x_i = m_system->getParticles()[i]->getPosition()[0];
            double y_i = m_system->getParticles()[i]->getPosition()[1];

            for (int k=0; k < N; k++) {

                if (k != i) {
                    double x_k = m_system->getParticles()[k]->getPosition()[0];
                    double y_k = m_system->getParticles()[k]->getPosition()[1];
                    double r_ki = sqrt( (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i) );

                    double factor = 1.0 / (1 + beta*r_ki);
                    jastrow_laplacian += ( (m_a(k,i)*factor*factor) / r_ki ) -
                                        2*m_a(k,i)*beta*factor*factor*factor;
                }
        }

    }

    // cross term
    double slaterJastrow = 0;
    for (int i=0; i < N; i++) {
        std::vector<double> gradSlater  = computeGradient(i);
        std::vector<double> gradJastrow = gradientJastrow(i);
        slaterJastrow += gradSlater[0]*gradJastrow[0];
        slaterJastrow += gradSlater[1]*gradJastrow[1];
    }

    return jastrow_laplacian + 2.0*slaterJastrow;
}

std::vector<double> SlaterWaveFunction::gradientJastrow(int i) {
    // compute graident jastrow ratio w.r.t. to particle number i

    // compute jastrow factor
    std::vector<double> ratioJastrow(2);
    ratioJastrow[0] = 0; ratioJastrow[1] = 0;
    double beta = m_parameters[2];
    double x_i = m_system->getParticles()[i]->getPosition()[0];
    double y_i = m_system->getParticles()[i]->getPosition()[1];

    for (int j=0; j < m_system->getNumberOfParticles(); j++) {
        if (j != i) {
            double x_j = m_system->getParticles()[j]->getPosition()[0];
            double y_j = m_system->getParticles()[j]->getPosition()[1];
            double r_ij = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
            double factor = 1.0 / ( r_ij*(1 + beta*r_ij)*(1 + beta*r_ij) );

            ratioJastrow[0] += (x_i - x_j)*m_a(i,j)*factor;
            ratioJastrow[1] += (y_i - y_j)*m_a(i,j)*factor;
        }
    }

    return ratioJastrow;
}

void SlaterWaveFunction::computeParametersDerivative() {

    int N = m_system->getNumberOfParticles();
    double beta  = m_parameters[2];

    double d_alphaUp = 0.0;
    double d_alphaDown = 0.0;
    double d_beta = 0.0;

    for(int i = 0; i < N/2; i++) {
        for(int j = 0; j < N/2; j++) {
            Particle* xUp = m_system->getParticles()[i];
            Particle* xDown = m_system->getParticles()[i+N/2];
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            d_alphaUp   += m_spwf->compute_dalpha(nx,ny, xUp)*m_SlaterInverseUp(i,j);
            d_alphaDown += m_spwf->compute_dalpha(nx,ny, xDown)*m_SlaterInverseDown(i,j);
        }
    }

    m_dAlpha = d_alphaUp+d_alphaDown;

    for(int i = 0; i < N; i++) {

        double xi = m_system->getParticles()[i]->getPosition()[0];
        double yi = m_system->getParticles()[i]->getPosition()[1];

        for(int j = i+1; j < N; j++) {

            double xj = m_system->getParticles()[j]->getPosition()[0];
            double yj = m_system->getParticles()[j]->getPosition()[1];

            double r_ij = sqrt( (xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) );

            d_beta += -r_ij*((m_a(i,j)*r_ij))/((1+beta*r_ij)*(1+beta*r_ij));

        }
    }

    m_dBeta = d_beta;

}

double SlaterWaveFunction::getAlphaDerivative() {
    return m_dAlpha;
}

double SlaterWaveFunction::getBetaDerivative() {
    return m_dBeta;
}

void SlaterWaveFunction::setAlpha(double alpha) {
    m_spwf->setAlpha(alpha);
    m_parameters[1] = alpha;
}

void SlaterWaveFunction::setBeta(double beta) {
    m_parameters[2] = beta;
}

