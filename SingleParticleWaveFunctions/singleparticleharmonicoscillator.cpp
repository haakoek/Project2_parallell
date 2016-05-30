#include "singleparticleharmonicoscillator.h"
#include "../particle.h"
#include "Math/hermitepolynomials.h"
#include "cmath"
#include "iostream"
#include "singleparticlewavefunctions.h"

using namespace std;

SingleParticleHarmonicOscillator::SingleParticleHarmonicOscillator(double omega, double alpha)
{
    m_omega = omega;
    m_alpha = alpha;
}

double SingleParticleHarmonicOscillator::evaluateSp(int nx, int ny, Particle* particle) {

    double x = particle->getPosition()[0];
    double y = particle->getPosition()[1];

    double Hnx = HermitePolynomials::evaluate(nx,x,m_omega, m_alpha);
    double Hny = HermitePolynomials::evaluate(ny,y,m_omega, m_alpha);

    return Hnx*Hny*exp(-m_omega*m_alpha*0.5*(x*x+y*y));
}

std::vector<double> SingleParticleHarmonicOscillator::computeGradient(int nx, int ny, Particle *particle) {
    std::vector<double> gradient(2);
    gradient[0] = compute_dx(nx,ny,particle);
    gradient[1] = compute_dy(nx,ny,particle);
    return gradient;
}

double SingleParticleHarmonicOscillator::compute_dx(int nx,int ny, Particle *particle) {

    const double x = particle->getPosition()[0];
    const double y = particle->getPosition()[1];

    const double Hnx = HermitePolynomials::evaluate(nx,x,m_omega,m_alpha);
    const double Hny = HermitePolynomials::evaluate(ny,y,m_omega,m_alpha);

    const double dHnx = HermitePolynomials::evaluateDerivative(nx,x,m_omega,m_alpha);

    double r2 = x*x+y*y;

    return exp(-0.5*m_omega*m_alpha*r2) * Hny *
            ( dHnx - Hnx*m_omega*m_alpha*x);

}

double SingleParticleHarmonicOscillator::compute_dy(int nx,int ny, Particle *particle) {

    const double x = particle->getPosition()[0];
    const double y = particle->getPosition()[1];

    const double Hnx = HermitePolynomials::evaluate(nx,x,m_omega,m_alpha);
    const double Hny = HermitePolynomials::evaluate(ny,y,m_omega,m_alpha);


    const double dHny = HermitePolynomials::evaluateDerivative(ny,y,m_omega,m_alpha);

    double r2 = x*x+y*y;

    return exp(-0.5*m_omega*m_alpha*r2) * Hnx *
            ( dHny - Hny*m_omega*m_alpha*y);
}

double SingleParticleHarmonicOscillator::computeLaplacian(int nx, int ny, Particle *particle) {

    double x = particle->getPosition()[0];
    double y = particle->getPosition()[1];

    const double Hnx = HermitePolynomials::evaluate(nx,x,m_omega,m_alpha);
    const double Hny = HermitePolynomials::evaluate(ny,y,m_omega,m_alpha);

    const double dHnx = HermitePolynomials::evaluateDerivative(nx,x,m_omega,m_alpha);
    const double dHny = HermitePolynomials::evaluateDerivative(ny,y,m_omega,m_alpha);

    const double ddHnx = HermitePolynomials::evaluateDoubleDerivative(nx,x,m_omega,m_alpha);
    const double ddHny = HermitePolynomials::evaluateDoubleDerivative(ny,y,m_omega,m_alpha);

    double r2 = x*x + y*y;

    double m_omegaAlpha = m_omega*m_alpha;

    return exp(-0.5*m_omegaAlpha*r2) *
               ( - 2*m_omegaAlpha*x*Hny*dHnx
                 - 2*m_omegaAlpha*y*Hnx*dHny
                 + m_omegaAlpha*Hnx*Hny*(m_omegaAlpha*r2-2)
                 + Hny*ddHnx
                 + Hnx*ddHny );
}

double SingleParticleHarmonicOscillator::compute_dalpha(int nx, int ny, Particle* particle) {

    double x = particle->getPosition()[0];
    double y = particle->getPosition()[1];
    double r2 = x*x+y*y;

    double Hnx = HermitePolynomials::evaluate(nx,x,m_omega,m_alpha);
    double Hny = HermitePolynomials::evaluate(ny,y,m_omega,m_alpha);
    double dalpha_Hnx = HermitePolynomials::evaluateAlphaDerivative(nx,x,m_omega,m_alpha);
    double dalpha_Hny = HermitePolynomials::evaluateAlphaDerivative(ny,y,m_omega,m_alpha);

    return exp(-0.5*m_omega*m_alpha*r2)*(dalpha_Hnx*Hny + Hnx*dalpha_Hny - 0.5*m_omega*r2*Hnx*Hny);

}
