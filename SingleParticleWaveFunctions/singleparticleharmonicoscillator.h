#pragma once
#include "singleparticlewavefunctions.h"
#include <iostream>

using namespace std;

class SingleParticleHarmonicOscillator : public SingleParticleWaveFunctions
{
public:
    SingleParticleHarmonicOscillator(double omega, double alpha);
    double evaluateSp(int nx, int ny, class Particle *particle);
    double computeLaplacian(int nx, int ny, class Particle *particle);
    double compute_dx(int nx, int ny, class Particle* particle);
    double compute_dy(int nx, int ny, class Particle* particle);
    std::vector<double> computeGradient(int nx, int ny, class Particle* particle);
    double compute_dalpha(int nx, int ny, class Particle* particle);
    void setAlpha(double alpha) {
        //cout << "balle" << endl;
        m_alpha = alpha;}
private:
    double m_omega;
    double m_alpha;
};

