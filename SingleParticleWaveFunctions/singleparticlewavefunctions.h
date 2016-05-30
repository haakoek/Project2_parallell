#pragma once
#include <vector>

class SingleParticleWaveFunctions
{
public:
    virtual double evaluateSp(int nx, int ny, class Particle* particle) = 0;
    virtual double computeLaplacian(int nx, int ny, class Particle* particle) = 0;
    virtual std::vector<double> computeGradient(int nx, int ny, class Particle* particle) = 0;
    virtual double compute_dalpha(int nx,int ny, class Particle* particle) = 0;
    virtual void setAlpha(double alpha) = 0;
    //virtual double compute_dx(int k, double omega, class Particle* particle) = 0;
    //virtual double compute_dy(int k, double omega, class Particle* particle) = 0;
};

