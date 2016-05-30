#pragma once
#include "wavefunction.h"
#include <armadillo>

using namespace arma;

class SlaterWaveFunction : public WaveFunction
{
public:
    SlaterWaveFunction(class System* system,double omega, double alpha, double beta, class SingleParticleWaveFunctions* spwf);
    void intialize();
    double computeRatio(int particle_nr);
    void updateInverse(int particle_nr);
    double computeLaplacian(int particle_nr);
    std::vector<double> computeGradient(int particle_nr);
    double compute_Jastrow(int particle_nr);
    std::vector<double> gradientJastrow(int i);
    std::vector<double> computeGradientImportance(int particle_nr);
    double compute_JastrowLaplacian();
    void computeParametersDerivative();
    double getAlphaDerivative();
    double getBetaDerivative();
    void setAlpha(double alpha);
    void setBeta(double beta);
    //void updateInverse()
    //double computeR()
    //double computeLaplacian()
    //std vector computeGradient()
private:
    class SingleParticleWaveFunctions* m_spwf = nullptr;
    mat m_SlaterInverseUp;
    mat m_SlaterInverseDown;
    mat m_SlaterMatrixUp;
    mat m_SlaterMatrixDown;
    mat m_a;
    mat m_quantumNumbers;
    double m_Ratio  = 1;
    double m_dAlpha = 0;
    double m_dBeta  = 0;
    int m_i = 0;
    bool m_firstStepLaplacian = true;
    double lapUp = 0.0;
    double lapDown = 0.0;
};



