#pragma once
#include <vector>

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual void intialize() = 0;
    virtual double computeRatio(int particle_nr) = 0;
    virtual void updateInverse(int particle_nr) = 0;
    virtual double computeLaplacian(int particle_nr) = 0;
    virtual std::vector<double> computeGradient(int particle_nr) = 0;
    virtual std::vector<double> computeGradientImportance(int particle_nr) = 0;
    virtual double compute_Jastrow(int particle_nr) = 0;
    virtual double compute_JastrowLaplacian() = 0;
    virtual void computeParametersDerivative() = 0;
    virtual double getAlphaDerivative() = 0;
    virtual double getBetaDerivative() = 0;
    virtual void setAlpha(double alpha) = 0;
    virtual void setBeta(double beta) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

