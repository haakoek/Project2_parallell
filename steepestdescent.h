#pragma once
class SteepestDescent
{
public:
    SteepestDescent(class System* system, int nSteps, int n_particles, int n_dim);
    double getOptAlpha();
    double getOptBeta();
    void optimize(double alpha);
    void altOptimize(std::vector<double> parameters);
    int getNrOfIters();
private:
    class System* m_system = nullptr;
    int m_nSteps = 0;
    int m_iterations = 0;
    int m_nparticles = 0;
    int m_numberOfDimensions = 0;
    double m_optimalAlpha = 0;
    double m_optimalBeta = 0;
};

