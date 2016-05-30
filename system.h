#pragma once
#include <vector>

class System {
public:
    System(bool importanceSampling, bool with_jastrow, bool writeToFile);
    bool metropolisStepSlaterDet();
    bool metropolisStepSlaterDetImportance();
    std::vector<double> quantumForce(int random_particle);
    double evaluateGreensFunction(std::vector<double> newPosition, std::vector<double> oldPosition,
                                          std::vector<double> quantumForce);
    void runMetropolisSteps                 (int numberOfMetropolisSteps);
    void setNumberOfParticles               (int numberOfParticles);
    void setNumberOfDimensions              (int numberOfDimensions);
    void setStepLength                      (double stepLength);
    void setTimeStepImportanceSampling      (double timeStep);
    void setEquilibrationFraction           (double equilibrationFraction);
    void setHamiltonian                     (class Hamiltonian* hamiltonian);
    void setWaveFunction                    (class WaveFunction* waveFunction);
    void setInitialState                    (class InitialState* initialState);
    class WaveFunction*                     getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*                      getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                          getSampler()        { return m_sampler; }
    std::vector<class Particle*>&           getParticles()      { return m_particles; }
    int getNumberOfParticles()              { return m_numberOfParticles; }
    int getNumberOfDimensions()             { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()        { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()       { return m_equilibrationFraction; }
    double getAlphaDerivativeEnergy();
    double getBetaDerivativeEnergy();
    void setAlpha(double alpha);
    void setBeta(double beta);
    void printToTerminal();
    void printParticlePositions();
    void setWriteToFile(bool writeToFile);
    double getEnergy();
    void setRank(int rank) {m_rank = rank;}
    void setSize(int size) {m_size = size;}
    int getRank() {return m_rank;}
    int getSize() {return m_size;}

private:
    int m_rank;
    int m_size;
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_timeStep   = 0.005;
    bool                            m_importanceSampling = false;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    bool m_with_jastrow = true;
    bool m_writeToFile  = false;
};

