#pragma once
#include <iostream>
#include <fstream>
using namespace std;


class Sampler {
public:
    Sampler(class System* system, bool writeToFile);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void setWriteToFile(bool writeToFile);
    double getEnergy()          { return m_energy; }
    double getDEDalpha()        { return m_dEdalpha; }
    double getDEDbeta()         { return m_dEbeta;}
    void reset();
    double getAcceptanceRate() {return m_acceptanceRate;}

private:

    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    int     m_accepted = 0;

    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_variance = 0;
    double  m_E2 = 0;

    double  m_dEdalpha = 0.0;
    double  m_dE1alpha = 0.0;
    double  m_dE2alpha = 0.0;

    double  m_dEbeta   = 0.0;
    double  m_dE1beta  = 0.0;
    double  m_dE2beta  = 0.0;
    double m_probR[600];
    bool m_writeToFile = false;

    double m_meanKineticEnergy = 0.0;
    double m_meanPotentialEnergy = 0.0;
    double m_meanDistance = 0.0;
    double m_acceptanceRate = 0;
    class System* m_system = nullptr;
    ofstream m_outfile;
    ofstream m_positionsfile;
    ofstream m_averages;
};
