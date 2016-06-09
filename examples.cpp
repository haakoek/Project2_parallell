#include <mpi/mpi.h>
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
#include <WaveFunctions/slaterwavefunction.h>
#include <WaveFunctions/wavefunction.h>
#include <SingleParticleWaveFunctions/singleparticleharmonicoscillator.h>
#include <SingleParticleWaveFunctions/singleparticlewavefunctions.h>
#include "Hamiltonians/simplequantumdothamiltonian.h"
#include <cstdlib>

using namespace std;

int Examples::TwoParticles(int argc, char* argv[]) {

    MPI_Init (&argc, &argv);	/* starts MPI */
    int rank, size;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */

    int numberOfDimensions = 2;
    int numberOfParticles  = 2;
    int numberOfSteps = (int) 1e7;
    double stepLength = 1.7;
    double timeStep = 0.0001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool imp_sampling = true;
    bool optimize     = false;
    bool writeToFile  = false;

    double omega, alpha, beta;
    int int_jastrow;
    bool jastrow;

    if(argc > 1) {
        omega = atof(argv[1]);
        alpha = atof(argv[2]);
        beta  = atof(argv[3]);
        int_jastrow = atoi(argv[4]);

        if(int_jastrow == 1) {
            jastrow = true;
        } else {
            jastrow = false;
        }

    }

    System* system =                     new System(imp_sampling, jastrow,writeToFile);
    system->setRank(rank);
    system->setSize(size);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);

    if(optimize) {

        SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

        if(jastrow) {
            std::vector<double> parameters(2);
            parameters[0] = alpha;
            parameters[1] = beta;
            sd->altOptimize(parameters);
        } else {
            sd->optimize(alpha);
        }

    } else {


        clock_t begin = clock();

        system->runMetropolisSteps(numberOfSteps);

        clock_t end   = clock();

        system->printToTerminal();

        cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    }

    MPI_Finalize();

    return 0;

}

int Examples::SixParticles(int argc, char* argv[]) {

    int numberOfDimensions = 2;
    int numberOfParticles  = 6;
    int numberOfSteps = (int) 1e5;
    double omega = 1.0;
    double alpha            = 0.841995;          // variational parameter 1
    double beta             = 0.663797;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.0004;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);

    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e5,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;

}

int Examples::TwelveParticles(int argc, char* argv[]) {


    int numberOfDimensions = 2;
    int numberOfParticles  = 12;
    int numberOfSteps = (int) 1e5;
    double omega = 1.0;
    double alpha            = 0.85;          // variational parameter 1
    double beta             = 0.65;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);


    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e5,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;


}

int Examples::TwentyParticles(int argc, char *argv[]) {

    int numberOfDimensions = 2;
    int numberOfParticles  = 20;
    int numberOfSteps = (int) 1e6;
    double omega = 1.0;
    double alpha            = 0.841677;          // variational parameter 1
    double beta             = 0.663795;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);

    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    //sd->altOptimize(parameters);

    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;


    return 0;

}
