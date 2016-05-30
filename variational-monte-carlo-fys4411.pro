TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt
QMAKE_CXXFLAGS += -std=c++11


SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    examples.cpp \
    steepestdescent.cpp \
    SingleParticleWaveFunctions/singleparticlewavefunctions.cpp \
    SingleParticleWaveFunctions/singleparticleharmonicoscillator.cpp \
    Math/hermitepolynomials.cpp \
    WaveFunctions/slaterwavefunction.cpp \
    Hamiltonians/simplequantumdothamiltonian.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    examples.h \
    steepestdescent.h \
    SingleParticleWaveFunctions/singleparticlewavefunctions.h \
    SingleParticleWaveFunctions/singleparticleharmonicoscillator.h \
    Math/hermitepolynomials.h \
    WaveFunctions/slaterwavefunction.h \
    Hamiltonians/simplequantumdothamiltonian.h

LIBS += -llapack -lblas -larmadillo

INCLUDEPATH += /usr/include/openmpi-x86_64
QMAKE_CXX = /usr/bin/mpicxx #/usr/lib64/openmpi/bin/mpicxx #mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = /usr/bin/mpicxx #/usr/lib64/openmpi/bin/mpicc #mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

