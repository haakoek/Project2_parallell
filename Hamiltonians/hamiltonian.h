#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system, bool jastrow);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual double computeKineticEnergy();
    virtual double getKineticEnergy(std::vector<class Particle*> particles) = 0;
    virtual double getPotentialEnergy(std::vector<class Particle*> particles) = 0;
protected:
    class System* m_system = nullptr;
private:
    bool m_jastrow = true;
};

