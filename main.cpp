#include "examples.h"
#include <cstdlib>

int main(int argc, char* argv[]) {

    int nr_Particles = atoi(argv[5]);

    if(nr_Particles == 2) {
        return Examples::TwoParticles(argc, argv);
    } else if(nr_Particles == 6) {
        return Examples::SixParticles(argc, argv);
    } else if(nr_Particles == 12) {
        return Examples::TwelveParticles(argc, argv);
    } else if(nr_Particles == 20) {
        return Examples::TwentyParticles(argc, argv);
    } else {
        return 0;
    }
}
