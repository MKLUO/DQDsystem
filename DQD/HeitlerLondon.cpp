#include <functional>
using namespace std::placeholders;

#include "HeitlerLondon.h"

double calculateJWithSetting_HL(const Setting & setting) {

    HilbertSpace hilbertSpace = HilbertSpace(100, 50);

    auto fn1 = bind(fockDarwin, _1, _2, setting);

    HilbertSpace::SingleParticleState FDState_left =
            hilbertSpace.createSingleParticleState(ScalarFunction(&fockDarwin));
}

