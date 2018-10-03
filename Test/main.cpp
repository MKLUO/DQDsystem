#include "HeitlerLondon.h"

int main() {
    Setting setting = Setting::defaultSetting();

    //double J = calculateJWithSetting_HL(setting);

    ScalarField fd_field(200, 400, 0.01, fockDarwin(setting, Orientation::Left));

    Complex norm = fd_field * fd_field;

    return 0;
}