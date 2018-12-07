#include "hl_hm.hpp"
#include "plot.hpp"

#include "fourier.hpp"

#include <iostream>
#include <iomanip>

#define PLOT_TRUNC 5

using HL_HM::Setting;
using HL_HM::Ansatz;
using HL_HM::Basis;

void test_HM(Setting);

int main() {

    fourier::info();

    Setting setting = Setting::defaultSetting();

    setting.a = 10E-9;
    setting.d = 50E-9;
    setting.B = 1.0;

    // test_HL();
    test_HM(setting);
    // test2();
    // test_FD(setting);
    // test_CONV();

    // test_complex();

    // test(setting);

    return 0;
}

void test_HM(Setting setting) {

    Matrix result = coulombMatrix(setting, Ansatz::HM, Basis::FD);

    plotter::printMatrix(result);
}