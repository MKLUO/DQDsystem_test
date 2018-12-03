#include "hl_hm.hpp"
#include "plot.hpp"

#include "fourier.hpp"

#include <iostream>
#include <iomanip>
#include <corecrt_math_defines.h>

#define PLOT_TRUNC 5

void test_HL();
void test_HM(Setting);
// void test2();
void test_FD(Setting);
// void test_CONV();
// void test_complex();
void test(Setting);



int main() {

    Setting setting = Setting::defaultSetting();

    setting.d = 100E-9;

    // test_HL();
    test_HM(setting);
    // test2();
    // test_FD(setting);
    // test_CONV();

    // test_complex();

    // test(setting);

    return 0;
}

void test_HL() {
    int cases = 10;
    double max_B = 4.;

    std::vector<double> Bs;

    Setting setting = Setting::defaultSetting();
    setting.d = 40E-9;
    for (int i = 0; i < cases; ++i)
        Bs.push_back(max_B / double(cases) * double(i));

    for (double B : Bs) {
        setting.B = B;
        std::cout << std::setprecision(std::numeric_limits<double>::digits10) << 
            B << " : " << calculateJWithSetting_HL(setting) << 
        std::endl;
    }
}

void test_HM(Setting setting) {

    Matrix result = coulombMatrix(setting, Ansatz::HM, Basis::FD);

    plotter::printMatrix(result);
}

void test2() {
    
    Setting setting = Setting::defaultSetting();

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    SPState scalarState =
        hilbertSpace.createSingleParticleState(scalar(1.));

    SPState gauss1 =
        hilbertSpace.createSingleParticleState(gaussian(10E-9));

    State state1 = (scalarState ^ scalarState).normalize();
    State state2 = (scalarState ^ scalarState).normalize() * 2;

    State state3 = (gauss1 ^ scalarState);


    Operator identityOp =
        hilbertSpace.createOperator(
            identity,
            identity);

    Operator laplaceOp1 =
        hilbertSpace.createOperator(
            laplacian,
            identity);

    ScalarField gauss_lap = laplacian(gauss1.getField());

    ComplexHighRes a = hilbertSpace.operatorValue(state2, identityOp + identityOp, state1);
    double b = hilbertSpace.expectationValue(state2, identityOp + identityOp);
    double c = hilbertSpace.expectationValue(state3, identityOp);
    double d = hilbertSpace.expectationValue(state3, laplaceOp1);

    plotter::plotField(gauss_lap, 10, "./FIELDCAR_gausslap");
    plotter::plotField(gauss1.getField(), 10, "./FIELDCAR_gauss");

    ScalarField wave1 = hilbertSpace.createScalarField(planeWave(   
        2. * M_PI * 1. / 200E-9,
        2. * M_PI * 0. / 2.));

    ScalarField wave1_lap = laplacian(wave1);

    plotter::plotField(wave1, 10, "./FIELDCAR_wave1");
    plotter::plotField(wave1_lap, 10, "./FIELDCAR_wave1_lap");
}

void test_FD(Setting setting) {

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    SPState left =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Left));

    SPState right =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Right));

    ScalarField left_sf = left.getField();
    ScalarField right_sf = right.getField();

    plotter::plotField(left_sf ^ left_sf, 
        PLOT_TRUNC, "./temp/FIELDCAR_FDL");
    plotter::plotField(right_sf, 
        PLOT_TRUNC, "./temp/FIELDCAR_FDR");
    plotter::plotField(left_sf ^ right_sf, 
        PLOT_TRUNC, "./temp/FIELDCAR_FDLR");         

    ComplexHighRes test1 = (left ^ right) * (left ^ right);

    Operator i_ts =
            hilbertSpace.createOperator(
                    identity_twoSite(setting));   
    
    // ComplexHighRes result = 
    //     hilbertSpace.operatorValue(
    //         left ^ right, 
    //         i_ts, 
    //         right ^ left
    //     );  

    return;                        
}

void test_CONV() {
    ScalarField c1 = ScalarField(130, 150, 1.0, scalar(1.0));
    ScalarField c2 = ScalarField(100, 100, 1.0, scalar(1.0));

    ScalarField c3 = fourier::convolution(c1, c2);

    plotter::plotField(c3, 10, "./temp/C3");
}

void test_complex() {
    ComplexHighRes zero;
    ComplexHighRes a(1.0 + 1.0i);

    auto b = a + a;
    auto c = (b + b + a + b * b) / a / b;
    auto d = c / c;
    auto e = c * c;
    auto f = e * e * e;
    auto g = f * f;

    auto h = g + g + g;

    return;
}

void test(Setting setting) {

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    SPState left_up    = OrthofockDarwinWithSpin(hilbertSpace, setting, Orientation::Left,  Spin::Up);
    SPState left_down  = OrthofockDarwinWithSpin(hilbertSpace, setting, Orientation::Left,  Spin::Down);
    SPState right_up   = OrthofockDarwinWithSpin(hilbertSpace, setting, Orientation::Right, Spin::Up);
    SPState right_down = OrthofockDarwinWithSpin(hilbertSpace, setting, Orientation::Right, Spin::Down);

    std::vector<State> states;

    states.push_back((left_up   ^ right_down).antisym()); // (up, down)
    states.push_back((left_down ^ right_up)  .antisym()); // (down, up)
    states.push_back((left_up   ^ left_down) .antisym()); // (2, 0)
    states.push_back((right_up  ^ right_down).antisym()); // (0, 2)

    ComplexHighRes a = left_up    * left_down;
    ComplexHighRes b = left_down  * left_up;
    ComplexHighRes c = right_up   * left_up;
    ComplexHighRes d = right_down * left_down;

    ComplexHighRes a2 = states[0] * states[1];
    ComplexHighRes b2 = states[1] * states[0];
    ComplexHighRes c2 = states[2] * states[3];
    ComplexHighRes d2 = states[3] * states[2];

    return;
}