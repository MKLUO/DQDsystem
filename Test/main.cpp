#include "HilbertSpace.h"
#include "Plot.h"
#include "Fourier.h"


int main() {

    HilbertSpace hilbertSpace = HilbertSpace(40, 40, 0.1);

    // Build state

    SPState spState1 =
            hilbertSpace.createSingleParticleState(
                    gaussian(.5));

    State dpState = (spState1 ^ spState1) + (spState1 ^ spState1);

    Complex c = dpState * dpState;

    fourier::test();

    ScalarField field(40, 40, 0.1, gaussian(.5));

    std::vector<Complex> field_data = field.getDatas();
    std::vector<Complex> field_FT =
            fourier::fft2d(field_data, field.getWidth(), field.getHeight());
    std::vector<Complex> field_FT_IFT =
            fourier::ifft2d(field_FT, field.getWidth(), field.getHeight());

    return 0;
}
