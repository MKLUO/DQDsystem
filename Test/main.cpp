#include "HilbertSpace.h"
#include "Plot.h"
#include "Fourier.h"

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;


int main()
{

	HilbertSpace hilbertSpace = HilbertSpace(40, 40, 0.1);

	// Build state

	SPState spState1 =
			hilbertSpace.createSingleParticleState(
					[](double x, double y) {
						return gaussian(x, y, 1.);
					});

	State dpState = (spState1 ^ spState1);

	Complex c = spState1 * spState1;

	return 0;
}
