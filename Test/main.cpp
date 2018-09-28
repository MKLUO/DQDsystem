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
			        gaussian(.5));

	State dpState = (spState1 ^ spState1) + (spState1 ^ spState1);

	Complex c = dpState * dpState;

	return 0;
}
