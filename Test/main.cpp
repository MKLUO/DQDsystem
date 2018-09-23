#include "HilbertSpace.h"

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::SingleParticleStatePair SPStatePair;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;

Complex gaussian(double x, double y) {
	return exp(-(x*x + y*y)/8.0);
}

int main()
{
	HilbertSpace hilbertSpace = HilbertSpace(20, 20, 0.1);

	// Build H-L singlet/triplet state

	SPState spState =
			hilbertSpace.createSingleParticleState(
					[](double x, double y) {
						return gaussian(x, y);
					});

	State dpState = State({spState ^ spState});

	return 0;
}
