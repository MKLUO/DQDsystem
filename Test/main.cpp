#include "HilbertSpace.h"

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::SingleParticleStatePair SPStatePair;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;

Complex gaussian(double x, double y) {
	return exp(-(x*x + y*y));
}

int main()
{
	HilbertSpace hilbertSpace = HilbertSpace(40, 40, 0.1);

	// Build H-L singlet/triplet state

	SPState spState =
			hilbertSpace.createSingleParticleState(
					[](double x, double y) {
						return gaussian(x, y);
					});

	State dpState = State({spState ^ spState}) + State({spState ^ spState}) * (-1.0) ;

	Complex c = dpState * dpState;

	return 0;
}
