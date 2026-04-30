#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace IPC
{

	namespace FEM {
		const static double PI = 3.1415926535897932;
		//const static double density = 1000;

		//const static double lengthRateLame = YoungModulus / (2 * (1 + PoissonRate));
		//const static double volumeRateLame = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
		//const static double lengthRate = 4 * lengthRateLame / 3;
		//const static double volumeRate = volumeRateLame + 5 * lengthRateLame / 6;
		//const static double friction = 0.5;
		//static double Hhat = 1e-6;
		//const static double IPC_dt = 0.01;
		//static double Kappa = 0;

		//const static double clothThicness = 1e-3;
		//const static double clothYoungModulus = 1e3;
		//const static double stretchStiffness = clothYoungModulus / (2 * (1 + PoissonRate));
		//const static double shearStiffness = stretchStiffness*0.3;
		//const static double clothDensity = 1e2;
	}
}

#endif // !FEM_PARAMETERS_H

