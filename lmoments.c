// Functions for computing the L-moments of a sample,
// plus a bunch of other auxiliar functions

// Function prototypes
void pwm(float series[], int n, float beta[], float A, float B);
void lMoments(float series[], int n, float lMoment[], float A, float B);

// pwm()
// Calculates the first three probability weighted moments of a sample,
// using either the direct method from the definition of pwm (when A=B=0)
// or a plotting position formula (when A<=B<=0).
void pwm(float series[], int n, float beta[], float A, float B) {

	int i;
	float acum[3], F;

	acum[0] = acum[1] = acum[2] = 0;
	if (A==0 && B==0) {
		for (i=1; i<=n; i++) {
			acum[0] = acum[0] + series[i-1];
			acum[1] = acum[1] + series[i-1] * (i-1) / (n-1);
			acum[2] = acum[2] + (series[i-1] * (i-1) * (i-2) / (n-1) / (n-2));
		}
	}
	else {
		for (i=1; i<=n; i++) {
			F = (i+A) / (n+B);
			acum[0] += series[i-1];
			acum[1] += series[i-1]*(1-F);
			acum[2] += series[i-1]*(1-F)*(1-F);
		}
	}
	beta[0] = acum[0] / n;
	beta[1] = acum[1] / n;
	beta[2] = acum[2] / n;
	//printf("\nbeta0: %.4f\n", beta[0]);
	//printf("beta1: %.4f\n", beta[1]);
	//printf("beta2: %.4f\n", beta[2]);
}

// lMoments()
// Estimates the first two L-moments of the sample
void lMoments(float series[], int n, float lMoment[], float A, float B) {

	int i, j, ordenLMom;
	float C, D, E, acum[3], beta[3];

	// Calculate the first three PWMs
	pwm(series, n, beta, A, B);

	// Obtain the first two L-moments
//	ordenLMom = 1;
//	lMoment[ordenLMom] = beta[0];
//	for (j=2; j<3; j++) {
//		acum[0] = 0;
//		ordenLMom = j;
//		for (i=0; i<ordenLMom; i++) {
//			C = pow(-1, ordenLMom-1-i);
//			D = 1;
//			E = 1;
//			if (i>0) {
//				if (i==ordenLMom-1) {
//					E = factorial(ordenLMom-1+i) / factorial(i) / factorial(ordenLMom-1);
//				}
//				else {
//					D = factorial(ordenLMom-1) / factorial(i) / factorial(ordenLMom-1-i);
//					E = factorial(ordenLMom-1+i) / factorial(i) / factorial(ordenLMom-1);
//				}
//			}
//			acum[0] = acum[0] + (C * D * E * beta[i]);
//			lMoment[ordenLMom] = acum[0];
//		}
//	}
//	if (A!=0 || B!=0) lMoment[2]=-lMoment[2];
	lMoment[1] = beta[0];
	lMoment[2] = beta[0] - 2*beta[1];
	lMoment[3] = beta[0] - 6*beta[1] + 6*beta[2];
	//lMoment[4] = beta[0] - 12*beta[1] + 30*beta[2] - 20*beta[3];
	//printf("\nL-moment order 1 = %.4f\n", lMoment[1]);
	//printf("L-moment order 2 = %.4f\n", lMoment[2]);
	//printf("L-moment order 3 = %.4f\n", lMoment[3]);
}

