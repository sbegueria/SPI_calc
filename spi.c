#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "auxiliary.c"
#include "lmoments.c"
#include "pdfs.c"

// Max size of raw rainfall and events matrices
#define NUMDATOSMAX 5000
#define NUMRESULTMAX 5000
#define NUMSEASONSMAX 12

// Define max() and min() functions
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

// Function prototypes
void spi(float dataSeries[], int n, int seasons, float spiSeries[]);

// Main program:
// Calculate the Standardized Precipitation Index
//int main(int argc, char **argv) {
int main(int argc, char **argv) {

	FILE *entrada,*salida;
	char pathOrigen[30],pathDestino[30],estacion[36],latitud[6];
	float rainSeries[NUMDATOSMAX],acumSeries[NUMDATOSMAX],
		spiSeries[NUMDATOSMAX];
	int anio,mes,seasonality,acumulated,numRegistros,acumRegistros,
		indice,jndice;

	// Initialize variables
	anio=mes=seasonality=acumulated=numRegistros=acumRegistros=indice=jndice=0;
	for (indice=0; indice<NUMDATOSMAX; indice++) rainSeries[indice]=spiSeries[indice]=0.0;

	// Read in-line arguments
	if (argc!=4) {
		printf("\nUsage:\tspi [acumulated] [source file] [result file]\n");
		exit(1);
	}
	sscanf(argv[1], "%u", &acumulated);
	sscanf(argv[2], "%s", pathOrigen);
	sscanf(argv[3], "%s", pathDestino);

	// Open input file
	if((entrada=fopen(pathOrigen,"rt"))==NULL)
		{
		printf("\nError: File can't be opened");
		exit(1);
		}
	// Read heading
	fgets (estacion, 36, entrada);
	fscanf(entrada, "%u;%u\n", &anio,&mes);
	fscanf(entrada, "%u\n", &seasonality);
	if(seasonality>NUMSEASONSMAX) {
		printf("\nError: Too many seasons. Maximum is %d", NUMSEASONSMAX);
		exit(1);
	}
	// Read data
	indice=0;
	while(!feof(entrada)) {
		if(indice==NUMDATOSMAX) {
			printf("\nError: Too many data in input file. Maximum is %d", NUMDATOSMAX);
			exit(1);
		}
		fscanf(entrada,"%f\n", &rainSeries[indice]);
		indice++;
	}
	numRegistros=indice;
	// Close file
	fclose(entrada);
	// Print metadata (just to check)
	printf("\nseries: %s", estacion);
	printf("initial date: %d/%d\n", mes, anio);
	printf("seasonality: %d\n", seasonality);
	printf("%d registers\n", numRegistros);
	printf("calculating SPI at %d month", acumulated);
	if (acumulated>1) printf("s");
	printf("\n");

	// Compute the cumulative series
	anio += (acumulated-1)/12;
	mes += acumulated-1;
	while (mes>12) mes-=12;
	acumRegistros = numRegistros-acumulated+1;
	for (indice=acumulated-1; indice<numRegistros; indice++) {
		for (jndice=0; jndice<acumulated; jndice++) {
			acumSeries[indice-acumulated+1] += rainSeries[indice-jndice];
		}
	}

	// Compute the SPI series
	spi(acumSeries, acumRegistros, seasonality, spiSeries);

	// Write results to file
	if((salida=fopen(pathDestino,"wt"))==NULL) {
		printf("\nError: Output file could not be opened");
		exit(1);
	}
	fprintf(salida,"%s%u;%u\n%u", estacion,anio,mes,seasonality);
	//for (jndice=1; jndice<=seasonality; jndice++) {
	//	fprintf(salida,"\n%f;%f;%f", logLogisticParams[jndice][0],
	//			logLogisticParams[jndice][1], logLogisticParams[jndice][2]);
	//}
	for (indice=0; indice<acumRegistros; indice++) {
		fprintf(salida,"\n%f", spiSeries[indice]);
	}
	fclose(salida);

	// Quit
	exit(0);
}

// spi()
// Calculates the Standardized Precipitation Index from a series. The
// SPI is the standardized value of the cumulative precipitation over
// a period defined by 'acumulated', computed following a Pearson III
// probability distribution.
void spi(float dataSeries[], int n, int seasons, float spiSeries[]) {

	int i, j, k, nSeason;
	float seasonSeries[NUMDATOSMAX], lMoment[4], gammaParams[NUMSEASONSMAX+1][2],
		  pearsonIIIParams[NUMSEASONSMAX+1][3];

	// Loop through all seasons defined by seasons
	for (j=1; j<=seasons; j++) {
		// Extract and sort the seasonal series
		k = 0;
		for (i=j-1; i<n; i+=seasons) {
			seasonSeries[k] = dataSeries[i];
			k++;
		}
		nSeason = k;
		upward(seasonSeries, nSeason);
		// Compute L-moments
		lMoments(seasonSeries, nSeason, lMoment, -0.35, 0);
		//lMoments(seasonSeries, nSeason, lMoment, 0, 0);
		//printf("\n%u season", j);
		//printf("\nL-moment order 1: %.4f", lMoment[1]);
		//printf("\nL-moment order 2: %.4f", lMoment[2]);
		// Fit a Gamma distribution
		//gammaFit(lMoment, gammaParams[j]);
		//printf("\nGamma distribution alfa param.: %.4f\n", gammaParams[j][0]);
		//printf("Gamma distribution beta param.: %.4f\n", gammaParams[j][1]);
		// Calculate the standardized values following a Gamma distribution
		//for (i=j-1; i<n; i+=seasons) {
		//	spiSeries[i] = gammaStandardize(dataSeries[i], gammaParams[j]);
		//}
		// Fit a Pearson III distribution
		pearsonIIIFit(lMoment, pearsonIIIParams[j]);
		//printf("PIII distribution origin param.: %.4f\n", pearsonIIIParams[j][0]);
		//printf("PIII distribution alpha param.: %.4f\n", pearsonIIIParams[j][1]);
		//printf("PIII distribution beta param.: %.4f\n", pearsonIIIParams[j][2]);
		// Calculate the standardized values following a Pearson III distribution
		for (i=j-1; i<n; i+=seasons) {
			spiSeries[i] = pearsonIIIStandardize(dataSeries[i], pearsonIIIParams[j]);
		}
	}
}
