#ifndef SOIL_PARAMS_H
#define SOIL_PARAMS_H

#include "constants.h"

typedef struct
{
    double theta_sat; 
    double b; 
    double psi_sat; 
    double k_sat; 	
	double Rs_unit; 
} SoilParams;

double calc_theta_sat(double PSAND, double PCLAY);
double calc_bcoef(double PSAND, double PCLAY);
double calc_psi_sat(double PSAND, double PCLAY);
double calc_k_sat(double PSAND, double PCLAY);
void calc_soil_params_from_texture(double PSAND, double PCLAY, SoilParams *params);

#endif /* SOIL_PARAMS_H */