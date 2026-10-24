#include "soil_params.h"
#include <stdio.h>
#include <math.h>

double calc_theta_sat(double PSAND, double PCLAY)
{
    double theta_sat = 0.505 - 0.142 * PSAND - 0.037 * PCLAY;

    if (theta_sat < 0.01)
        theta_sat = 0.01;
    if (theta_sat > 0.90)
        theta_sat = 0.90;

    return theta_sat;
}

double calc_bcoef(double PSAND, double PCLAY)
{
    double bcoef = 3.10 + 15.7 * PCLAY - 0.3 * PSAND;

    if (bcoef < 1.0)
        bcoef = 1.0;
    if (bcoef > 20.0)
        bcoef = 20.0;

    return bcoef;
}

double calc_psi_sat(double PSAND, double PCLAY)
{
    double PSILT = 1.0 - PSAND - PCLAY;
    double psi_sat = -0.01 * pow(10.0, 1.54 - 0.95 * PSAND + 0.63 * PSILT);

    if (psi_sat > -0.001)
        psi_sat = -0.001;
    if (psi_sat < -10.0)
        psi_sat = -10.0;

    return psi_sat;
}

double calc_k_sat(double PSAND, double PCLAY)
{
    double k_sat = 0.0254 * pow(10.0, -0.6 + 1.26 * PSAND - 0.64 * PCLAY) / 3600.0;
    if (k_sat < 1e-12)
        k_sat = 1e-12;
    if (k_sat > 1e-3)
        k_sat = 1e-3;
    return k_sat;
}

void calc_soil_params_from_texture(double PSAND, double PCLAY, SoilParams *params)
{
    if (PSAND < 0 || PSAND > 1 || PCLAY < 0 || PCLAY > 1)
    {
        fprintf(stderr, "Warning: PSAND and PCLAY should be in 0-1 range\n");
    }

    params->theta_sat = calc_theta_sat(PSAND, PCLAY);
    params->b = calc_bcoef(PSAND, PCLAY);
    params->psi_sat = calc_psi_sat(PSAND, PCLAY);
    params->k_sat = calc_k_sat(PSAND, PCLAY);
	
    double PSILT = 1.0 - PSAND - PCLAY;   
    if (PSILT < 0.0) PSILT = 0.0;         
    params->Rs_unit = 0.25 * PSAND + 1.3 * PSILT + 3.8 * PCLAY;
}