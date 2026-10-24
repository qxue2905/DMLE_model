#include "model.h"
#include "solver.h"
#include "constants.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

double esat(double T)
{
    return 611.0 * exp(17.27 * T / (T + 237.3));
}

double qsat(double T, double P)
{
    double e_s = esat(T);
    return EPSILON * e_s / (P - 0.378 * e_s);
}

double Delta_func(double T, double qsat)
{
    return (17.502 * 240.97 * qsat) / ((240.97 + T) * (240.97 + T));
}

double calc_hu(double psi, double Ta)
{
    double Tk = Ta + T0;
    double hu = exp(psi * G / (RV * Tk));
    return hu;
}

double calc_latent_heat(double Ta)
{
    return LAMBDA_0 - LAMBDA_COEF * Ta;
}

double calc_air_density(double P, double Ta, double qa)
{
    double R_moist = R_DRY * (1 + 0.608 * qa);
    return P / (R_moist * (Ta + T0));
}

double calc_aerodynamic_resistance(double U, double h_flx)
{
    double U_eff = U;
    if (U_eff < U_MIN)
        U_eff = U_MIN;

    double h = h_flx;

    double Ra = (1.0 / (VON_KARM * VON_KARM)) *
                log(h / Z0_BARESOIL) * log(h / Z0_BARESOIL) / U_eff;

    if (Ra < 1.0)
        Ra = 1.0;

    return Ra;
}

double expint_E1(double x)
{
    const double gamma = 0.5772156649015328606;  
    const double EPS   = 1e-15;
    const double FPMIN = 1e-300;

    if (x <= 0.0) return 1e300;   

    if (x < 1.0) {
        double sum  = -gamma - log(x);
        double term = 1.0;                
        for (int k = 1; k <= 200; k++) {
            term *= -x / (double)k;
            double add = -term / (double)k; 
            sum += add;
            if (fabs(add) < fabs(sum) * EPS) break;
        }
        return sum;
    } else {
        double b = x + 1.0;
        double c = 1.0 / FPMIN;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 200; i++) {
            double a = -(double)i * (double)i;
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            double del = c * d;
            h *= del;
            if (fabs(del - 1.0) < EPS) break;
        }
        return exp(-x) * h;
    }
}

double Ei_neg(double x)
{
    return -expint_E1(x);
}

double Pi_e(double film_thickness, double T)
{
    double inner = PI * KB * T / (E_CHARGE * Z_VALENCE);
    return 0.5 * EPS_R * EPS_0 * inner * inner
           / (film_thickness * film_thickness);
}

double Pi_m(double film_thickness)
{
    return -A_SVL / (6.0 * PI
                     * film_thickness * film_thickness * film_thickness);
}

double Pi_total(double film_thickness, double T)
{
    return Pi_e(film_thickness, T) + Pi_m(film_thickness);
}

double solve_film_thickness(double hm, double T)
{
    double target = -RHO_O * G * hm;   

    double lo = 1.0e-12;  
    double hi = 1.0e-5;    

    double flo = Pi_total(lo, T) - target;
    double fhi = Pi_total(hi, T) - target;

    if (flo * fhi > 0.0)
        return (fabs(flo) < fabs(fhi)) ? lo : hi;

    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        double fm  = Pi_total(mid, T) - target;
        if (fabs(fm) < 1e-12 * target + 1e-30) return mid;
        if (flo * fm <= 0.0) { hi = mid; fhi = fm; }
        else                 { lo = mid; flo = fm; }
    }
    return 0.5 * (lo + hi);
}

double compute_film_thickness(double h_m)
{
    if (h_m >= 0.0) return 0.0;
    double numerator = -A_SVL;  
    double denominator = 6.0 * PI * RHO_O * G * h_m;
    double f_cubed = -numerator / denominator;
    if (f_cubed <= 0.0) return 0.0;
    return pow(f_cubed, 1.0 / 3.0);
}

double compute_B(double film_thickness, double T)
{
    double a = 1.621e-7 / T;          
    double d = film_thickness;       
    double x = a / d;                 

    double E1x = expint_E1(x);
    double bracket = (4.0 - 5.0*x - x*x) * exp(-x)
                   + (6.0*x*x + x*x*x) * E1x;
    double B = d * d * d * bracket;
    if (B < 0.0) B = 0.0;

    return B;
}

double compute_De(double theta_m, double n, double hm_m)
{
    double coef     = A_SVL / (6.0 * PI * RHO_O * G * hm_m);
    double abs_coef = fabs(coef);
    double inner    = pow(abs_coef, -1.0/3.0);
    double denom    = (theta_m / (1.0 - n)) * inner;
    return 6.0 / denom;
}

double compute_ka(double film_thickness, double De, double n, double T)
{
    double common = RHO_O * G / (PI * MU_O * De) * (1.0 - n);

    double ka;
    if (film_thickness >= FILM_CUT) {
        ka = 4.0 * common
             * film_thickness * film_thickness * film_thickness;
    } else {
        ka = common * compute_B(film_thickness, T);
    }

    if (ka < 0.0) ka = 0.0;

    return ka;
}

double compute_ka_with_formula(double film_thickness, double De, double n, double T,
                               char *formula_out, int formula_size)
{
    double common = RHO_O * G / (PI * MU_O * De) * (1.0 - n);

    double ka;
    if (film_thickness >= FILM_CUT) {
        double d3 = film_thickness * film_thickness * film_thickness;
        ka = 4.0 * common * d3;

        if (formula_out && formula_size > 0) {
            snprintf(formula_out, formula_size,
                     "ka=4*rho*g/(pi*mu*De)*(1-n)*d^3; "
                     "rho=%.2f; g=%.5f; mu=%.2e; De=%.3e; n=%.4f; d=%.3e; "
                     "4*%.2f*%.5f/(pi*%.2e*%.3e)*(1-%.4f)*(%.3e)^3=%.6e",
                     RHO_O, G, MU_O, De, n, film_thickness,
                     RHO_O, G, MU_O, De, n, film_thickness, ka);
        }
    } else {
        double a = 1.621e-7 / T;
        double x = a / film_thickness;
        double E1x = expint_E1(x);
        double bracket = (4.0 - 5.0*x - x*x) * exp(-x)
                       + (6.0*x*x + x*x*x) * E1x;
        double B = film_thickness * film_thickness * film_thickness * bracket;
        if (B < 0.0) B = 0.0;
        ka = common * B;

        if (formula_out && formula_size > 0) {
            snprintf(formula_out, formula_size,
                     "ka=rho*g/(pi*mu*De)*(1-n)*B; "
                     "rho=%.2f; g=%.5f; mu=%.2e; De=%.3e; n=%.4f; d=%.3e; "
                     "a=1.621e-7/T=%.3e; x=a/d=%.4f; B=%.6e; "
                     "%.2f*%.5f/(pi*%.2e*%.3e)*(1-%.4f)*%.6e=%.6e",
                     RHO_O, G, MU_O, De, n, film_thickness,
                     a, x, B,
                     RHO_O, G, MU_O, De, n, B, ka);
        }
    }

    if (ka < 0.0) ka = 0.0;
    return ka;
}

void init_water_table(WaterTable *wt, double depth, double coupling_factor, int is_active)
{
    wt->depth = depth;
    wt->coupling_factor = coupling_factor;
    wt->is_active = is_active;
    wt->wtd_max_depth = WTD_MAX_DEPTH;

    if (wt->depth <= 0)
    {
        wt->depth = 0.001;
    }

    if (wt->coupling_factor < 0)
        wt->coupling_factor = 0;
    if (wt->coupling_factor > 1)
        wt->coupling_factor = 1;
}

double get_bottom_node_depth(const SoilLayers *layers)
{
    int n = layers->n_layers;

    if (n == 1)
    {
        return layers->dz[0] / 2.0;
    }
    else
    {
        return (layers->dg[n - 1] + layers->dg[n - 2]) / 2.0;
    }
}

double calculate_water_table_flux(const SoilParams *soil_params, SoilState *state,
                                  const SoilLayers *layers, const WaterTable *wt,
                                  double psi_last, double theta_last, double dt)
{
    int n = layers->n_layers;
    int last_layer = n - 1;

    if (!wt->is_active)
    {
        return 0.0;
    }

    double theta_ice = state->theta_ice[last_layer];
    double w_sat_eff = soil_params[last_layer].theta_sat - theta_ice;
    if (w_sat_eff < THETA_MIN)
        w_sat_eff = THETA_MIN;

    double s = theta_last / w_sat_eff;
    if (s > 1)
        s = 1;

    double k_ratio = theta_last / w_sat_eff;
    if (k_ratio > 1)
        k_ratio = 1;

    double frozen_frac = theta_ice / (theta_last + theta_ice + 1e-10);

    double k_last = soil_params[last_layer].k_sat * pow(k_ratio, 2 * soil_params[last_layer].b + 3);
    k_last *= exp(-6 * frozen_frac);

    if (wt->depth > wt->wtd_max_depth)
    {
        return -k_last;
    }

    double wtd_depth = wt->depth;

    if (wtd_depth < layers->dg[last_layer])
    {
        wtd_depth = layers->dg[last_layer];
    }

    double z_node = get_bottom_node_depth(layers);
    double s_wtd = 1.0;
    double psi_wtd = soil_params[last_layer].psi_sat * pow(s_wtd, -soil_params[last_layer].b);
    double dz = wtd_depth - z_node;
    double head = (psi_last - psi_wtd) / dz;
    double flux = -k_last * wt->coupling_factor * head - k_last;

    return flux;
}

void compute_hydraulic_params(const SoilParams *soil_params, SoilState *state,
                              const SoilLayers *layers, double *psi, double *K)
{
    int n = layers->n_layers;

    for (int i = 0; i < n; i++)
    {
        double theta_total = state->theta_liq[i] + state->theta_ice[i];
        double frozen_frac = state->theta_ice[i] / (theta_total + 1e-10);

        double w_sat_eff = soil_params[i].theta_sat - state->theta_ice[i];
        if (w_sat_eff < THETA_MIN)
            w_sat_eff = THETA_MIN;

        double s = state->theta_liq[i] / w_sat_eff;
        if (s > 1)
            s = 1;

        psi[i] = soil_params[i].psi_sat * pow(s, -soil_params[i].b);

        double k_ratio = state->theta_liq[i] / w_sat_eff;
        if (k_ratio > 1)
            k_ratio = 1;

        K[i] = soil_params[i].k_sat * pow(k_ratio, 2 * soil_params[i].b + 3);
        K[i] *= exp(-6 * frozen_frac);
    }
}

double calc_infiltration_capacity(const SoilParams *soil_params, double theta)
{
    double theta_sat = soil_params[0].theta_sat;
    double b = soil_params[0].b;
    double psi_sat = soil_params[0].psi_sat;
    double k_sat = soil_params[0].k_sat;
    double dz = 0.01;

    double I_max = k_sat * (b * psi_sat * (theta / theta_sat - 1.0) / dz + 1.0);

    if (I_max < 0.0) I_max = 0.0;

    return I_max;
}

void calculate_baresoil_evaporation(const SoilParams *soil_params, SoilState *state,
                                    const SoilLayers *layers, const Forcing *forcing,
                                    const WaterTable *wt, double dt, FluxResult *result)
{
    double Ta = forcing->Ta;
    double qa = forcing->qa;
    double U = forcing->U;
    double P = forcing->P;
    double E_eff = forcing->E_eff;
    double obs_LE = forcing->obs_LE;
    double obs_H = forcing->obs_H;
    double theta1 = state->theta_liq[0];
    double theta_ice1 = state->theta_ice[0];

    if (isnan(Ta) || isnan(qa) || isnan(U) || isnan(P) ||
        isnan(E_eff) || isnan(obs_LE) || isnan(obs_H))
    {
        result->LE = 0.0;
        result->Delta = 0.0;
        result->rho_a = 1.2;
        result->Cp = 1005.0;
        result->Lv = 2.45e6;
        result->Ra = 50.0;
        result->qs = 0.0;
        result->hu = 0.5;
        result->theta1 = state->theta_liq[0];
        result->psi = -0.5;
        result->Rs = 0.0;
        result->ka = 0.0;
        result->film_thickness = 0.0;
        result->De = 0.0;
        result->n_theta = 0.0;
        result->theta_m = 0.0;
		result->ka_formula[0] = '\0';
        return;
    }

    double Lv = calc_latent_heat(Ta);
    double R_moist = R_DRY * (1 + 0.608 * qa);
    double Tk = Ta + T0;
    double rho_a = P / (R_moist * Tk);
    double e_s = esat(Ta);
    double q_s = EPSILON * e_s / (P - 0.378 * e_s);
    double Delta = (17.502 * 240.97 * q_s) / ((240.97 + Ta) * (240.97 + Ta));
    double Cp = CP_DRY + CP_VAPOR * qa;
    double Ra = calc_aerodynamic_resistance(U, wt->h_flx);
    double w_sat_eff = soil_params[0].theta_sat - theta_ice1;
    if (w_sat_eff < THETA_MIN)
        w_sat_eff = THETA_MIN;
    double s = theta1 / w_sat_eff;
    if (s > 1)
        s = 1;

    double psi = soil_params[0].psi_sat * pow(s, -soil_params[0].b);
    double hu = exp(psi * G / (RV * Tk));
    double hm_m = -1.0e3;                              
    double theta_m = pow(hm_m / soil_params[0].psi_sat, -1 / soil_params[0].b)
                     * soil_params[0].theta_sat;
    double n_theta = soil_params[0].theta_sat;        

    double De = compute_De(theta_m, n_theta, hm_m);
    double film_thickness = compute_film_thickness(psi);
	char ka_formula[512] = {0};
	double ka = compute_ka_with_formula(film_thickness, De, n_theta, Tk,
                                    ka_formula, sizeof(ka_formula));
	
	double Rs_unit_full = layers->dz[0] * 1000 * 1000 * soil_params[0].Rs_unit * theta1;
	double Rs = Rs_unit_full * rho_a / RHO_O * 2.56 * 0.00001 / ka * q_s * hu * 9.8 / 461.524 / Tk;

	double term1 = Delta * E_eff;
    double term2 = rho_a * Cp / Ra * (q_s - qa / hu);
    double numerator = term1 + term2;
	double Ra_cor = Ra * (1 + 10 / log(wt->h_flx / Z0_BARESOIL));
    double term3 = Cp / Lv / hu * ((Ra_cor + Rs) / Ra);
    double denominator = Delta + term3;          
    double LE = 0.0;
    if (!isnan(denominator) && fabs(denominator) > 1e-10)
    {
        LE = numerator / denominator;
    }

    result->LE = LE;
    result->Delta = Delta;
    result->rho_a = rho_a;
    result->Cp = Cp;
    result->Lv = Lv;
    result->Ra = Ra;
    result->qs = q_s;
    result->hu = hu;
    result->theta1 = theta1;
    result->psi = psi;
    result->Rs = Rs;
    result->ka = ka;
    result->film_thickness = film_thickness;
    result->De = De;
    result->n_theta = n_theta;
    result->theta_m = theta_m;
	snprintf(result->ka_formula, sizeof(result->ka_formula),
	         "%s", ka_formula);
}

void update_soil_moisture_wt(const SoilParams *soil_params, SoilState *state,
                             const SoilLayers *layers, double Pg, double Eg,
                             const WaterTable *wt, double dt,
                             double *horton_out)
{
    int n = layers->n_layers;
    if (n <= 0) return;

    if (horton_out) *horton_out = 0.0;

    double I_max = calc_infiltration_capacity(soil_params, state->theta_liq[0]);

    double horton = 0.0;
    double Pg_effective = Pg;

    if (Pg > I_max && I_max > 0.0) {
        horton = Pg - I_max;
        Pg_effective = I_max;
    }

    if (horton_out) *horton_out = horton;

    if (Pg_effective < 0.0) Pg_effective = 0.0;

    double *psi         = (double *)malloc(n * sizeof(double));
    double *K           = (double *)malloc(n * sizeof(double));
    double *flux        = (double *)malloc(n * sizeof(double));

    if (!psi || !K || !flux) {
        free(psi); free(K); free(flux);
        return;
    }

    compute_hydraulic_params(soil_params, state, layers, psi, K);

    for (int i = 0; i < n - 1; i++) {
        double K_int = sqrt(K[i] * K[i + 1]);
        double head  = (psi[i] - psi[i + 1]) / layers->dzdif[i];
        flux[i] = -K_int * head - K_int;
    }

    int last_layer = n - 1;

    if (wt->is_active) {
        flux[last_layer] = calculate_water_table_flux(soil_params, state, layers, wt,
                                                     psi[last_layer],
                                                     state->theta_liq[last_layer], dt);
    } else {
        flux[last_layer] = -K[last_layer];
    }
	
    double *dF_dtheta = (double *)calloc(n * n, sizeof(double));

    for (int i = 0; i < n - 1; i++) {
        double K_int = sqrt(K[i] * K[i + 1]);
        double head  = (psi[i] - psi[i + 1]) / layers->dzdif[i];

        double theta_i = state->theta_liq[i];
        if (theta_i < THETA_MIN) theta_i = THETA_MIN;
        double theta_next = state->theta_liq[i + 1];
        if (theta_next < THETA_MIN) theta_next = THETA_MIN;

        double dKint_dtheta_i    = (2 * soil_params[i].b + 3) / 2 * K_int / theta_i;
        double dKint_dtheta_next = (2 * soil_params[i + 1].b + 3) / 2 * K_int / theta_next;

        double dpsi_dtheta_i    = -soil_params[i].b * psi[i] / theta_i;
        double dpsi_dtheta_next = -soil_params[i + 1].b * psi[i + 1] / theta_next;

        double dhead_dtheta_i    =  dpsi_dtheta_i / layers->dzdif[i];
        double dhead_dtheta_next = -dpsi_dtheta_next / layers->dzdif[i];

        dF_dtheta[i * n + i]     = -(dKint_dtheta_i    * head + K_int * dhead_dtheta_i)    - dKint_dtheta_i;
        dF_dtheta[i * n + i + 1] = -(dKint_dtheta_next * head + K_int * dhead_dtheta_next) - dKint_dtheta_next;
    }

    double theta_n = state->theta_liq[last_layer];
    if (theta_n < THETA_MIN) theta_n = THETA_MIN;

    double dK_dtheta_n   = (2 * soil_params[last_layer].b + 3) * K[last_layer] / theta_n;
    double dpsi_dtheta_n = -soil_params[last_layer].b * psi[last_layer] / theta_n;

    double wtd_depth = wt->depth;
    if (wtd_depth < layers->dg[last_layer]) wtd_depth = layers->dg[last_layer];
    double z_node = get_bottom_node_depth(layers);
    double dz = wtd_depth - z_node;
    if (fabs(dz) < 1e-6) dz = 1e-6;

    double psi_wtd = soil_params[last_layer].psi_sat;
    double head_bottom = (psi[last_layer] - psi_wtd) / dz;

    if (wt->is_active) {
        double dF_dtheta_bottom = -dK_dtheta_n * wt->coupling_factor * head_bottom
                                  - K[last_layer] * wt->coupling_factor * dpsi_dtheta_n / dz
                                  - dK_dtheta_n;
        dF_dtheta[last_layer * n + last_layer] = dF_dtheta_bottom;
    } else {
        dF_dtheta[last_layer * n + last_layer] = -dK_dtheta_n;
    }

    double theta_w = THETA_W;
    double *a            = (double *)calloc(n, sizeof(double));
    double *b            = (double *)calloc(n, sizeof(double));
    double *c            = (double *)calloc(n, sizeof(double));
    double *d            = (double *)calloc(n, sizeof(double));
    double *delta_theta  = (double *)calloc(n, sizeof(double));

    if (!a || !b || !c || !d || !delta_theta) {
        free(psi); free(K); free(flux); free(dF_dtheta);
        free(a); free(b); free(c); free(d); free(delta_theta);
        return;
    }

    a[0] = 0.0;
    b[0] = 1.0 - theta_w * dt / layers->dz[0] * dF_dtheta[0 * n + 0];
    c[0] = -theta_w * dt / layers->dz[0] * dF_dtheta[0 * n + 1];
    d[0] = dt / layers->dz[0] * (flux[0] + Pg_effective - Eg);

    for (int i = 1; i < n - 1; i++) {
        double dzi = layers->dz[i];
        double dF_net_dtheta_im1 = dF_dtheta[i * n + i - 1] - dF_dtheta[(i - 1) * n + i - 1];
        double dF_net_dtheta_i   = dF_dtheta[i * n + i]     - dF_dtheta[(i - 1) * n + i];
        double dF_net_dtheta_ip1 = dF_dtheta[i * n + i + 1];

        a[i] = -theta_w * dt / dzi * dF_net_dtheta_im1;
        b[i] = 1.0 - theta_w * dt / dzi * dF_net_dtheta_i;
        c[i] = -theta_w * dt / dzi * dF_net_dtheta_ip1;
        d[i] = dt / dzi * (flux[i] - flux[i - 1]);
    }

    int i_last = n - 1;
    double dzi_last = layers->dz[i_last];
    double dF_net_dtheta_im1 = dF_dtheta[i_last * n + i_last - 1] - dF_dtheta[(i_last - 1) * n + i_last - 1];
    double dF_net_dtheta_i   = dF_dtheta[i_last * n + i_last]     - dF_dtheta[(i_last - 1) * n + i_last];

    a[i_last] = -theta_w * dt / dzi_last * dF_net_dtheta_im1;
    b[i_last] = 1.0 - theta_w * dt / dzi_last * dF_net_dtheta_i;
    c[i_last] = 0.0;
    d[i_last] = dt / dzi_last * (flux[i_last] - flux[i_last - 1]);

    int success = solve_tridiagonal(a, b, c, d, delta_theta, n);

    if (success != 0) {
        delta_theta[0] = dt / layers->dz[0] * (flux[0] + Pg_effective - Eg);
        for (int i = 1; i < n; i++) {
            delta_theta[i] = dt / layers->dz[i] * (flux[i] - flux[i - 1]);
        }
    }

    for (int i = 0; i < n; i++) {
        state->theta_liq[i] += delta_theta[i];
    }

    for (int i = 0; i < n - 1; i++) {
        double w_sat_eff_i = soil_params[i].theta_sat - state->theta_ice[i];
        if (w_sat_eff_i < THETA_MIN) w_sat_eff_i = THETA_MIN;

        if (state->theta_liq[i] > w_sat_eff_i) {
            double excess = (state->theta_liq[i] - w_sat_eff_i) * layers->dz[i] / dt;
            state->theta_liq[i] = w_sat_eff_i;
            state->theta_liq[i + 1] += excess * dt / layers->dz[i + 1];
        }
    }

    {
        double w_sat_eff_last = soil_params[last_layer].theta_sat - state->theta_ice[last_layer];
        if (w_sat_eff_last < THETA_MIN) w_sat_eff_last = THETA_MIN;
        if (state->theta_liq[last_layer] > w_sat_eff_last) {
            state->theta_liq[last_layer] = w_sat_eff_last;
        }
    }

    for (int i = 0; i < n - 1; i++) {
        if (state->theta_liq[i] < THETA_MIN) {
            double deficit = (THETA_MIN - state->theta_liq[i]) * layers->dz[i] / dt;
            state->theta_liq[i] = THETA_MIN;
            state->theta_liq[i + 1] -= deficit * dt / layers->dz[i + 1];
        }
    }

    if (state->theta_liq[last_layer] < THETA_MIN) {
        state->theta_liq[last_layer] = THETA_MIN;
    }

    if (wt->is_active) {
        for (int i = 0; i < n; i++) {
            double layer_bottom = layers->dg[i];
            if (wt->depth < layer_bottom) {
                double w_sat_eff = soil_params[i].theta_sat - state->theta_ice[i];
                if (w_sat_eff < THETA_MIN) w_sat_eff = THETA_MIN;
                state->theta_liq[i] = w_sat_eff;
            }
        }
    }

    free(psi);
    free(K);
    free(flux);
    free(dF_dtheta);
    free(a);
    free(b);
    free(c);
    free(d);
    free(delta_theta);
}