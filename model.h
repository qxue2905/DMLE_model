#ifndef MODEL_H
#define MODEL_H

#include "soil_params.h"
#include "soil_layers.h"

typedef struct
{
    double Ta;          
    double qa;          
    double U;           
    double P;           
    double E_eff;       
    double obs_LE;      
    double obs_H;       
    double rain;        
    long long datetime; 
} DataPoint;

typedef struct
{
    DataPoint *points; 
    int n;             
    int capacity;      
} Dataset;

typedef struct
{
    double depth; 
    double psi; 
    double coupling_factor;
    int is_active; 
    double wtd_max_depth; 
    double h_flx;         
} WaterTable;

typedef struct
{
    double Ta;     
    double qa;     
    double U;      
    double P;      
    double E_eff;  
    double obs_LE; 
    double obs_H;  
    double rain;   
} Forcing;

typedef struct
{
    double LE;              
    double Delta;           
    double rho_a;         
    double Cp;              
    double Lv;              
    double Ra;              
    double qs;              
    double hu;              
    double theta1;         
    double psi;           
    double Rs;              
    double ka;              
    double film_thickness;  
    double De;              
    double n_theta;         
    double theta_m;         
	char   ka_formula[512];   
} FluxResult;

typedef struct
{
    double *LE_sim;             
    double *Delta_sim;          
    double *rho_a_sim;          
    double *Cp_sim;            
    double *Lv_sim;             
    double *Ra_sim;            
    double *qs_sim;           
    double *hu_sim;             
    double *theta1_sim;         
    double *psi_sim;            
    double *Rs_sim;            
    double *ka_sim;           
    double *film_thickness_sim; 
    double *De_sim;             
    double *n_theta_sim;        
    double *theta_m_sim;        
	char (*ka_formula_sim)[512];
    int n;                     
} ModelResult;

double esat(double T);
double qsat(double T, double P);
double Delta_func(double T, double qsat);
double calc_hu(double psi, double Ta);
double calc_latent_heat(double Ta);
double calc_air_density(double P, double Ta, double qa);
double calc_aerodynamic_resistance(double U, double h_flx);
double expint_E1(double x);
double Ei_neg(double x);
double Pi_e(double film_thickness, double T);
double Pi_m(double film_thickness);
double Pi_total(double film_thickness, double T);
double solve_film_thickness(double hm, double T);
double compute_film_thickness(double h_m);
double compute_B(double film_thickness, double T);
double compute_De(double theta_m, double n, double hm_m);
double compute_ka(double film_thickness, double De, double n, double T);
double compute_ka_with_formula(double film_thickness, double De, double n, double T,
                               char *formula_out, int formula_size);
void init_water_table(WaterTable *wt, double depth, double coupling_factor, int is_active);
double get_bottom_node_depth(const SoilLayers *layers);
double calculate_water_table_flux(const SoilParams *soil_params, SoilState *state,
                                  const SoilLayers *layers, const WaterTable *wt,
                                  double psi_last, double theta_last, double dt);
void compute_hydraulic_params(const SoilParams *soil_params, SoilState *state,
                              const SoilLayers *layers, double *psi, double *K);
double calc_infiltration_capacity(const SoilParams *soil_params, double theta);
void calculate_baresoil_evaporation(const SoilParams *soil_params, SoilState *state,
                                    const SoilLayers *layers, const Forcing *forcing,
                                    const WaterTable *wt, double dt, FluxResult *result);
void update_soil_moisture_wt(const SoilParams *soil_params, SoilState *state,
                             const SoilLayers *layers, double Pg, double Eg,
                             const WaterTable *wt, double dt,
                             double *horton_out);
int run_model_with_params(const Dataset *ds, double dt, int n_layers,
                          const WaterTable *wt, const SoilParams *soil_params,
                          ModelResult **result_out);
void free_model_result(ModelResult *result);

#endif /* MODEL_H */