#ifndef SOIL_LAYERS_H
#define SOIL_LAYERS_H

#define SOIL_LAYER_DEPTHS {0.01, 0.04, 0.10, 0.20, 0.40, 0.60, 0.80, \
                           1.00, 1.50, 2.00, 3.00, 5.00, 8.00, 12.0} /* 累计深度 */
#define SOIL_N_LAYERS 14

typedef struct
{
    int n_layers;  
    double *dz;    
    double *dg;    
    double *dzdif; 
} SoilLayers;

typedef struct
{
    int n_layers;
    double *theta_liq; 
    double *theta_ice; 
} SoilState;

SoilLayers *create_soil_layers(int n_layers);
void free_soil_layers(SoilLayers *layers);
SoilState *create_soil_state(int n_layers);
void free_soil_state(SoilState *state);
void init_soil_state(SoilState *state, double init_theta);
void copy_soil_state(SoilState *dest, const SoilState *src, int n_layers);

#endif /* SOIL_LAYERS_H */