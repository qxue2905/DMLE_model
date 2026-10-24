#include "soil_layers.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static const double ISBA_DEPTHS[SOIL_N_LAYERS] = {
    0.01, 0.04, 0.10, 0.20, 0.40, 0.60, 0.80,
    1.00, 1.50, 2.00, 3.00, 5.00, 8.00, 12.0};

SoilLayers *create_soil_layers(int n_layers)
{
    if (n_layers < 1 || n_layers > 14)
    {
        n_layers = 14; 
    }

    SoilLayers *layers = (SoilLayers *)malloc(sizeof(SoilLayers));
    if (!layers)
        return NULL;

    layers->n_layers = n_layers;
    layers->dz = (double *)calloc(n_layers, sizeof(double));
    layers->dg = (double *)calloc(n_layers, sizeof(double));
    layers->dzdif = (double *)calloc(n_layers, sizeof(double));

    if (!layers->dz || !layers->dg || !layers->dzdif)
    {
        free(layers->dz);
        free(layers->dg);
        free(layers->dzdif);
        free(layers);
        return NULL;
    }

    for (int i = 0; i < n_layers; i++)
    {
        layers->dg[i] = ISBA_DEPTHS[i];
    }

    layers->dz[0] = layers->dg[0];
    for (int i = 1; i < n_layers; i++)
    {
        layers->dz[i] = layers->dg[i] - layers->dg[i - 1];
    }

    if (n_layers > 1)
    {
        double *z_node = (double *)malloc(n_layers * sizeof(double));

        z_node[0] = layers->dg[0] / 2.0; 
        for (int i = 1; i < n_layers; i++)
        {
            z_node[i] = (layers->dg[i - 1] + layers->dg[i]) / 2.0; 
        }

        for (int i = 0; i < n_layers - 1; i++)
        {
            layers->dzdif[i] = z_node[i + 1] - z_node[i];
        }

        layers->dzdif[n_layers - 1] = layers->dz[n_layers - 1];

        free(z_node);
    }

    return layers;
}

void free_soil_layers(SoilLayers *layers)
{
    if (layers)
    {
        free(layers->dz);
        free(layers->dg);
        free(layers->dzdif);
        free(layers);
    }
}

SoilState *create_soil_state(int n_layers)
{
    SoilState *state = (SoilState *)malloc(sizeof(SoilState));
    if (!state)
        return NULL;

    state->n_layers = n_layers;
    state->theta_liq = (double *)calloc(n_layers, sizeof(double));
    state->theta_ice = (double *)calloc(n_layers, sizeof(double));

    if (!state->theta_liq || !state->theta_ice)
    {
        free(state->theta_liq);
        free(state->theta_ice);
        free(state);
        return NULL;
    }

    return state;
}

void free_soil_state(SoilState *state)
{
    if (state)
    {
        free(state->theta_liq);
        free(state->theta_ice);
        free(state);
    }
}

void init_soil_state(SoilState *state, double init_theta)
{
    if (!state)
        return;

    for (int i = 0; i < state->n_layers; i++)
    {
        state->theta_liq[i] = init_theta; 
        state->theta_ice[i] = 0.0;        
    }
}

void copy_soil_state(SoilState *dest, const SoilState *src, int n_layers)
{
    if (!dest || !src)
        return;

    memcpy(dest->theta_liq, src->theta_liq, n_layers * sizeof(double));
    memcpy(dest->theta_ice, src->theta_ice, n_layers * sizeof(double));
}