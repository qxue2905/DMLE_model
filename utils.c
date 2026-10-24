#define _POSIX_C_SOURCE 200809L 

#include "utils.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <strings.h> 

Dataset *create_dataset(int capacity)
{
    Dataset *ds = (Dataset *)malloc(sizeof(Dataset));
    if (!ds)
        return NULL;

    ds->points = (DataPoint *)malloc(capacity * sizeof(DataPoint));
    if (!ds->points)
    {
        free(ds);
        return NULL;
    }

    ds->n = 0;
    ds->capacity = capacity;
    return ds;
}

void free_dataset(Dataset *ds)
{
    if (ds)
    {
        free(ds->points);
        free(ds);
    }
}

int add_data_point(Dataset *ds, const DataPoint *point)
{
    if (!ds || !point)
        return -1;

    if (ds->n >= ds->capacity)
    {
        int new_capacity = ds->capacity * 2;
        DataPoint *new_points = (DataPoint *)realloc(ds->points, new_capacity * sizeof(DataPoint));
        if (!new_points)
            return -1;
        ds->points = new_points;
        ds->capacity = new_capacity;
    }

    ds->points[ds->n++] = *point;
    return 0;
}

static char *trim(char *str)
{
    char *end;

    while (isspace((unsigned char)*str))
        str++;

    if (*str == 0)
        return str;

    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end))
        end--;

    *(end + 1) = '\0';
    return str;
}

int read_csv_data(const char *filename, Dataset *ds)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Error: 无法打开文件 %s\n", filename);
        return -1;
    }

    char line[2048];
    if (!fgets(line, sizeof(line), fp))
    {
        fclose(fp);
        return -1;
    }

    char *header = strdup(line);       
    char *token = strtok(header, ","); 
    int col_Ta = -1, col_qa = -1, col_U = -1, col_P = -1;
    int col_E_eff = -1, col_obs_LE = -1, col_obs_H = -1, col_rain = -1;
    int col_time = -1;
    int col_idx = 0;

    while (token)
    {
        char *p = trim(token);
        if (strcasecmp(p, "Ta") == 0 || strcasecmp(p, "T") == 0)
        {
            col_Ta = col_idx;
        }
        else if (strcasecmp(p, "qa") == 0 || strcasecmp(p, "Qa") == 0)
        {
            col_qa = col_idx;
        }
        else if (strcasecmp(p, "U") == 0 || strcasecmp(p, "Wind") == 0)
        {
            col_U = col_idx;
        }
        else if (strcasecmp(p, "P") == 0 || strcasecmp(p, "Press") == 0)
        {
            col_P = col_idx;
        }
        else if (strcasecmp(p, "E_eff") == 0 || strcasecmp(p, "Eeff") == 0)
        {
            col_E_eff = col_idx;
        }
        else if (strcasecmp(p, "obs_LE") == 0 || strcasecmp(p, "LE") == 0)
        {
            col_obs_LE = col_idx;
        }
        else if (strcasecmp(p, "obs_H") == 0 || strcasecmp(p, "H") == 0)
        {
            col_obs_H = col_idx;
        }
        else if (strcasecmp(p, "rain") == 0 || strcasecmp(p, "Rain") == 0)
        {
            col_rain = col_idx;
        }
        else if (strcasecmp(p, "time") == 0 || strcasecmp(p, "Time") == 0 ||
                 strcasecmp(p, "datetime") == 0 || strcasecmp(p, "Datetime") == 0)
        {
            col_time = col_idx;
        }

        token = strtok(NULL, ",");
        col_idx++;
    }
    free(header);

    if (col_Ta < 0 || col_qa < 0 || col_U < 0 || col_P < 0 ||
        col_E_eff < 0 || col_obs_LE < 0 || col_obs_H < 0)
    {
        fprintf(stderr, "Error: 缺少必要的列\n");
        fprintf(stderr, "需要的列: Ta, qa, U, P, E_eff, obs_LE, obs_H\n");
        fclose(fp);
        return -1;
    }

    int line_num = 0;
    while (fgets(line, sizeof(line), fp))
    {
        line_num++;
        DataPoint point;
        memset(&point, 0, sizeof(DataPoint));
        point.datetime = 0;

        char *line_copy = strdup(line);
        char *saveptr;
        char *token_data = strtok_r(line_copy, ",", &saveptr);
        int col = 0;

        while (token_data)
        {
            char *p = trim(token_data);

            double val = 0.0;
            char *endptr;
            val = strtod(p, &endptr);

            if (endptr == p || val == -9999.0 || val == -9999)
            {
                val = NAN;
            }

            if (col == col_Ta)
                point.Ta = val;
            else if (col == col_qa)
                point.qa = val;
            else if (col == col_U)
                point.U = val;
            else if (col == col_P)
                point.P = val;
            else if (col == col_E_eff)
                point.E_eff = val;
            else if (col == col_obs_LE)
                point.obs_LE = val;
            else if (col == col_obs_H)
                point.obs_H = val;
            else if (col == col_rain)
                point.rain = val;
            else if (col == col_time)
            {
                point.datetime = (long long)val;
            }

            token_data = strtok_r(NULL, ",", &saveptr);
            col++;
        }
        free(line_copy);

        if (!isnan(point.Ta) && !isnan(point.qa) && !isnan(point.U) &&
            !isnan(point.P) && !isnan(point.E_eff) &&
            !isnan(point.obs_LE) && !isnan(point.obs_H))
        {
            if (isnan(point.rain))
                point.rain = 0.0;
            add_data_point(ds, &point);
        }
    }

    fclose(fp);
    return 0;
}

int write_results_csv(const char *filename, const Dataset *ds, const ModelResult *result)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
        return -1;

    fprintf(fp, "datetime,Ta,qa,U,P,E_eff,obs_LE,obs_H,"
                "LE_sim,Delta,rho_a,Cp,Lv,Ra,q_s,hu,theta1,psi,Rs,"
                "ka,film_thickness,De,n_theta,theta_m,ka_formula\n");

    for (int i = 0; i < ds->n && i < result->n; i++)
    {
        fprintf(fp, "%lld,%.4f,%.6f,%.3f,%.1f,%.2f,%.2f,%.2f,"
                    "%.4f,%.4f,%.6f,%.4f,%.2f,%.4f,%.6f,%.6f,%.6f,%.6f,%.4f,"
                    "%.6e,%.6e,%.6e,%.6f,%.6f,\"%s\"\n",
                ds->points[i].datetime,
                ds->points[i].Ta,
                ds->points[i].qa,
                ds->points[i].U,
                ds->points[i].P,
                ds->points[i].E_eff,
                ds->points[i].obs_LE,
                ds->points[i].obs_H,
                result->LE_sim[i],
                result->Delta_sim[i],
                result->rho_a_sim[i],
                result->Cp_sim[i],
                result->Lv_sim[i],
                result->Ra_sim[i],
                result->qs_sim[i],
                result->hu_sim[i],
                result->theta1_sim[i],
                result->psi_sim[i],
                result->Rs_sim[i],
                result->ka_sim[i],
                result->film_thickness_sim[i],
                result->De_sim[i],
                result->n_theta_sim[i],
                result->theta_m_sim[i],
                result->ka_formula_sim[i]);   
    }

    fclose(fp);
    return 0;
}
