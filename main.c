#include "constants.h"
#include "soil_params.h"
#include "soil_layers.h"
#include "model.h"
#include "utils.h"
#include "plot.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

ModelResult *create_model_result(int n);
void free_model_result(ModelResult *result);

#define N_LAYERS_DEFAULT 14

static int read_param_file(const char *filename, double *depth, double *h_flx,
                           double *coupling_factor)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: 无法打开参数文件 %s\n", filename);
        return -1;
    }

    char line[256];
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Error: 参数文件为空\n");
        fclose(fp);
        return -1;
    }
    fclose(fp);

    double d = 0.0, h = 0.0, cf = 0.0;
    if (sscanf(line, "%lf %lf %lf", &d, &h, &cf) != 3) {
        fprintf(stderr, "Error: 参数文件格式错误（应为 3 个数）\n");
        return -1;
    }

    if (d <= 0) { fprintf(stderr, "Error: 地下水位深度必须为正数\n"); return -1; }
    if (h <= 0) { fprintf(stderr, "Error: 通量测量高度必须为正数\n"); return -1; }
    if (cf < 0 || cf > 1) {
        fprintf(stderr, "Warning: 耦合系数应在 [0,1] 范围内，使用默认值 0.8\n");
        cf = 0.8;
    }

    *depth = d;
    *h_flx = h;
    *coupling_factor = cf;

    return 0;
}

static int read_soil_texture(const char *filename, double *sand, double *clay, int n_layers)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Error: 无法打开土壤质地文件 %s\n", filename);
        return -1;
    }

    int count = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp) && count < n_layers)
    {
        double s, c;
        if (sscanf(line, "%lf %lf", &s, &c) == 2)
        {
            if (s < 0 || s > 1 || c < 0 || c > 1)
            {
                fprintf(stderr, "Error: 土壤质地值超出[0,1]范围（行 %d）\n", count + 1);
                fclose(fp);
                return -1;
            }
            sand[count] = s;
            clay[count] = c;
            count++;
        }
    }
    fclose(fp);

    if (count != n_layers)
    {
        fprintf(stderr, "Error: 土壤质地文件应包含 %d 行，实际读取 %d 行\n", n_layers, count);
        return -1;
    }
    return 0;
}

static void print_usage(const char *prog)
{
    printf("\n用法: %s -soil <soil_file> [-param <param_file>] [data_file]\n", prog);
    printf("  -soil <file>    : 土壤质地文件，每行两个数：砂粒 黏粒 (0~1)，共14行\n");
    printf("  -param <file>   : 参数文件 (默认 param.txt)，3 个数：地下水位深度 测量高度 耦合系数\n");
    printf("  data_file       : 观测数据CSV文件 (可选，默认 data/flux_data.csv)\n");
    printf("\n");
}

static void print_evaluation_panel(const Dataset *ds, const ModelResult *result)
{
    int n = 0;
    double sum_obs = 0.0;

    for (int i = 0; i < ds->n && i < result->n; i++)
    {
        double obs = ds->points[i].obs_LE;
        double sim = result->LE_sim[i];
        if (!isnan(obs) && !isnan(sim) && isfinite(obs) && isfinite(sim))
        {
            sum_obs += obs;
            n++;
        }
    }

    if (n < 2)
    {
        printf("有效数据点不足，无法计算评估指标\n");
        return;
    }

    double mean_obs = sum_obs / n;

    double sse = 0.0; 
    double sst = 0.0;  

    for (int i = 0; i < ds->n && i < result->n; i++)
    {
        double obs = ds->points[i].obs_LE;
        double sim = result->LE_sim[i];
        if (!isnan(obs) && !isnan(sim) && isfinite(obs) && isfinite(sim))
        {
            sse += (obs - sim) * (obs - sim);
            sst += (obs - mean_obs) * (obs - mean_obs);
        }
    }

    double rmse = sqrt(sse / n);
    double r2   = (sst > 1e-10) ? (1.0 - sse / sst) : 0.0;

    printf("\n========================================\n");
    printf("  模型评估结果（决定系数定义）\n");
    printf("========================================\n");
    printf("  数据量 N   = %d\n", n);
    printf("  R²         = %.4f\n", r2);
    printf("  RMSE       = %.4f W/m^2\n", rmse);
    printf("========================================\n\n");
}

int run_model_with_params(const Dataset *ds, double dt, int n_layers,
                          const WaterTable *wt, const SoilParams *soil_params,
                          ModelResult **result_out)
{
    SoilLayers *layers = create_soil_layers(n_layers);
    if (!layers)
        return -1;

    SoilState *state = create_soil_state(n_layers);
    if (!state)
    {
        free_soil_layers(layers);
        return -1;
    }
    init_soil_state(state, 0.10);   

    ModelResult *result = create_model_result(ds->n);
    if (!result)
    {
        free_soil_layers(layers);
        free_soil_state(state);
        return -1;
    }

    double actual_dt = (dt > 0) ? dt : DEFAULT_DT;

    for (int i = 0; i < ds->n; i++)
    {
        DataPoint *dp = &ds->points[i];

        Forcing forcing;
        forcing.Ta     = dp->Ta;
        forcing.qa     = dp->qa;
        forcing.U      = dp->U;
        forcing.P      = dp->P;
        forcing.E_eff  = dp->E_eff;
        forcing.obs_LE = dp->obs_LE;
        forcing.obs_H  = dp->obs_H;
        forcing.rain   = dp->rain;

        FluxResult flux_result;
        calculate_baresoil_evaporation(soil_params, state, layers, &forcing,
                                       wt, actual_dt, &flux_result);
        double horton = 0.0;
        update_soil_moisture_wt(soil_params, state, layers, forcing.rain,
                                flux_result.LE / (flux_result.Lv * RHO_W),
                                wt, actual_dt, &horton);
        result->LE_sim[i]             = flux_result.LE;
        result->Delta_sim[i]          = flux_result.Delta;
        result->rho_a_sim[i]          = flux_result.rho_a;
        result->Cp_sim[i]             = flux_result.Cp;
        result->Lv_sim[i]             = flux_result.Lv;
        result->Ra_sim[i]             = flux_result.Ra;
        result->qs_sim[i]             = flux_result.qs;
        result->hu_sim[i]             = flux_result.hu;
        result->theta1_sim[i]         = state->theta_liq[0];  
        result->psi_sim[i]            = flux_result.psi;
        result->Rs_sim[i]             = flux_result.Rs;
        result->ka_sim[i]             = flux_result.ka;
        result->film_thickness_sim[i] = flux_result.film_thickness;
        result->De_sim[i]             = flux_result.De;
        result->n_theta_sim[i]        = flux_result.n_theta;
        result->theta_m_sim[i]        = flux_result.theta_m;
		snprintf(result->ka_formula_sim[i],
		         sizeof(result->ka_formula_sim[i]),
		         "%s", flux_result.ka_formula);
    }

    free_soil_layers(layers);
    free_soil_state(state);

    *result_out = result;
    return 0;
}

ModelResult *create_model_result(int n)
{
    ModelResult *result = (ModelResult *)malloc(sizeof(ModelResult));
    if (!result)
        return NULL;

    result->n = n;
    result->LE_sim             = (double *)calloc(n, sizeof(double));
    result->Delta_sim          = (double *)calloc(n, sizeof(double));
    result->rho_a_sim          = (double *)calloc(n, sizeof(double));
    result->Cp_sim             = (double *)calloc(n, sizeof(double));
    result->Lv_sim             = (double *)calloc(n, sizeof(double));
    result->Ra_sim             = (double *)calloc(n, sizeof(double));
    result->qs_sim             = (double *)calloc(n, sizeof(double));
    result->hu_sim             = (double *)calloc(n, sizeof(double));
    result->theta1_sim         = (double *)calloc(n, sizeof(double));
    result->psi_sim            = (double *)calloc(n, sizeof(double));
    result->Rs_sim             = (double *)calloc(n, sizeof(double));
    result->ka_sim             = (double *)calloc(n, sizeof(double));
    result->film_thickness_sim = (double *)calloc(n, sizeof(double));
    result->De_sim             = (double *)calloc(n, sizeof(double));
    result->n_theta_sim        = (double *)calloc(n, sizeof(double));
    result->theta_m_sim        = (double *)calloc(n, sizeof(double));
	result->ka_formula_sim = (char (*)[512])calloc(n, 512);

    if (!result->LE_sim || !result->Delta_sim || !result->rho_a_sim ||
        !result->Cp_sim || !result->Lv_sim || !result->Ra_sim ||
        !result->qs_sim || !result->hu_sim || !result->theta1_sim ||
        !result->psi_sim || !result->Rs_sim || !result->ka_sim ||
        !result->film_thickness_sim || !result->De_sim ||
        !result->n_theta_sim || !result->theta_m_sim ||
        !result->ka_formula_sim)
    {
        free_model_result(result);
        return NULL;
    }

    return result;
}

void free_model_result(ModelResult *result)
{
    if (result)
    {
        free(result->LE_sim);
        free(result->Delta_sim);
        free(result->rho_a_sim);
        free(result->Cp_sim);
        free(result->Lv_sim);
        free(result->Ra_sim);
        free(result->qs_sim);
        free(result->hu_sim);
        free(result->theta1_sim);
        free(result->psi_sim);
        free(result->Rs_sim);
        free(result->ka_sim);
        free(result->film_thickness_sim);
        free(result->De_sim);
        free(result->n_theta_sim);
        free(result->theta_m_sim);
		free(result->ka_formula_sim); 
        free(result);
    }
}

int main(int argc, char *argv[])
{
    char *soil_file = NULL;
    char *data_file = "data/flux_data.csv";   
    char *param_file = "param.txt";          

    int i = 1;
    while (i < argc)
    {
        if (strcmp(argv[i], "-soil") == 0 && i + 1 < argc)
            soil_file = argv[++i];
        else if (strcmp(argv[i], "-param") == 0 && i + 1 < argc)
            param_file = argv[++i];
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            print_usage(argv[0]);
            return 0;
        }
        else if (argv[i][0] != '-')
            data_file = argv[i];
        else
        {
            print_usage(argv[0]);
            return -1;
        }
        i++;
    }

    if (!soil_file)
    {
        print_usage(argv[0]);
        return -1;
    }

    double wtd_depth = 100.0;
    double h_flx = 5.0;
    double coupling_factor = 0.8;

    read_param_file(param_file, &wtd_depth, &h_flx, &coupling_factor);

    const int N_LAYERS = N_LAYERS_DEFAULT;
    double *sand = (double *)malloc(N_LAYERS * sizeof(double));
    double *clay = (double *)malloc(N_LAYERS * sizeof(double));
    if (!sand || !clay)
    {
        free(sand);
        free(clay);
        return -1;
    }

    if (read_soil_texture(soil_file, sand, clay, N_LAYERS) != 0)
    {
        free(sand);
        free(clay);
        return -1;
    }

    SoilParams *soil_params = (SoilParams *)malloc(N_LAYERS * sizeof(SoilParams));
    if (!soil_params)
    {
        free(sand);
        free(clay);
        return -1;
    }
    for (int layer = 0; layer < N_LAYERS; layer++)
        calc_soil_params_from_texture(sand[layer], clay[layer], &soil_params[layer]);

    WaterTable wt;
    init_water_table(&wt, wtd_depth, coupling_factor, 1);
    wt.h_flx = h_flx;

    Dataset *ds = create_dataset(10000);
    if (!ds)
    {
        free(sand);
        free(clay);
        free(soil_params);
        return -1;
    }

    if (read_csv_data(data_file, ds) != 0 || ds->n == 0)
    {
        free_dataset(ds);
        free(sand);
        free(clay);
        free(soil_params);
        return -1;
    }

    ModelResult *result = NULL;
    int ret = run_model_with_params(ds, DEFAULT_DT, N_LAYERS, &wt,
                                    soil_params, &result);
    if (ret != 0 || !result)
    {
        free_dataset(ds);
        free(sand);
        free(clay);
        free(soil_params);
        return -1;
    }

    print_evaluation_panel(ds, result);

    ret = write_results_csv("output/rs_results.csv", ds, result);
    if (ret == 0)
        printf("结果已保存到 output/rs_results.csv\n");

    ret = plot_results(ds, result, "output");
    if (ret == 0)
        printf("图像已保存到 output/combined_plot.png\n");

    free_dataset(ds);
    free_model_result(result);
    free(sand);
    free(clay);
    free(soil_params);

    return 0;
}