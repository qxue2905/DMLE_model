#include "plot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int plot_results(const Dataset *ds, const ModelResult *result, const char *output_dir)
{
    if (!ds || !result || ds->n == 0)
        return -1;

    char cmd[1024];
#ifdef _WIN32
    snprintf(cmd, sizeof(cmd), "if not exist %s mkdir %s", output_dir, output_dir);
#else
    snprintf(cmd, sizeof(cmd), "mkdir -p %s", output_dir);
#endif
    system(cmd);

    int n_valid = 0;
    for (int i = 0; i < ds->n && i < result->n; i++)
    {
        double obs = ds->points[i].obs_LE;
        double sim = result->LE_sim[i];
        if (!isnan(obs) && !isnan(sim) && isfinite(obs) && isfinite(sim))
            n_valid++;
    }
    if (n_valid < 5)
        return -1;   

    char data_file[512];
    snprintf(data_file, sizeof(data_file), "%s/plot_data.tmp", output_dir);

    FILE *fp = fopen(data_file, "w");
    if (!fp)
        return -1;

    for (int i = 0; i < ds->n && i < result->n; i++)
    {
        double obs = ds->points[i].obs_LE;
        double sim = result->LE_sim[i];
        if (!isnan(obs) && !isnan(sim) && isfinite(obs) && isfinite(sim))
        {
            fprintf(fp, "%d %.4f %.4f %.4f %.6e\n",
                    i,
                    obs,
                    sim,
                    result->theta1_sim[i] * 100.0,   
                    result->Rs_sim[i]);              
        }
    }
    fclose(fp);

    char plot_file[512];
    snprintf(plot_file, sizeof(plot_file), "%s/combined_plot.gnu", output_dir);

    fp = fopen(plot_file, "w");
    if (!fp)
    {
        remove(data_file);
        return -1;
    }

    fprintf(fp, "set terminal pngcairo size 1600,1200 font 'Arial,12'\n");
    fprintf(fp, "set output '%s/combined_plot.png'\n", output_dir);
    fprintf(fp, "set multiplot layout 2,2 title 'Bare Soil Evaporation: Model vs Observation' font 'Arial,16'\n");
    fprintf(fp, "set title 'Time Series: LE' font 'Arial,13'\n");
    fprintf(fp, "set xlabel 'Time Step'\n");
    fprintf(fp, "set ylabel 'LE (W/m^2)'\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set key top right\n");
    fprintf(fp, "plot '%s' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Observed LE', \\\n", data_file);
    fprintf(fp, "     '%s' using 1:3 with lines lw 2 lc rgb '#1f77b4' title 'Simulated LE'\n", data_file);
    fprintf(fp, "set title 'Scatter: Simulated vs Observed' font 'Arial,13'\n");
    fprintf(fp, "set xlabel 'Observed LE (W/m^2)'\n");
    fprintf(fp, "set ylabel 'Simulated LE (W/m^2)'\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set key top left\n");
    fprintf(fp, "stats '%s' using 2:3 nooutput\n", data_file);
    fprintf(fp, "f(x) = a*x + b\n");
    fprintf(fp, "fit f(x) '%s' using 2:3 via a,b\n", data_file);
    fprintf(fp, "set label 1 sprintf('N = %%d', STATS_records) at graph 0.05, graph 0.95 font 'Arial,11'\n");
    fprintf(fp, "set label 2 sprintf('R^2 = %%.4f', STATS_correlation**2) at graph 0.05, graph 0.90 font 'Arial,11'\n");
    fprintf(fp, "set label 3 sprintf('RMSE = %%.2f W/m^2', sqrt(STATS_sse/STATS_records)) at graph 0.05, graph 0.85 font 'Arial,11'\n");
    fprintf(fp, "plot x with lines lt 1 lc rgb '#ff0000' lw 2 title '1:1', \\\n");
    fprintf(fp, "     f(x) with lines lt 2 lc rgb '#0000ff' lw 2 title sprintf('y=%%.2fx+%%.2f', a, b), \\\n");
    fprintf(fp, "     '%s' using 2:3 with points pt 7 ps 1.2 lc rgb '#2ca02c' title 'Data'\n", data_file);
    fprintf(fp, "set title 'Surface Resistance Rs' font 'Arial,13'\n");
    fprintf(fp, "set xlabel 'Time Step'\n");
    fprintf(fp, "set ylabel 'Rs (s/m)'\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set key top right\n");
    fprintf(fp, "set yrange [0:*]\n");
    fprintf(fp, "plot '%s' using 1:5 with lines lw 2 lc rgb '#9467bd' title 'Rs'\n", data_file);
    fprintf(fp, "set title 'Surface Soil Moisture theta1' font 'Arial,13'\n");
    fprintf(fp, "set xlabel 'Time Step'\n");
    fprintf(fp, "set ylabel 'theta1 (%%)'\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set key top right\n");
    fprintf(fp, "plot '%s' using 1:4 with lines lw 2 lc rgb '#ff7f0e' title 'theta1'\n", data_file);

    fprintf(fp, "unset multiplot\n");
    fprintf(fp, "unset output\n");
    fprintf(fp, "quit\n");
    fclose(fp);

    int gnuplot_available;
#ifdef _WIN32
    gnuplot_available = system("gnuplot --version > nul 2>&1");
#else
    gnuplot_available = system("which gnuplot > /dev/null 2>&1");
#endif
    if (gnuplot_available != 0)
    {
        remove(data_file);
        remove(plot_file);
        return -1;
    }

    char cmd_plot[1024];
#ifdef _WIN32
    snprintf(cmd_plot, sizeof(cmd_plot), "gnuplot \"%s\" 2>nul", plot_file);
#else
    snprintf(cmd_plot, sizeof(cmd_plot), "gnuplot \"%s\" 2>/dev/null", plot_file);
#endif
    int ret = system(cmd_plot);

    remove(data_file);
    remove(plot_file);

    return (ret == 0) ? 0 : -1;
}