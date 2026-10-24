#ifndef PLOT_H
#define PLOT_H

#include "model.h"
#include "utils.h"

int plot_results(const Dataset *ds, const ModelResult *result,
                 const char *output_dir);

#endif /* PLOT_H */