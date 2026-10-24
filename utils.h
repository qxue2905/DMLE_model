#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"

Dataset *create_dataset(int capacity);
void free_dataset(Dataset *ds);
int add_data_point(Dataset *ds, const DataPoint *point);
int read_csv_data(const char *filename, Dataset *ds);
int write_results_csv(const char *filename, const Dataset *ds, const ModelResult *result);

#endif /* UTILS_H */