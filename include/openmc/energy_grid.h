#ifndef OPENMC_ENERGY_GRID_H
#define OPENMC_ENERGY_GRID_H

#include <algorithm>
#include <execution>

#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/vector.h"
#include "openmc/tensor.h"

namespace openmc {
  struct EnergyGrid {
    vector<int> grid_index;
    vector<double> energy;

    void insert_grid(const vector<double> &other);
    
    void thin_grid(double tolerance);
    void update_dix_and_bound();
  };
} //namespace openmc

#endif //OPENMC_ENERGY_GRID_H
