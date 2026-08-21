#ifndef OPENMC_ENERGY_GRID_H
#define OPENMC_ENERGY_GRID_H

#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

namespace openmc {
  struct EnergyGrid {
    vector<int> grid_index;
    vector<double> energy;
  };
  namespace data {
      extern bool use_ueg;
      extern std::shared_ptr<EnergyGrid> union_e_grid;
    } //namespace data
} //namespace openmc

#endif //OPENMC_ENERGY_GRID_H
