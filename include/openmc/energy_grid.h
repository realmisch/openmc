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

  template<typename T>
  struct GridVector : vector<T> {
    T& operator [](int index) {
      if (data::use_ueg) return *data::union_e_grid;
      return this->data()[index];
    }

    const T& operator[](int index) const {
      if (data::use_ueg) return *data::union_e_grid;
      return this->data()[index];  
    }
  };
} //namespace openmc

#endif //OPENMC_ENERGY_GRID_H
