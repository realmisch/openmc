//! \file ueg.h
//! \brief Unionized Energy Grid Implementation

#ifndef OPENMC_UEG_H
#define OPENMC_UEG_H

#include<vector>

#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"

namespace openmc {
  struct UnionEnergyGrid {
    // init method
    void init();
    // data members
    std::vector<int> grid_index;
    std::vector<double> energy;

  };

  namespace data {
    extern std::unique_ptr<UnionEnergyGrid> union_e_grid;
  }

} // namespace openmc

#endif // OPENMC_UEG_H
