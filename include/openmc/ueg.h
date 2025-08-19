//! \file ueg.h
//! \brief Unionized Energy Grid Implementation

#ifndef OPENMC_UEG_H
#define OPENMC_UEG_H

#include<vector>

#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"

namespace openmc {
  namespace data {
    extern std::shared_ptr<Nuclide::EnergyGrid> union_e_grid;
  }

  void create_union_energy_grid();
  void thin_union_energy_grid();
  void create_union_energy_xs();
} // namespace openmc

#endif // OPENMC_UEG_H
