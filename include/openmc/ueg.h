//! \file ueg.h
//! \brief Unionized Energy Grid Implementation

#ifndef OPENMC_UEG_H
#define OPENMC_UEG_H

#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/memory.h"

namespace openmc {
  struct XsUpdateMap {
    int nuc_idx;
    int rxn_idx;
    int t;
  };

  void create_union_energy_grid();
  double unionize_nuclides();
  void unionize_nuclide_idx();
} // namespace openmc

#endif // OPENMC_UEG_H
