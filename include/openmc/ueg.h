//! \file ueg.h
//! \brief Unionized Energy Grid Implementation

#ifndef OPENMC_UEG_H
#define OPENMC_UEG_H

#include "openmc/vector.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/memory.h"

namespace openmc {
  void create_union_energy_grid();
  void thin_union_energy_grid();
  void create_union_energy_xs();
} // namespace openmc

#endif // OPENMC_UEG_H
