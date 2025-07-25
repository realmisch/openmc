#include<fstream>

#include<fmt.core>

#include "openmc/ueg.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

namespace openmc {
  namespace data {
    unique_ptr<UnionEnergyGrid> union_e_grid;
  } // namespace data

  void create_union_energy_grid() {
    data::union_e_grid = make_unique<UnionEnergyGrid>();
    vector<double> important_e_grid {};
    vector<double>& ueg = data::union_e_grid->energy;
    for (const auto& nuclide : data::nuclides) {
       vector<double>& energies = nuclide->grid_[0].energy;
       ueg->insert(ueg->end(), energies.begin(), energies.end());

       if (nuclide->urr_present_) {
        const auto& urr_energies = nuclide->urr_data_[0].energy_;
        important_e_grid.insert(important_e_grid.end(), urr_energies->begin(), urr_energies->end());

    } 
  }
  }

  void thin_union_energy_grid() {
    vector<double>& ueg = data::union_e_grid->energy;
    double tau = settings::ueg_grid_cutoff;
    int grid_size = 0;
    for (int i = 0; i < ueg->size() - 1; i++) {
      double current_e = ueg[i];
      double next_e = ueg[i + 1];
      if ((current_e - next_e)/current_e) < tau) {
        ueg[grid_size] = 0.5*(current_e + next_e);
      } else {
        ueg[grid_size] = current_e;
        grid_size++;
      }
    } 
    ueg->resize(grid_size + 1);
    ueg->shrink_to_fit();
  }
} // namespace openmc
