#include<fstream>

#include<fmt/core.h>

#include "openmc/ueg.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"

namespace openmc {
  namespace data {
std::shared_ptr<Nuclide::EnergyGrid> union_e_grid;
  } // namespace data

  void create_union_energy_grid() {
    data::union_e_grid = std::make_shared<Nuclide::EnergyGrid>();
    vector<double> important_e_grid {};
    vector<double>& ueg = data::union_e_grid->energy;
    vector<int>& ueg_index = data::union_e_grid->grid_index;

    for (const auto& nuclide : data::nuclides) {
      for (int t = 0; t < nuclide->kTs_.size(); t++) {
         vector<double>& energies = nuclide->grid_[t].energy;
         ueg.insert(ueg.end(), energies.begin(), energies.end());

         if (nuclide->urr_present_) {
          const auto& urr_energies = nuclide->urr_data_[t].energy_;
          important_e_grid.insert(important_e_grid.end(), urr_energies.begin(), urr_energies.end());
        }
      }
      std::cout << "UEG Size : " << ueg.size() << std::endl;
    }
    std::sort(ueg.begin(), ueg.end());
    thin_union_energy_grid();
    
    ueg.insert(ueg.end(), important_e_grid.begin(), important_e_grid.end());
    std::sort(ueg.begin(), ueg.end());
    ueg.erase(std::unique(ueg.begin(), ueg.end()), ueg.end());

    int neutron = static_cast<int>(ParticleType::neutron);
    int M = settings::n_log_bins;

    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    double pseudo_spacing = std::log(E_max / E_min);
    auto umesh = xt::linspace(0.0, pseudo_spacing, M + 1);

    int ueg_size = ueg.size();
    ueg_index.resize(M + 1);
    int j = 0;

    for (int k = 0; k <= M; k++) {
      while (std::log(ueg[j + 1] / E_min) <= umesh[k]) {
       if (j + 2 == ueg_size) break;
       ++j;
    }
    ueg_index[k] = j;
    }
    std::cout << "UEG FINAL Size : " << ueg.size() << std::endl;
    create_union_energy_xs();
  }

  void thin_union_energy_grid() {
    vector<double>& ueg = data::union_e_grid->energy;
    double tau = settings::ue_grid_cutoff;
    int grid_size = 0;
    for (int i = 0; i < ueg.size() - 1; i++) {
      double current_e = ueg[i];
      double next_e = ueg[i + 1];
      if ((next_e - current_e) < tau*current_e) {
        ueg[grid_size] = 0.5*(current_e + next_e);
      } else {
        ueg[grid_size] = current_e;
        grid_size++;
      }
    } 
    ueg.resize(grid_size + 1);
    ueg.shrink_to_fit();
  }

  void create_union_energy_xs() {
    Nuclide::EnergyGrid & ueg = *data::union_e_grid;
    for (auto & nuclide : data::nuclides) {
      auto & xs = nuclide->xs_;
      auto & grid = nuclide->grid_;
    for (int t = 0; t < nuclide->kTs_.size(); t++) {
        auto & energy = grid[t].energy;
        array<size_t, 2> new_shape {ueg.energy.size(), 5};
        xs[t].resize(new_shape);
    
        auto temp_xs = xt::interp(xt::adapt(ueg.energy), xt::adapt(grid[t].energy), xs[t]);
        std::copy(temp_xs.begin(), temp_xs.end(), xs[t].begin());
        grid[t] = ueg;
        //majorant = xt::maximum(xs[t], majorant);
    }
  } 
}
} // namespace openmc
