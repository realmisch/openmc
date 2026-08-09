
#include "openmc/energy_grid.h"
#include "openmc/nuclide.h"

namespace openmc {
  void EnergyGrid::insert_grid(const vector<double> &other, const bool sort_result) {
    energy.insert(energy.end(), other.begin(), other.end());
    if (sort_result)
      std::inplace_merge(std::execution::par_unseq,
                         energy.begin(),
                         energy.begin() + (energy.size() - other.size()),
                         energy.end());
    energy.erase(std::unique(std::execution::par_unseq, 
                             energy.begin(), 
                             energy.end()), 
                 energy.end());
  }

  

  void EnergyGrid::thin_grid(double tolerance) {
    int grid_size = 0;
    for (int i = 0; i < energy.size() - 1; ++i) {
      double current_erg = energy[i];
      double next_erg = energy[i + 1];

      if ((next_erg - current_erg) < tolerance * current_erg)
        energy[grid_size] = 0.5 * (current_erg + next_erg);
      else {
        energy[grid_size] = current_erg;
        grid_size++;
      }
    }
    energy.resize(grid_size + 1);
  }

  void EnergyGrid::update_dix_and_bound() {
    int neutron = ParticleType::neutron().transport_index();
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    auto min_it = energy.begin();
    auto max_it = --energy.end();

    while (*min_it < E_min) min_it++;
    while (*max_it > E_max) max_it--;

    if (max_it + 1 != energy.end())
      energy.erase(max_it - 1, energy.end());
    if (min_it != energy.begin())
      energy.erase(energy.begin(), min_it + 1);

    int M = settings::n_log_bins;
    auto umesh = tensor::linspace(0.0, std::log(E_max / E_min), M + 1);

    int grid_size = energy.size();
    grid_index.resize(M + 1);
    
    int j = 0;
    for (int k = 0; k <= M; ++k) {
      while (std::log(energy[j + 1] / E_min) <= umesh[k]) {
        if (j + 2 == grid_size) break;
        ++j;
      }
      grid_index[k] = j;
    }
  }
}

