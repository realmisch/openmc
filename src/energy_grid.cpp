
#include "openmc/energy_grid.h"

namespace openmc {
  void EnergyGrid::insert_grid(const &EnergyGrid other) {
    vector<double> & other_energy = other.energy;

    energy.insert(energy.end(), other_energy.begin(), other_energy.end());
    energy.erase(std::unique(std::execution_par_unseq, energy.begin(), energy.end()), energy.end());
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

  void EnergyGrid::update_double_index() {
    int neutron = ParticleType::neutron().transport_index();
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

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

