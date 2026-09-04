#include<fmt/core.h>
#include<execution>
#include<set>

#include "openmc/energy_grid.h"
#include "openmc/ueg.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/search.h"
#include "openmc/memory.h"
#include "openmc/vector.h"
#include "openmc/tensor.h"

#include "openmc/reaction.h"
#include "openmc/message_passing.h"

namespace openmc {
  namespace data {
    bool use_ueg = false;
    std::shared_ptr<EnergyGrid> ue_grid;
  } // namespace data

  void create_union_energy_grid() {

    int neutron = ParticleType::neutron().transport_index();
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    int M = settings::n_log_bins;
    double log_range = std::log(E_max / E_min);
    auto log_mesh = tensor::linspace(0.0, log_range, M + 1);

    settings::energy_cutoff[0] = std::max(E_min, settings::energy_cutoff[0]);
    
    vector<double>& ueg = data::ue_grid->energy;
    vector<int>& ueg_index = data::ue_grid->grid_index;

    double mem_size = unionize_nuclides();
    ueg_index.resize(M + 1); 

    vector<double> bin_energy(M + 1);
    for (int k = 0; k <= M; ++k)
      bin_energy[k] = E_min * std::exp(log_mesh[k]);

    int j = 0;
    for (int k = 0; k <= M; ++k) {
      while (ueg[j + 1] <= bin_energy[k]) {
        if (j + 2 == ueg.size()) break;
        ++j;
      }
      ueg_index[k] = j;
    }
    write_message("Global Unionized Energy Grid: {} grid points - {:.3f} GB of memory", ueg.size(), mem_size);
    if (mem_size > 10)
      warning(fmt::format("{} GB required for Unionized Energy Grid cross sections", mem_size));
    data::use_ueg = true;
  }
  
  double unionize_nuclides() 
  {
    int neutron = ParticleType::neutron().transport_index();
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    auto& ueg = data::ue_grid->energy;
    //imp_e_grid will contain energy points that should not be thinned (URR and Sab energies)
    vector<double> imp_e_grid {E_min, E_max};

    int total_energies = 0.0;
    for (const auto& nuc : data::nuclides)
      for (int t = 0; t < nuc->kTs_.size(); t++)
        total_energies += nuc->grid_[t].energy.size();

    ueg.reserve(total_energies);

    //Get energies from all nuclides
    for (const auto& nuc : data::nuclides) {
      for (int t = 0; t < nuc->kTs_.size(); t++) {
        const vector<double>& energies = nuc->grid_[t].energy;
        ueg.insert(ueg.end(), energies.begin(), energies.end());

        /*
        //Add URR energies to important energy grid
        if (nuclide->urr_present_) {
        const auto& urr_energies = nuclide->urr_data_[t].energy_;
        imp_e_grid.insert(imp_e_grid.end(), urr_energies.begin(), urr_energies.end());
        }
        */
        //Add threshold energies to important energy grid
        for (auto& rxn : nuc->reactions_) {
          imp_e_grid.insert(imp_e_grid.end(), energies[rxn->xs_[t].threshold]);
          imp_e_grid.insert(imp_e_grid.end(), energies.back());
        }
      }
    }
    
    std::sort(std::execution::par_unseq, ueg.begin(), ueg.end());
    ueg.erase(std::unique(std::execution::par_unseq, ueg.begin(), ueg.end()), ueg.end());
    //Sort ueg energy and thin redundant points according to the thinning cutoff parameter
    double tau = settings::ue_grid_cutoff;
    int grid_size = 0;
    int i = 0;
    while (i < ueg.size()) {
      int end = i;
      while (end + 1 < ueg.size() && (ueg[end + 1] - ueg[end]) < tau * ueg[end])
        ++end;
      ueg[grid_size++] = (end > i) ? 0.5 * (ueg[i] + ueg[end]) : ueg[i];
      i = end + 1;
    } 
    ueg.resize(grid_size + 1);
    ueg.shrink_to_fit();

    //Insert important grid points
    ueg.insert(ueg.end(), imp_e_grid.begin(), imp_e_grid.end());
    std::sort(std::execution::par_unseq, ueg.begin(), ueg.end());
    
    auto min_it = ueg.begin();
    auto max_it = ueg.end() - 1;

    while (*min_it < E_min) min_it++;
    while (*max_it > E_max) max_it--;

    ueg.erase(max_it + 1, ueg.end());
    ueg.erase(ueg.begin(), min_it);

    ueg.erase(std::unique(std::execution::par_unseq, ueg.begin(), ueg.end()), ueg.end());
    //std::sort(std::execution::par_unseq, ueg.begin(), ueg.end());
    
    const tensor::View<const double> e(ueg.data(), {ueg.size()}, {1});

    struct XsUpdateMap {
      int nuc_idx;
      int rxn_idx;
      int t;
    };

    vector<XsUpdateMap> tasks;
    int num_temps = 0;
    for (int n = 0; n < data::nuclides.size(); ++n) {
      auto& nuc = data::nuclides[n];
      num_temps += (nuc->reactions_.size() * nuc->kTs_.size());
      for (int rx = 0; rx < nuc->reactions_.size(); ++rx)
        for (int t = 0; t < nuc->kTs_.size(); ++t)
          tasks.push_back({n, rx, t});
    }

    #pragma omp parallel for
    for (int i_task = 0; i_task < tasks.size(); ++i_task) {
      const auto& task = tasks[i_task];
      auto& nuc = data::nuclides[task.nuc_idx];
      auto& rxn = nuc->reactions_[task.rxn_idx];
      auto& grid = nuc->grid_;
      const int t = task.t;

      auto& xs = rxn->xs_[t];
      const size_t n_energies = grid[t].energy.size();

      const tensor::View<const double> grid_data(grid[t].energy.data(), {n_energies}, {1});
      auto ep = grid_data.slice(tensor::range(xs.threshold, n_energies));

      const tensor::View<const double> xsp(xs.value.data(), {xs.value.size()}, {1});

      if (xs.threshold != 0)
        xs.threshold = lower_bound_index(ueg.begin(), ueg.end(), ep[0]);

      auto e_grid = e.slice(tensor::range(xs.threshold, e.size()));
      auto rxn_xs = tensor::interp(e_grid, ep, xsp, 0.0, xs.value.back());
      xs.value = vector<double>(rxn_xs.cbegin(), rxn_xs.cend());
    }

    for (auto& nuc : data::nuclides) {
      nuc->grid_.clear();
      nuc->create_ue_derived(nuc->prompt_photons_.get(), nuc->delayed_photons_.get(), ueg);
    }

    write_message("Num points : {} | Num temps : {} | sizeof double : {}", std::to_string(ueg.size()), std::to_string(num_temps), std::to_string(sizeof(double)));
    double mem_size = (double)(ueg.size()*num_temps)*sizeof(double)*BYTES_TO_GIGABYTES;
    return mem_size;
  }
} // namespace openmc w
