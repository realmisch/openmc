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

#include "openmc/constants.h"
#include "openmc/reaction.h"

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
    
    double mem_size = 0.0;
    if (settings::ue_grid_method == UnionizationMethod::MATERIAL) {
      int max_grid_size = 1;
      for (auto& mat : model::materials) {
        vector<double>& ueg = mat->ue_grid_.energy;
        vector<int>& ueg_index = mat->ue_grid_.grid_index;

        vector<Nuclide *> temp_nuclides;
        for (auto& nuc : mat->nuclide_)
          temp_nuclides.push_back(data::nuclides[nuc].get());

        mem_size += unionize_nuclides(temp_nuclides, ueg);
        ueg_index.resize(M + 1); 
        int j = 0;
        for (int k = 0; k <= M; ++k) {
          while (std::log(ueg[j + 1] / E_min) <= log_mesh[k]) {
            if (j + 2 == ueg.size()) break;
            ++j;
          }
          ueg_index[k] = j;
        }
        max_grid_size = ueg.size() > max_grid_size ? ueg.size() : max_grid_size;
      }
      write_message("Material-Wise Energy Grid Unionization performed with up to {} grid points", std::to_string(max_grid_size));
    } else if (settings::ue_grid_method == UnionizationMethod::GLOBAL) {
      vector<double>& ueg = data::ue_grid->energy;
      vector<int>& ueg_index = data::ue_grid->grid_index;
        
      vector<Nuclide *> temp_nuclides;
      for (auto& nuc : data::nuclides)
        temp_nuclides.push_back(nuc.get());

      mem_size = unionize_nuclides(temp_nuclides, ueg);
      ueg_index.resize(M + 1); 
      int j = 0;
      for (int k = 0; k <= M; ++k) {
        while (std::log(ueg[j + 1] / E_min) <= log_mesh[k]) {
          if (j + 2 == ueg.size()) break;
          ++j;
        }
        ueg_index[k] = j;
      }
      write_message("Global Energy Grid performed with {} grid points", ueg.size());
    }
    if (mem_size > 10)
      warning(fmt::format("{} GB required for Unionized Energy Grid cross sections", mem_size));
    data::use_ueg = true;
  }
  
  double unionize_nuclides(vector<Nuclide *>& nuclides, vector<double>& ueg) 
  {
    int neutron = ParticleType::neutron().transport_index();
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    //imp_e_grid will contain energy points that should not be thinned (URR and Sab energies)
    vector<double> imp_e_grid {E_min, E_max};

    int num_temps;
    //Get energies from all nuclides
    for (const auto& nuc : nuclides) {
      num_temps += nuc->kTs_.size();
      for (int t = 0; t < nuc->kTs_.size(); t++) {
         vector<double>& energies = nuc->grid_[t].energy;
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
    //Sort ueg energy and thin redundant points according to the thinning cutoff parameter
    double tau = settings::ue_grid_cutoff;
    int grid_size = 0;
    for (int i = 0; i < ueg.size() - 1; i++) {
      double current_e = ueg[i];
      double next_e = ueg[i + 1];
      if ((next_e - current_e) < tau * current_e) {
        ueg[grid_size] = 0.5*(current_e + next_e);
      } else {
        ueg[grid_size] = current_e;
        grid_size++;
      }
    }
    ueg.resize(grid_size + 1);
    ueg.shrink_to_fit();

    //Insert important grid points
    ueg.insert(ueg.end(), imp_e_grid.begin(), imp_e_grid.end());
    std::sort(std::execution::par_unseq, ueg.begin(), ueg.end());
    
    auto min_it = ueg.begin();
    auto max_it = --ueg.end();

    while (*min_it < E_min) min_it++;
    while (*max_it > E_max) max_it--;

    ueg.erase(ueg.begin(), min_it + 1);
    ueg.erase(max_it - 1, ueg.end());
    ueg.erase(std::unique(std::execution::par_unseq, ueg.begin(), ueg.end()), ueg.end());
    std::sort(std::execution::par_unseq, ueg.begin(), ueg.end());
    
    const auto e = tensor::Tensor<double>(ueg.data(), ueg.size());

    //Iterate through all nuclides to update XS
    for (auto & nuc : nuclides) {
      auto & grid = nuc->grid_;
      //Interpolate XS for each nuclide temperature and reaction
      for (auto& rxn : nuc->reactions_) {
        for (int t = 0; t < nuc->kTs_.size(); t++) {
          auto & xs = rxn->xs_[t];
          auto n_energies = grid[t].energy.size();
          auto grid_data = tensor::Tensor<double>(grid[t].energy.data(), n_energies);
          auto ep = grid_data.slice(tensor::range(xs.threshold, n_energies));
          
          auto xsp = tensor::Tensor<double>(xs.value.data(), xs.threshold + xs.value.size());
          if (xs.threshold != 0)
            xs.threshold = lower_bound_index(ueg.begin(), ueg.end(), ep[0]);

          auto e_grid = e.slice(tensor::range(xs.threshold, e.size()));

          auto rxn_xs = tensor::interp(e_grid, ep, xsp, 0.0, xs.value.back());
          xs.value = vector<double>(rxn_xs.cbegin(), rxn_xs.cend());
        }
      }

      grid.erase(grid.begin(), grid.end());
      nuc->create_ue_derived(nuc->prompt_photons_.get(), nuc->delayed_photons_.get(), ueg);
    }


    double mem_size = (double)(ueg.size()*num_temps) * 12.0 / 1073741824.0;
    return mem_size;
  }
} // namespace openmc
