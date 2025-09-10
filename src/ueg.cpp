#include<fstream>

#include<fmt/core.h>

#include "openmc/ueg.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/search.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xnoalias.hpp"

#include "xtensor/xnpy.hpp"

#include "openmc/reaction.h"

namespace openmc {
  namespace data {
    std::shared_ptr<Nuclide::EnergyGrid> union_e_grid;
  } // namespace data

  void create_union_energy_grid() {
    data::union_e_grid = std::make_shared<Nuclide::EnergyGrid>();
    int neutron = static_cast<int>(ParticleType::neutron);
    double E_min = data::energy_min[neutron];
    double E_max = data::energy_max[neutron];

    settings::energy_cutoff[0] = std::max(E_min, settings::energy_cutoff[0]);

    //imp_e_grid will contain energy points that should not be thinned (URR and Sab energies)
    vector<double> imp_e_grid {E_min, E_max};
    vector<double>& ueg = data::union_e_grid->energy;
    vector<int>& ueg_index = data::union_e_grid->grid_index;

    //Get energies from all nuclides
    for (const auto& nuclide : data::nuclides) {
      for (int t = 0; t < nuclide->kTs_.size(); t++) {
         vector<double>& energies = nuclide->grid_[t].energy;
         ueg.insert(ueg.end(), energies.begin(), energies.end());

        //Add URR energies to important energy grid
         if (nuclide->urr_present_) {
          const auto& urr_energies = nuclide->urr_data_[t].energy_;
          imp_e_grid.insert(imp_e_grid.end(), urr_energies.begin(), urr_energies.end());
        }
      
        //Add threshold energies to important energy grid
        for (auto& rxn : nuclide->reactions_)
          imp_e_grid.insert(imp_e_grid.end(), energies[rxn->xs_[t].threshold]);
      }
    }
    //Sort ueg energy and thin redundant points according to the thinning cutoff parameter
    std::sort(ueg.begin(), ueg.end());
    thin_union_energy_grid();

    auto min_it = ueg.begin();
    auto max_it = --ueg.end();

    while (*min_it < E_min) min_it++;
    while (*max_it > E_max) max_it--;

    ueg.erase(ueg.begin(), min_it + 1);
    ueg.erase(max_it - 1, ueg.end());
    
    ueg.insert(ueg.end(), imp_e_grid.begin(), imp_e_grid.end());
    std::sort(ueg.begin(), ueg.end());
    ueg.erase(std::unique(ueg.begin(), ueg.end()), ueg.end());

    //Generate logarithmic bin indices for double indexing
    //Identical algorithm to Nuclide::init_grid
    int M = settings::n_log_bins;

    double pseudo_spacing = std::log(E_max / E_min);
    auto umesh = xt::linspace(0.0, pseudo_spacing, M + 1);

    int ueg_size = ueg.size();
    ueg_index.resize(M + 1);
    int j = 0;

    for (int k = 0; k <= M; ++k) {
      while (std::log(ueg[j + 1] / E_min) <= umesh[k]) {
       if (j + 2 == ueg_size) break;
       ++j;
      }
      ueg_index[k] = j;
    }
    create_union_energy_xs();
  }

  void thin_union_energy_grid() {
    vector<double>& ueg = data::union_e_grid->energy;
    double tau = settings::ue_grid_cutoff;
    int grid_size = 0;
    for (int i = 0; i < ueg.size() - 1; i++) {
      double current_e = ueg[i];
      double next_e = ueg[i + 1];
      //Replace energy point with an average if its neighbor is
      //within relative difference of the cutoff
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
    auto e = xt::adapt(ueg.energy);
    //xt::dump_npy("orig.npy", xt::adapt(data::nuclides[0]->grid_[0].energy));
    //xt::dump_npy("orig_xs.npy", data::nuclides[0]->xs_[0]);
    //Iterate through all nuclides to update XS
    for (auto & nuclide : data::nuclides) {
      write_message("Processing {}", nuclide->name_);
      auto & grid = nuclide->grid_;
      //Interpolate XS for each nuclide temperature and reaction
      for (auto& rxn : nuclide->reactions_) {
        for (int t = 0; t < nuclide->kTs_.size(); t++) {
          auto & xs = rxn->xs_[t];
          auto ep = xt::view(xt::adapt(grid[t].energy), xt::range(xs.threshold, grid[t].energy.size()));
          
          if (xs.threshold != 0)
            xs.threshold = lower_bound_index(ueg.energy.begin(), ueg.energy.end(), ep[0]);

          auto xsp = xt::view(xt::adapt(xs.value), xt::all());
          auto rxn_xs = xt::interp(xt::view(e, xt::range(xs.threshold, e.size())), ep, xsp, 0.0, xs.value.back());
          
          xs.value = vector<double>(rxn_xs.begin(), rxn_xs.end());
        }
      }
      for (int t = 0; t < nuclide->kTs_.size(); t++) grid[t] = ueg;
      auto temp = xt::interp(xt::adapt(ueg.energy), xt::adapt(nuclide->energy_0K_), xt::adapt(nuclide->elastic_0K_));
      nuclide->energy_0K_ = ueg.energy;
      nuclide->elastic_0K_ = vector<double>(temp.begin(), temp.end()); 
      nuclide->create_ue_derived(nuclide->prompt_photons_.get(), nuclide->delayed_photons_.get());
    }
  //xt::dump_npy("erg.npy", xt::adapt(data::union_e_grid->energy));
  //xt::dump_npy("xs.npy", data::nuclides[0]->xs_[0]);
  }
} // namespace openmc
