// -*- mode: c++ ; -*-
/** \file falaise/snemo/reconstruction/gamma_clustering_module.h
 * Author(s) :    Xavier Garrido <garrido@lal.in2p3.fr>
 * Creation date: 2012-10-07
 * Last modified: 2014-02-28
 *
 * Copyright (C) 2011-2014 Xavier Garrido <garrido@lal.in2p3.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 *
 * Description:
 *
 *   Module for gamma clustering
 *
 * History:
 *
 */

#ifndef FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_MODULE_H
#define FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_MODULE_H 1

// Third party:
// - Boost:
#include <boost/scoped_ptr.hpp>
// - Bayuex/dpp :
#include <dpp/base_module.h>
// - Bayeux/geomtools:
#include <geomtools/geom_id.h>

// This project:
#include <falaise/snemo/datamodels/calibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

namespace snemo {

  namespace datamodel {
    class particle_track_data;
  }

  namespace geometry {
    class locator_plugin;
  }

  namespace reconstruction {

    /// \brief The data processing module for the gamma clustering
    class gamma_clustering_module : public dpp::base_module
    {

    public:
      /// Typedef for list of geom ids
      typedef std::vector<geomtools::geom_id> gid_list_type;

      /// Typedef for time ordered calorimeter hits aka. gamma cluster
      typedef std::map<double, const snemo::datamodel::calibrated_data::calorimeter_hit_handle_type> cluster_type;

      /// Typedef for collection of clusters
      typedef std::vector<cluster_type> cluster_collection_type;

      /// Constructor
      gamma_clustering_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);

      /// Destructor
      virtual ~gamma_clustering_module();

      /// Initialization
      virtual void initialize(const datatools::properties  & setup_,
                              datatools::service_manager   & service_manager_,
                              dpp::module_handle_dict_type & module_dict_);

      /// Reset
      virtual void reset();

      /// Data record processing
      virtual process_status process(datatools::things & data_);

    protected:

      /// Special method to process and generate particle track data
      void _process(snemo::datamodel::particle_track_data & ptd_);

      /// Give default values to specific class members.
      void _set_defaults ();

      /// Get calorimeter neighbours given teh current calorimeter hit
      void _get_geometrical_neighbours(const snemo::datamodel::calibrated_data::calorimeter_hit_handle_type & hit_,
                                       const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & hits_,
                                       cluster_type & cluster_,
                                       gid_list_type & registered_calos_) const;

      void _get_time_neighbours(cluster_type & cluster_, cluster_collection_type & clusters_) const;

      double _get_probability(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                              const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit);

      bool _are_on_same_wall(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                              const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit);

    private:

      std::string _PTD_label_;                                  //!< The label of the input/output  data bank
      const snemo::geometry::locator_plugin * _locator_plugin_; //!< The locator plugin

      double _cluster_time_range_;     //!< The time condition for clustering
      std::string _cluster_grid_mask_; //!< The spatial condition for clustering
      double _min_prob_;     //!< The minimal probability required between clusters
      double _sigma_good_calo_;     //!< The minimal probability required between clusters

      // Macro to automate the registration of the module :
      DPP_MODULE_REGISTRATION_INTERFACE(gamma_clustering_module);
    };

  } // end of namespace reconstruction

} // end of namespace snemo

#include <datatools/ocd_macros.h>

// Declare the OCD interface of the module
DOCD_CLASS_DECLARATION(snemo::reconstruction::gamma_clustering_module)

#endif // FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_MODULE_H
