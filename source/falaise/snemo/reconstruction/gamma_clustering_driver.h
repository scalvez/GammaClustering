/** \file falaise/snemo/reconstruction/gamma_clustering_driver.h
 * Author(s)     : Xavier Garrido <garrido@lal.in2p3.fr>
 * Creation date : 2012-10-07
 * Last modified : 2014-02-09
 *
 * Copyright (C) 2012-2014 Xavier Garrido <garrido@lal.in2p3.fr>
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
 * Description:
 *
 *   A driver class that wraps the GammaClustering algorithm.
 *
 * History:
 *
 */

#ifndef FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_DRIVER_H
#define FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_DRIVER_H 1

// Standard library:
#include <string>
#include <list>
#include <map>
#include <vector>

// Third party:
// - Boost:
#include <boost/scoped_ptr.hpp>

// - Bayeux/datatools:
#include <datatools/logger.h>
#include <datatools/properties.h>

// This project:
#include <falaise/snemo/datamodels/calibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

namespace geomtools {
  class manager;
}

namespace snemo {

  namespace datamodel {
    class particle_track_data;
  }

  namespace geometry {
    class locator_plugin;
  }

  namespace reconstruction {

    /// Driver for the gamma clustering algorithms
    class gamma_clustering_driver
    {
    public:
      /// Typedef for list of geom ids
      typedef std::vector<geomtools::geom_id> gid_list_type;

      /// Typedef for time ordered calorimeter hits aka. gamma cluster
      typedef std::map<double, const snemo::datamodel::calibrated_data::calorimeter_hit_handle_type> cluster_type;

      /// Typedef for collection of clusters
      typedef std::vector<cluster_type> cluster_collection_type;

      static const std::string & gamma_clustering_id();

      /// Initialization flag
      void set_initialized(const bool initialized_);

      /// Getting initialization flag
      bool is_initialized() const;

      /// Setting logging priority
      void set_logging_priority(const datatools::logger::priority priority_);

      /// Getting logging priority
      datatools::logger::priority get_logging_priority() const;

      /// Check the geometry manager
      bool has_geometry_manager() const;

      /// Address the geometry manager
      void set_geometry_manager(const geomtools::manager & gmgr_);

      /// Return a non-mutable reference to the geometry manager
      const geomtools::manager & get_geometry_manager() const;

      /// Constructor
      gamma_clustering_driver();

      /// Destructor
      virtual ~gamma_clustering_driver();

      /// Initialize the gamma tracker through configuration properties
      virtual void initialize(const datatools::properties & setup_);

      /// Reset the clusterizer
      virtual void reset();

      /// Data record processing
      virtual int process(snemo::datamodel::particle_track_data & ptd_);

    protected:

      /// Special method to process and generate particle track data
      void _process(snemo::datamodel::particle_track_data & ptd_);

      /// Give default values to specific class members.
      void _set_defaults ();

      /// Get calorimeter neighbours given teh current calorimeter hit
      virtual void _get_geometrical_neighbours(const snemo::datamodel::calibrated_calorimeter_hit & hit_,
                                               const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & hits_,
                                               cluster_type & cluster_,
                                               gid_list_type & registered_calos_) const;

      /// Split calorimeter cluster given a cluster time range value
      virtual void _get_time_neighbours(cluster_type & cluster_, cluster_collection_type & clusters_) const;

      /// Associate clusters given Time-Of-Flight calculation
      virtual void _get_tof_association(const cluster_collection_type & from_clusters_,
                                        cluster_collection_type & to_clusters_) const;

      virtual double _get_tof_probability(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                                          const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit) const;

      virtual bool _are_on_same_wall(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                                     const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit) const;

    private:

      bool _initialized_;                                       //!< Initialize flag
      datatools::logger::priority _logging_priority_;           //!< Logging priority
      const geomtools::manager * _geometry_manager_;            //!< The SuperNEMO geometry manager
      const snemo::geometry::locator_plugin * _locator_plugin_; //!< The SuperNEMO locator plugin
      double _cluster_time_range_;     //!< The time condition for clustering
      std::string _cluster_grid_mask_; //!< The spatial condition for clustering
      double _min_prob_;     //!< The minimal probability required between clusters
      double _sigma_good_calo_;     //!< The minimal probability required between clusters
      std::string _PTD_label_;                       //!< The label of the input/output  data bank

      datatools::properties _gc_setup_;                         //!< The Gamma Clustering parameters
      // for members
    };

  }  // end of namespace reconstruction

}  // end of namespace snemo


// Declare the OCD interface of the module
#include <datatools/ocd_macros.h>
DOCD_CLASS_DECLARATION(snemo::reconstruction::gamma_clustering_driver)

#endif // FALAISE_GAMMA_CLUSTERING_PLUGIN_SNEMO_RECONSTRUCTION_GAMMA_CLUSTERING_DRIVER_H

/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** tab-width: 2 --
** End: --
*/
