/// \file falaise/snemo/reconstruction/gamma_clustering_module.cc

// Ourselves:
#include <snemo/reconstruction/gamma_clustering_module.h>

// Standard library:
#include <stdexcept>
#include <sstream>

// Third party:
// - Bayeux/datatools:
#include <datatools/service_manager.h>
// - Bayeux/geomtools:
#include <geomtools/geometry_service.h>

// This project:
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/processing/services.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>

namespace snemo {

  namespace reconstruction {

    // Registration instantiation macro :
    DPP_MODULE_REGISTRATION_IMPLEMENT(gamma_clustering_module,
                                      "snemo::reconstruction::gamma_clustering_module");

    void gamma_clustering_module::_set_defaults()
    {
      _locator_plugin_ = 0;
      _PTD_label_ = snemo::datamodel::data_info::default_tracker_clustering_data_label();
      return;
    }

    // Initialization :
    void gamma_clustering_module::initialize(const datatools::properties  & setup_,
                                             datatools::service_manager   & service_manager_,
                                             dpp::module_handle_dict_type & /* module_dict_ */)
    {
      DT_THROW_IF(is_initialized(),
                  std::logic_error,
                  "Module '" << get_name() << "' is already initialized ! ");

      dpp::base_module::_common_initialize(setup_);

      if (setup_.has_key("PTD_label")) {
        _PTD_label_ = setup_.fetch_string("PTD_label");
      }

      // Geometry manager :
      std::string geo_label = snemo::processing::service_info::default_geometry_service_label();
      if (setup_.has_key("Geo_label")) {
        geo_label = setup_.fetch_string("Geo_label");
      }
      DT_THROW_IF(geo_label.empty(), std::logic_error,
                  "Module '" << get_name() << "' has no valid '" << "Geo_label" << "' property !");
      DT_THROW_IF(! service_manager_.has(geo_label) ||
                  ! service_manager_.is_a<geomtools::geometry_service>(geo_label),
                  std::logic_error,
                  "Module '" << get_name() << "' has no '" << geo_label << "' service !");
      geomtools::geometry_service & Geo
        = service_manager_.get<geomtools::geometry_service>(geo_label);

      // Get geometry locator plugin
      const geomtools::manager & geo_mgr = Geo.get_geom_manager();
      std::string locator_plugin_name;
      if (setup_.has_key ("locator_plugin_name")) {
        locator_plugin_name = setup_.fetch_string ("locator_plugin_name");
      } else {
        // If no locator plugin name is set, then search for the first one
        const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins ();
        for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
             ip != plugins.end();
             ++ip) {
          const std::string & plugin_name = ip->first;
          if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin> (plugin_name)) {
            DT_LOG_DEBUG (get_logging_priority(), "Find locator plugin with name = " << plugin_name);
            locator_plugin_name = plugin_name;
            break;
          }
        }
      }
      // Access to a given plugin by name and type :
      DT_THROW_IF(! geo_mgr.has_plugin(locator_plugin_name) ||
                  ! geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name),
                  std::logic_error,
                  "Found no locator plugin named '" << locator_plugin_name << "'");
      _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);

      _set_initialized(true);
      return;
    }

    void gamma_clustering_module::reset()
    {
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is not initialized !");
      _set_initialized(false);
      _set_defaults();
      return;
    }

    // Constructor :
    gamma_clustering_module::gamma_clustering_module(datatools::logger::priority logging_priority_)
      : dpp::base_module(logging_priority_)
    {
      _set_defaults();
      return;
    }

    // Destructor :
    gamma_clustering_module::~gamma_clustering_module()
    {
      if (is_initialized()) gamma_clustering_module::reset();
      return;
    }

    // Processing :
    dpp::base_module::process_status gamma_clustering_module::process(datatools::things & data_record_)
    {
      DT_THROW_IF(! is_initialized(), std::logic_error,
                  "Module '" << get_name() << "' is not initialized !");

      /*********************************
       * Check particle track data     *
       *********************************/
      snemo::datamodel::particle_track_data * ptr_particle_track_data = 0;
      if (! data_record_.has(_PTD_label_)) {
        ptr_particle_track_data
          = &(data_record_.add<snemo::datamodel::particle_track_data>(_PTD_label_));
      } else {
        ptr_particle_track_data
          = &(data_record_.grab<snemo::datamodel::particle_track_data>(_PTD_label_));
      }
      snemo::datamodel::particle_track_data & ptd = *ptr_particle_track_data;

      /********************
       * Process the data *
       ********************/

      // Main processing method :
      _process(ptd);

      return dpp::base_module::PROCESS_SUCCESS;
    }

    void gamma_clustering_module::_process(snemo::datamodel::particle_track_data  & ptd_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      const snemo::datamodel::calibrated_calorimeter_hit::collection_type & cch
        = ptd_.get_non_associated_calorimeters();

      // Registered calorimeter hits (to be skipped)
      gid_list_type registered_calos;

      // Getting gamma clusters
      std::vector<cluster_type> the_reconstructed_clusters;
      for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
             icalo = cch.begin(); icalo != cch.end(); ++icalo) {
        const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = icalo->get();

        const geomtools::geom_id & a_gid = a_calo_hit.get_geom_id();
        const double a_time              = a_calo_hit.get_time();

        // If already clustered then skip it
        if (std::find(registered_calos.begin(), registered_calos.end(), a_gid)
            != registered_calos.end())
          continue;

        {
          // Insert empty cluster
          DT_LOG_TRACE(get_logging_priority(), "Insert new gamma cluster (#"
                       << the_reconstructed_clusters.size() << ")");
          cluster_type dummy;
          the_reconstructed_clusters.push_back(dummy);
        }
        cluster_type & a_cluster = the_reconstructed_clusters.back();
        a_cluster.insert(std::make_pair(a_time, *icalo));

        // Get neighbours given the current geom id
        _get_new_neighbours(*icalo, cch, a_cluster, registered_calos);
      }

      if (get_logging_priority() >= datatools::logger::PRIO_TRACE) {
        for (size_t i = 0; i < the_reconstructed_clusters.size(); ++i) {
          const cluster_type & a_cluster = the_reconstructed_clusters.at(i);
          DT_LOG_TRACE(get_logging_priority(), "New gamma cluster #" << i
                       << " (" << a_cluster.size() << " associated calorimeters)");
          for (cluster_type::const_iterator j = a_cluster.begin(); j != a_cluster.end(); ++j) {
            const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = j->second.get();
            a_calo_hit.tree_dump();
          }
        }
      }

      for (size_t i = 0; i < the_reconstructed_clusters.size(); ++i) {
        DT_LOG_TRACE(get_logging_priority(), "Adding a new clustered gamma");
        snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
        hPT.grab().set_track_id(ptd_.get_number_of_particles());
        hPT.grab().set_charge(snemo::datamodel::particle_track::neutral);
        ptd_.add_particle(hPT);

        const cluster_type & a_cluster = the_reconstructed_clusters.at(i);
        for (cluster_type::const_iterator j = a_cluster.begin(); j != a_cluster.end(); ++j) {
          const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = j->second.get();;
          hPT.grab().grab_associated_calorimeter_hits().push_back(j->second);

          const geomtools::geom_id & a_gid = a_calo_hit.get_geom_id();
          // Build vertex
          snemo::datamodel::particle_track::handle_spot hBS(new geomtools::blur_spot);
          hPT.grab().grab_vertices().push_back(hBS);
          geomtools::blur_spot & spot = hBS.grab();
          spot.set_hit_id(a_calo_hit.get_hit_id());
          spot.set_geom_id(a_gid);

          const snemo::geometry::calo_locator & calo_locator   = _locator_plugin_->get_calo_locator();
          const snemo::geometry::xcalo_locator & xcalo_locator = _locator_plugin_->get_xcalo_locator();
          const snemo::geometry::gveto_locator & gveto_locator = _locator_plugin_->get_gveto_locator();

          geomtools::vector_3d position;
          std::string label;
          if (calo_locator.is_calo_block_in_current_module(a_gid)) {
            calo_locator.get_block_position(a_gid, position);
            label = snemo::datamodel::particle_track::vertex_on_main_calorimeter_label();
          } else if (xcalo_locator.is_calo_block_in_current_module(a_gid)) {
            xcalo_locator.get_block_position(a_gid, position);
            label = snemo::datamodel::particle_track::vertex_on_x_calorimeter_label();
          } else if (gveto_locator.is_calo_block_in_current_module(a_gid)) {
            gveto_locator.get_block_position(a_gid, position);
            label = snemo::datamodel::particle_track::vertex_on_gamma_veto_label();
          } else {
            DT_THROW_IF(true, std::logic_error,
                        "Current geom id '" << a_gid << "' does not match any scintillator block !");
          }
          spot.grab_auxiliaries().store(snemo::datamodel::particle_track::vertex_type_key(),
                                        label);
          spot.set_blur_dimension(geomtools::blur_spot::dimension_three);
          spot.set_position(position);
        }
      }

      DT_LOG_DEBUG(get_logging_priority(), "Number of clusters : " << the_reconstructed_clusters.size());
      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return;
    }

    void gamma_clustering_module::_get_new_neighbours(const snemo::datamodel::calibrated_data::calorimeter_hit_handle_type & hit_,
                                                      const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & hits_,
                                                      cluster_type & cluster_,
                                                      gid_list_type & registered_calos_) const
    {
      const snemo::datamodel::calibrated_calorimeter_hit & a_hit = hit_.get();
      const geomtools::geom_id & a_gid = a_hit.get_geom_id();

      // If already clustered then skip it
      if (std::find(registered_calos_.begin(), registered_calos_.end(), a_gid)
          != registered_calos_.end())
        return;

      // Store the current calorimeter as registered one
      registered_calos_.push_back(a_gid);
      if (get_logging_priority() >= datatools::logger::PRIO_TRACE) {
        for (gid_list_type::const_iterator i = registered_calos_.begin();
             i != registered_calos_.end(); ++i) {
          DT_LOG_TRACE(get_logging_priority(), "Registered geom id = " << *i);
        }
      }

      gid_list_type the_neighbours;
      const snemo::geometry::calo_locator & calo_locator   = _locator_plugin_->get_calo_locator();
      const snemo::geometry::xcalo_locator & xcalo_locator = _locator_plugin_->get_xcalo_locator();
      const snemo::geometry::gveto_locator & gveto_locator = _locator_plugin_->get_gveto_locator();

      if (calo_locator.is_calo_block_in_current_module(a_gid)) {
        calo_locator.get_neighbours_ids(a_gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);
      } else if (xcalo_locator.is_calo_block_in_current_module(a_gid)) {
        xcalo_locator.get_neighbours_ids(a_gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);
      } else if (gveto_locator.is_calo_block_in_current_module(a_gid)) {
        gveto_locator.get_neighbours_ids(a_gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);
      } else {
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << a_gid << "' does not match any scintillator block !");
      }

      for (gid_list_type::const_iterator ineighbour = the_neighbours.begin();
           ineighbour != the_neighbours.end(); ++ineighbour) {
        const geomtools::geom_id & a_gid = *ineighbour;

        // Find of the eight neighbours belong to calibrated calo. hits
        geomtools::base_hit::has_geom_id_predicate hit_pred(a_gid);
        datatools::mother_to_daughter_predicate<geomtools::base_hit,
                                                snemo::datamodel::calibrated_calorimeter_hit> pred_M2D(hit_pred);
        datatools::handle_predicate<snemo::datamodel::calibrated_calorimeter_hit> pred_via_handle(pred_M2D);
        snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
          found = std::find_if(hits_.begin(), hits_.end(), pred_via_handle);

        if (found != hits_.end()) {
          if (std::find(registered_calos_.begin(), registered_calos_.end(), a_gid)
              == registered_calos_.end()) {
            registered_calos_.push_back(a_gid);
            cluster_.insert(std::make_pair(found->get().get_time(), *found));
            _get_new_neighbours(*found, hits_, cluster_, registered_calos_);
          }
        }
        // 2015/02/06 SC : using lambda function in C++14 (really easier)
        // if (std::find_if(cch.begin(), cch.end(), [ineighbour] (auto icalo)
        //                  {return ineighbour == icalo.get().get_geom_id();}) != cch.end())
      }
      return;
    }

  } // end of namespace reconstruction

} // end of namespace snemo

/* OCD support */
#include <datatools/object_configuration_description.h>
DOCD_CLASS_IMPLEMENT_LOAD_BEGIN(snemo::reconstruction::gamma_clustering_module, ocd_)
{
  ocd_.set_class_name("snemo::reconstruction::gamma_clustering_module");
  ocd_.set_class_description("A module that performs the gamma clustering using the Gamma_Clustering algorithms");
  ocd_.set_class_library("Falaise_GammaClustering");
  ocd_.set_class_documentation("This module uses the Gamma Clustering algorithms for.   \n"
                               "gamma involved in non associated calorimeter hits. See also OCD   \n"
                               "support for the ``snemo::reconstruction::gamma_clustering_driver`` class. \n");

  dpp::base_module::common_ocd(ocd_);

  {
    // Description of the 'PTD_label' configuration property :
    datatools::configuration_property_description & cpd
      = ocd_.add_property_info();
    cpd.set_name_pattern("PTD_label")
      .set_terse_description("The label/name of the 'particle track data' bank")
      .set_traits(datatools::TYPE_STRING)
      .set_mandatory(false)
      .set_long_description("This is the name of the bank to be used as    \n"
                            "the source of calorimeter hits and reconstructed vertices. \n"
                            )
      .set_default_value_string(snemo::datamodel::data_info::default_particle_track_data_label())
      .add_example("Use an alternative name for the 'particle track data' bank:: \n"
                   "                                  \n"
                   "  PTD_label : string = \"PTD2\"   \n"
                   "                                  \n"
                   )
      ;
  }

  {
    // Description of the 'Geo_label' configuration property :
    datatools::configuration_property_description & cpd
      = ocd_.add_property_info();
    cpd.set_name_pattern("Geo_label")
      .set_terse_description("The label/name of the geometry service")
      .set_traits(datatools::TYPE_STRING)
      .set_mandatory(false)
      .set_long_description("This is the name of the service to be used as the \n"
                            "geometry service.                                 \n"
                            "This property is only used if no geometry manager \n"
                            "as been provided to the module.                   \n"
                            )
      .set_default_value_string(snemo::processing::service_info::default_geometry_service_label())
      .add_example("Use an alternative name for the geometry service:: \n"
                   "                                     \n"
                   "  Geo_label : string = \"geometry2\" \n"
                   "                                     \n"
                   )
      ;
  }

  // Additionnal configuration hints :
  ocd_.set_configuration_hints("Here is a full configuration example in the      \n"
                               "``datatools::properties`` ASCII format::         \n"
                               "                                         \n"
                               "  PTD_label : string = \"PTD\"           \n"
                               "  Geo_label : string = \"geometry\"      \n"
                               "                                         \n"
                               "Additional specific parameters are used to configure         \n"
                               "the embedded ``Gamma_Clustering`` driver itself; see the OCD support \n"
                               "of the ``snemo::reconstruction::gamma_clustering_driver`` class.     \n"
                               );

  ocd_.set_validation_support(true);
  ocd_.lock();
  return;
}
DOCD_CLASS_IMPLEMENT_LOAD_END() // Closing macro for implementation
DOCD_CLASS_SYSTEM_REGISTRATION(snemo::reconstruction::gamma_clustering_module,
                               "snemo::reconstruction::gamma_clustering_module")
