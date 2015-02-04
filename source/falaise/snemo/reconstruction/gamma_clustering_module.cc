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
#include <geomtools/manager.h>

// This project:
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/processing/services.h>
#include <snemo/geometry/locator_plugin.h>
#include <snemo/geometry/calo_locator.h>
#include <snemo/geometry/xcalo_locator.h>
#include <snemo/geometry/gveto_locator.h>

namespace snemo {

  namespace reconstruction {

    // Registration instantiation macro :
    DPP_MODULE_REGISTRATION_IMPLEMENT(gamma_clustering_module,
                                      "snemo::reconstruction::gamma_clustering_module");

    const geomtools::manager & gamma_clustering_module::get_geometry_manager() const
    {
      return *_geometry_manager_;
    }

    void gamma_clustering_module::set_geometry_manager(const geomtools::manager & gmgr_)
    {
      DT_THROW_IF (is_initialized(), std::logic_error,
                   "Module '" << get_name() << "' is already initialized ! ");
      _geometry_manager_ = &gmgr_;

      // Check setup label:
      const std::string & setup_label = _geometry_manager_->get_setup_label();
      DT_THROW_IF (setup_label != "snemo::demonstrator" &&
                   setup_label != "snemo::tracker_commissioning",
                   std::logic_error,
                   "Setup label '" << setup_label << "' is not supported !");
      return;
    }

    void gamma_clustering_module::_set_defaults()
    {
      _geometry_manager_ = 0;
      _PTD_label_ = snemo::datamodel::data_info::default_tracker_clustering_data_label();
      return;
    }

    // Initialization :
    void gamma_clustering_module::initialize(const datatools::properties  & setup_,
                                             datatools::service_manager   & service_manager_,
                                             dpp::module_handle_dict_type & /* module_dict_ */)
    {
      DT_THROW_IF (is_initialized(),
                   std::logic_error,
                   "Module '" << get_name() << "' is already initialized ! ");

      dpp::base_module::_common_initialize(setup_);

      if (setup_.has_key("PTD_label")) {
        _PTD_label_ = setup_.fetch_string("PTD_label");
      }

      // Geometry manager :
      if (_geometry_manager_ == 0) {
        std::string geo_label = snemo::processing::service_info::default_geometry_service_label();
        if (setup_.has_key("Geo_label")) {
          geo_label = setup_.fetch_string("Geo_label");
        }
        DT_THROW_IF (geo_label.empty(), std::logic_error,
                     "Module '" << get_name() << "' has no valid '" << "Geo_label" << "' property !");
        DT_THROW_IF (! service_manager_.has(geo_label) ||
                     ! service_manager_.is_a<geomtools::geometry_service>(geo_label),
                     std::logic_error,
                     "Module '" << get_name() << "' has no '" << geo_label << "' service !");
        geomtools::geometry_service & Geo
          = service_manager_.get<geomtools::geometry_service>(geo_label);
        set_geometry_manager(Geo.get_geom_manager());


      // Get geometry locator plugin
        const geomtools::manager & geo_mgr = Geo.get_geom_manager();
        std::string locator_plugin_name;
        if (setup_.has_key ("locator_plugin_name"))
          {
            locator_plugin_name = setup_.fetch_string ("locator_plugin_name");
          }
        else
          {
            // If no locator plugin name is set, then search for the first one
            const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins ();
            for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin ();
                 ip != plugins.end ();
                 ++ip) {
              const std::string & plugin_name = ip->first;
              if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin> (plugin_name)) {
                DT_LOG_DEBUG (get_logging_priority (), "Find locator plugin with name = " << plugin_name);
                locator_plugin_name = plugin_name;
                break;
              }
            }
          }
        // Access to a given plugin by name and type :
        DT_THROW_IF (! geo_mgr.has_plugin (locator_plugin_name) ||
                     ! geo_mgr.is_plugin_a<snemo::geometry::locator_plugin> (locator_plugin_name),
                     std::logic_error,
                     "Found no locator plugin named '" << locator_plugin_name << "'");
        _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin> (locator_plugin_name);

      }
      _set_initialized(true);
      return;
    }

    void gamma_clustering_module::reset()
    {
      DT_THROW_IF (! is_initialized(),
                   std::logic_error,
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
      DT_THROW_IF (! is_initialized(), std::logic_error,
                   "Module '" << get_name() << "' is not initialized !");

      const snemo::datamodel::calibrated_calorimeter_hit::collection_type * the_calos = 0;

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
        the_calos = &ptr_particle_track_data->get_non_associated_calorimeters();
      }
      snemo::datamodel::particle_track_data & the_particle_track_data = *ptr_particle_track_data;

      /********************
       * Process the data *
       ********************/


      // Sanity check
      if (! the_calos) {
        DT_LOG_WARNING(get_logging_priority(), "No calorimeter hits to be processed !");
        return dpp::base_module::PROCESS_ERROR;
      }

      gamma_dict_type clustered_gammas;

      // Main processing method :
      _process(the_particle_track_data, clustered_gammas);

      return dpp::base_module::PROCESS_SUCCESS;
    }

    void gamma_clustering_module::_get_new_neighbours(geomtools::geom_id gid,
                                                      const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & cch,
                                                      std::vector<geomtools::geom_id>  & ccl,
                                                      std::vector<geomtools::geom_id>  & a_cluster)
    {
      if(std::find(ccl.begin(), ccl.end(),gid)==ccl.end())
        ccl.push_back(gid);
      else
        return;

      std::vector<geomtools::geom_id>  the_neighbours;
      std::vector<geomtools::geom_id>  the_calib_neighbours;

      const snemo::geometry::calo_locator & calo_locator
        = _locator_plugin_->get_calo_locator();
      const snemo::geometry::xcalo_locator & xcalo_locator
        = _locator_plugin_->get_xcalo_locator();
      const snemo::geometry::gveto_locator & gveto_locator
        = _locator_plugin_->get_gveto_locator();

      if (calo_locator.is_calo_block_in_current_module(gid))
        calo_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

      if (xcalo_locator.is_calo_block_in_current_module(gid))
        xcalo_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

      if (gveto_locator.is_calo_block_in_current_module(gid))
        gveto_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

      // for (auto ineighbour : the_neighbours)
      for (std::vector<geomtools::geom_id>::const_iterator ineighbour = the_neighbours.begin();
           ineighbour != the_neighbours.end(); ++ineighbour) {

        // if (std::find_if(cch.begin(), cch.end(), [ineighbour] (auto icalo)
        //                  {return ineighbour == icalo.get().get_geom_id();}) != cch.end())
        //find if the eight ineighbours belong also to cch
        for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
               icalo = cch.begin(); icalo != cch.end(); ++icalo) {
          const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = icalo->get();

          if (*ineighbour == a_calo_hit.get_geom_id())
            if(std::find(ccl.begin(), ccl.end(), *ineighbour) == ccl.end()) {
              the_calib_neighbours.push_back(*ineighbour);
              ccl.push_back(*ineighbour);
              a_cluster.push_back(*ineighbour);
            }
        }
      }
      // for(auto i_calib_neighbour : the_calib_neighbours)
      for (std::vector<geomtools::geom_id>::const_iterator i_calib_neighbour = the_calib_neighbours.begin();
           i_calib_neighbour != the_calib_neighbours.end(); ++i_calib_neighbour)
        _get_new_neighbours(*i_calib_neighbour, cch, ccl, a_cluster);
    }


    void gamma_clustering_module::_process(snemo::datamodel::particle_track_data  & ptd_, gamma_dict_type & clustered_gammas_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      const snemo::datamodel::calibrated_calorimeter_hit::collection_type & cch = ptd_.get_non_associated_calorimeters();

      std::vector<geomtools::geom_id>  ccl;

      size_t number_of_clusters = 0;

      std::vector<std::vector<geomtools::geom_id> >  the_reconstructed_clusters;

      // for (auto icalo : cch) {

      for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
             icalo = cch.begin(); icalo != cch.end(); ++icalo) {

        // const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = icalo->get();

        const geomtools::geom_id & gid = icalo->get().get_geom_id();

        std::vector<geomtools::geom_id> a_cluster;
        a_cluster.push_back(gid);

        if(std::find(ccl.begin(), ccl.end(),gid)!=ccl.end())
          continue;

        _get_new_neighbours(gid, cch, ccl, a_cluster);

        the_reconstructed_clusters.push_back(a_cluster);

        number_of_clusters++;
      }

      std::vector<std::map<double, geomtools::geom_id> >  the_ordered_reconstructed_clusters;

      // for(auto icluster : the_reconstructed_clusters)
      for(std::vector<std::vector<geomtools::geom_id> >::const_iterator icluster = the_reconstructed_clusters.begin();
          icluster != the_reconstructed_clusters.end(); ++icluster)
        {
          std::map<double,geomtools::geom_id> a_cluster;

          // for(auto igid : icluster)
          //   for(auto ihit : cch)
          for (std::vector<geomtools::geom_id>::const_iterator igid = icluster->begin();
               igid != icluster->end(); ++igid)
            for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
                   ihit = cch.begin(); ihit != cch.end(); ++ihit)
              if(*igid == ihit->get().get_geom_id())
                a_cluster.insert( std::pair<double,geomtools::geom_id >(ihit->get().get_time(),*igid) );

          the_ordered_reconstructed_clusters.push_back(a_cluster);
        }

      std::sort(the_ordered_reconstructed_clusters.begin(), the_ordered_reconstructed_clusters.end());

      int track_id = 0;

      //      for(auto icluster : the_ordered_reconstructed_clusters)
      for(std::vector<std::map<double, geomtools::geom_id> >::const_iterator icluster = the_ordered_reconstructed_clusters.begin();
          icluster != the_ordered_reconstructed_clusters.end(); ++icluster)
        {
          track_id++;

          if(icluster->size() < 2)
            {
              // for (auto ipair : icluster)
              for(std::map<double, geomtools::geom_id>::const_iterator ipair = icluster->begin();
                  ipair != icluster->end(); ++ipair)
              clustered_gammas_[track_id].insert(ipair->second);
              continue;
            }

          double t0 = 0; // not ideal
          double t1 = 0;

          // for (auto ipair : icluster)
          for(std::map<double, geomtools::geom_id>::const_iterator ipair = icluster->begin();
              ipair != icluster->end(); ++ipair)
            {
              t0 = t1;
              t1 = ipair->first;
              // std::cout << " " << ipair.first  << "   " << std::next(&ipair)->first << std::endl;
              // std::cout << " " << t0  << "   " << t1 << std::endl;

              if(t0!=0 && t1!=0 && t1-t0 > 2.5 /*ns*/)
                {
                  number_of_clusters++;
                  track_id++;
                }

              clustered_gammas_[track_id].insert(ipair->second);
            }
        }

      std::cout << "Number of clusters : " << number_of_clusters << std::endl;
      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
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
