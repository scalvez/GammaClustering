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

//- GSL:
#include <gsl/gsl_cdf.h>

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
      _PTD_label_ = snemo::datamodel::data_info::default_tracker_clustering_data_label();
      _locator_plugin_ = 0;

      _cluster_time_range_ = 2.5 * CLHEP::ns;
      _cluster_grid_mask_ = "first";
      _min_prob_ = 0.99;
      _sigma_good_calo_ = 1.5 * CLHEP::ns;
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

      std::string key;
      if (setup_.has_key(key = "cluster.time_range")) {
        _cluster_time_range_ = setup_.fetch_real(key);
        if (! setup_.has_explicit_unit(key)) {
          _cluster_time_range_ *= CLHEP::ns;
        }
      }

      if (setup_.has_key(key = "cluster.grid_mask")) {
        _cluster_grid_mask_ = setup_.fetch_string(key);
      }

      if (setup_.has_key(key = "cluster.min_prob")) {
        _min_prob_ = setup_.fetch_real(key);
      }

      if (setup_.has_key(key = "cluster.sigma_good_calo")) {
        _sigma_good_calo_ = setup_.fetch_real(key);
        if (! setup_.has_explicit_unit(key)) {
          _sigma_good_calo_ *= CLHEP::ns;
        }
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

    void gamma_clustering_module::_process(snemo::datamodel::particle_track_data & ptd_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      const snemo::datamodel::calibrated_calorimeter_hit::collection_type & cch
        = ptd_.get_non_associated_calorimeters();

      // Registered calorimeter hits (to be skipped)
      gid_list_type registered_calos;

      // Getting gamma clusters
      cluster_collection_type the_reconstructed_clusters;
      for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
             icalo = cch.begin(); icalo != cch.end(); ++icalo) {
        const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = icalo->get();

        const geomtools::geom_id & a_gid = a_calo_hit.get_geom_id();
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
        a_cluster.insert(std::make_pair(a_calo_hit.get_time(), *icalo));

        // Get geometrical neighbours given the current geom id
        _get_geometrical_neighbours(a_calo_hit, cch, a_cluster, registered_calos);

        // Ensure all calorimeter hits within a cluster are in time
        _get_time_neighbours(a_cluster, the_reconstructed_clusters);
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

      /*****  Associate clusters from TOF callculations  *****/

      // if( GC* )
      //It holds the indices of the two clusters to be later concatenated
      std::map <size_t,size_t> merge_indices;

      for (size_t i = 0; i < the_reconstructed_clusters.size(); ++i) {

        // Dirty fix, otherwise bug
        if(i == the_reconstructed_clusters.size()-1)
          break;

        // The candidate for association to next_cluster
        const cluster_type & a_cluster = the_reconstructed_clusters.at(i);

        // Retrieve the last calorimeter for check in quality
        cluster_type::const_iterator it_head = a_cluster.end();
        --it_head;

        // The algo tries to find a tail to a head

        if(a_cluster.size()>1)   // if a cluster is just one calo, no criterion on its quality
          {
            bool found_good_calo_head = false;

            while(!found_good_calo_head) {
              if(it_head->second.get().get_sigma_time() < _sigma_good_calo_)
                found_good_calo_head = true;

              if(it_head == a_cluster.begin())
                break;
              else
                --it_head;
            }
            // if no good calorimeter has been found, by default take the last in time
            if(!found_good_calo_head) {
              it_head = a_cluster.end();
              --it_head;
            }
          }

        const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit = it_head->second.get();

        // Holds the probability from head->tail and the index of the tail
        std::map <double,size_t> possible_clusters_association;

        for (size_t j = i+1; j < the_reconstructed_clusters.size(); ++j) {

          const cluster_type & next_cluster = the_reconstructed_clusters.at(j);
          cluster_type::const_iterator it_tail = next_cluster.begin();

          if(_are_on_same_wall(head_end_calo_hit,it_tail->second.get()))
            continue;

          if(next_cluster.size()>1)
            {
              bool found_good_calo_tail = false;

              while(!found_good_calo_tail) {
                const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit = it_tail->second.get();

                if(it_tail->second.get().get_sigma_time() < _sigma_good_calo_) {
                  found_good_calo_tail = true;
                  break;
                }

                ++it_tail;

                if(it_tail == next_cluster.end())
                  break;
              }

              if(!found_good_calo_tail)
                it_tail = next_cluster.begin();
            }

          const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit = it_tail->second.get();

          double proba = _get_probability(head_end_calo_hit,tail_begin_calo_hit);

          if(proba>_min_prob_)
            possible_clusters_association.insert(std::pair<double,size_t>(proba,j));  // keep all possible solutions
        }

        if(possible_clusters_association.size()==0)
          continue;

        std::map<double,size_t>::const_iterator it_best_proba = possible_clusters_association.end();
        --it_best_proba;

        bool tail_already_associated = false;

        //the probability distribution is flat so above P=50%, there are as much chances for a 51% pair and a 99% to be the correct pair
        //Here, it arbitrarily chooses the first pair built
        for(std::map<size_t,size_t>::const_iterator it_indices=merge_indices.begin(); it_indices!=merge_indices.end(); ++it_indices)
          if(it_indices->second==it_best_proba->second) {
            tail_already_associated = true;
            break;
          }

        if(!tail_already_associated)
          merge_indices.insert(std::pair<size_t,size_t>(i,it_best_proba->second));
      }

      cluster_collection_type the_reconstructed_gammas;

      std::vector <size_t> cluster_to_be_considered;

      for(size_t i = 0;i<the_reconstructed_clusters.size();++i) // initialize with all the clusters in the event
        cluster_to_be_considered.push_back(i);

      for(std::map<size_t,size_t>::const_iterator i_pair = merge_indices.begin(); i_pair!=merge_indices.end(); ++i_pair)
        {
          size_t i_cluster = i_pair->first;
          size_t i_next_cluster;

          if(std::find(cluster_to_be_considered.begin(),cluster_to_be_considered.end(),i_cluster)==cluster_to_be_considered.end())
            continue;  // skip the cluster if it has already been involved in an association before

          cluster_type new_cluster;

          //Fills a new cluster made of the concatenation of successive clusters
          while(merge_indices.find(i_cluster)!=merge_indices.end())
            {
              i_next_cluster = merge_indices.at(i_cluster);

              if(std::find(cluster_to_be_considered.begin(),cluster_to_be_considered.end(),i_cluster)!=cluster_to_be_considered.end())
                cluster_to_be_considered.erase(std::find(cluster_to_be_considered.begin(),cluster_to_be_considered.end(),i_cluster));
              if(std::find(cluster_to_be_considered.begin(),cluster_to_be_considered.end(),i_next_cluster)!=cluster_to_be_considered.end())
                cluster_to_be_considered.erase(std::find(cluster_to_be_considered.begin(),cluster_to_be_considered.end(),i_next_cluster));

              if(new_cluster.size() == 0)
                new_cluster.insert(the_reconstructed_clusters.at(i_cluster).begin(),the_reconstructed_clusters.at(i_cluster).end());

              new_cluster.insert(the_reconstructed_clusters.at(i_next_cluster).begin(),the_reconstructed_clusters.at(i_next_cluster).end());

              i_cluster = i_next_cluster;
            }

          the_reconstructed_gammas.push_back(new_cluster);
        }

      //Add the remaining isolated clusters
      for(std::vector<size_t>::const_iterator i_solo = cluster_to_be_considered.begin(); i_solo != cluster_to_be_considered.end(); ++i_solo) {
        cluster_type & new_cluster = the_reconstructed_clusters.at(*i_solo);
        the_reconstructed_gammas.push_back(new_cluster);
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

      if (get_logging_priority() >= datatools::logger::PRIO_TRACE) {
        std::cout << "=============================================" << std::endl;
        for (size_t i = 0; i < the_reconstructed_gammas.size(); ++i) {
          const cluster_type & a_gamma = the_reconstructed_gammas.at(i);
          DT_LOG_TRACE(get_logging_priority(), "New gamma cluster #" << i
                       << " (" << a_gamma.size() << " associated calorimeters)");
          for (cluster_type::const_iterator j = a_gamma.begin(); j != a_gamma.end(); ++j) {
            const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = j->second.get();
            a_calo_hit.tree_dump();
          }
        }
      }

      /*********  End of GC*  *********************/

      // Set new particles within 'particle track data' container
      if (ptd_.has_non_associated_calorimeters()) {
        DT_LOG_DEBUG(get_logging_priority(), "Removing non-associated calorimeters");
        ptd_.reset_non_associated_calorimeters();
      }
      // for (size_t i = 0; i < the_reconstructed_clusters.size(); ++i) {
      for (size_t i = 0; i < the_reconstructed_gammas.size(); ++i) {
        DT_LOG_TRACE(get_logging_priority(), "Adding a new clustered gamma");
        snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
        hPT.grab().set_track_id(ptd_.get_number_of_particles());
        hPT.grab().set_charge(snemo::datamodel::particle_track::neutral);
        ptd_.add_particle(hPT);

        // const cluster_type & a_cluster = the_reconstructed_clusters.at(i);
        const cluster_type & a_cluster = the_reconstructed_gammas.at(i);
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
      DT_LOG_DEBUG(get_logging_priority(), "Number of gammas : " << the_reconstructed_gammas.size());
      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return;
    }

    void gamma_clustering_module::_get_geometrical_neighbours(const snemo::datamodel::calibrated_calorimeter_hit & hit_,
                                                              const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & hits_,
                                                              cluster_type & cluster_,
                                                              gid_list_type & registered_calos_) const
    {
      const geomtools::geom_id & a_gid = hit_.get_geom_id();

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

      uint8_t mask = snemo::geometry::utils::NEIGHBOUR_NONE;
     if (_cluster_grid_mask_ == "first") {
        mask = snemo::geometry::utils::NEIGHBOUR_FIRST;
      } else if (_cluster_grid_mask_ == "second") {
        mask = snemo::geometry::utils::NEIGHBOUR_SECOND;
      } else if (_cluster_grid_mask_ == "diagonal") {
        mask = snemo::geometry::utils::NEIGHBOUR_DIAG;
      } else if (_cluster_grid_mask_ == "side") {
        mask = snemo::geometry::utils::NEIGHBOUR_SIDE;
      } else if (_cluster_grid_mask_ == "none") {
        mask = snemo::geometry::utils::NEIGHBOUR_NONE;
      } else {
        DT_THROW_IF(true, std::logic_error, "Unkown neighbour mask '" << _cluster_grid_mask_ << "' !")
          }
      if (calo_locator.is_calo_block_in_current_module(a_gid)) {
        calo_locator.get_neighbours_ids(a_gid, the_neighbours, mask);
      } else if (xcalo_locator.is_calo_block_in_current_module(a_gid)) {
        xcalo_locator.get_neighbours_ids(a_gid, the_neighbours, mask);
      } else if (gveto_locator.is_calo_block_in_current_module(a_gid)) {
        gveto_locator.get_neighbours_ids(a_gid, the_neighbours, mask);
      } else {
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << a_gid << "' does not match any scintillator block !");
      }

      for (gid_list_type::const_iterator ineighbour = the_neighbours.begin();
           ineighbour != the_neighbours.end(); ++ineighbour) {
        const geomtools::geom_id & a_gid = *ineighbour;
        if (std::find(registered_calos_.begin(), registered_calos_.end(), a_gid)
            != registered_calos_.end()) {
          continue;
        }

        // 2015/02/06 SC : using lambda function in C++14 (really easier)
        // if (std::find_if(cch.begin(), cch.end(), [ineighbour] (auto icalo)
        //                  {return ineighbour == icalo.get().get_geom_id();}) != cch.end())
        // Find if the eight neighbours belong to calibrated calo. hits
        geomtools::base_hit::has_geom_id_predicate hit_pred(a_gid);
        datatools::mother_to_daughter_predicate<geomtools::base_hit,
                                                snemo::datamodel::calibrated_calorimeter_hit> pred_M2D(hit_pred);
        datatools::handle_predicate<snemo::datamodel::calibrated_calorimeter_hit> pred_via_handle(pred_M2D);
        snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
          found = std::find_if(hits_.begin(), hits_.end(), pred_via_handle);
        if (found == hits_.end()) {
          continue;
        }

        registered_calos_.push_back(a_gid);
        cluster_.insert(std::make_pair(found->get().get_time(), *found));
        _get_geometrical_neighbours(found->get(), hits_, cluster_, registered_calos_);
      }
      return;
    }

    void gamma_clustering_module::_get_time_neighbours(cluster_type & cluster_,
                                                       cluster_collection_type & clusters_) const
    {
      if (cluster_.size() < 2) return;

      cluster_type::iterator it = cluster_.begin();
      for (;it != cluster_.end(); ++it) {
        const double current_time = it->first;
        const double next_time = boost::next(it)->first;
        const double delta_time = next_time - current_time;
        if (delta_time > _cluster_time_range_) {
          DT_LOG_TRACE(get_logging_priority(),
                       "Delta time > " << _cluster_time_range_/CLHEP::ns << " ns !!");
          break;
        }
      }
      if (it == cluster_.end()) return;

      // Create a new cluster to be checked later. To be done in this way
      // i.e. create first a new empty cluster and then add it to the collection
      // of clusters otherwise clusters_ stays in a frozen state if the cluster
      // is added first and then filled.
      cluster_type a_cluster;
      a_cluster.insert(boost::next(it), cluster_.end());
      cluster_.erase(boost::next(it), cluster_.end());
      clusters_.push_back(a_cluster);

      _get_time_neighbours(a_cluster, clusters_);
      return;
    }

    double gamma_clustering_module::_get_probability(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                                                     const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit)
    {
      geomtools::vector_3d head_position;
      geomtools::vector_3d tail_position;
      const geomtools::geom_id & head_gid = head_end_calo_hit.get_geom_id();
      const geomtools::geom_id & tail_gid = tail_begin_calo_hit.get_geom_id();

      const snemo::geometry::calo_locator & calo_locator   = _locator_plugin_->get_calo_locator();
      const snemo::geometry::xcalo_locator & xcalo_locator = _locator_plugin_->get_xcalo_locator();
      const snemo::geometry::gveto_locator & gveto_locator = _locator_plugin_->get_gveto_locator();


      if (calo_locator.is_calo_block_in_current_module(head_gid))
        calo_locator.get_block_position(head_gid, head_position);
      else if (xcalo_locator.is_calo_block_in_current_module(head_gid))
        xcalo_locator.get_block_position(head_gid, head_position);
      else if (gveto_locator.is_calo_block_in_current_module(head_gid))
        gveto_locator.get_block_position(head_gid, head_position);
      else
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << head_gid << "' does not match any scintillator block !");

      if (calo_locator.is_calo_block_in_current_module(tail_gid))
        calo_locator.get_block_position(tail_gid, tail_position);
      else if (xcalo_locator.is_calo_block_in_current_module(tail_gid))
        xcalo_locator.get_block_position(tail_gid, tail_position);
      else if (gveto_locator.is_calo_block_in_current_module(tail_gid))
        gveto_locator.get_block_position(tail_gid, tail_position);
      else
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << tail_gid << "' does not match any scintillator block !");

      const double track_length = (head_position-tail_position).mag();
      const double t1 = head_end_calo_hit.get_time();
      const double t2 = tail_begin_calo_hit.get_time();
      const double t_th = track_length / CLHEP::c_light;
      const double sigma_l = 120; //mm
      const double sigma_exp = pow(head_end_calo_hit.get_sigma_time(),2) + pow(tail_begin_calo_hit.get_sigma_time(),2) + pow(sigma_l/100,2);
      const double chi2 = pow(std::abs(t1 - t2) - t_th,2)/sigma_exp;

      return gsl_cdf_chisq_Q(chi2, 1);
    }

    bool gamma_clustering_module::_are_on_same_wall(const snemo::datamodel::calibrated_calorimeter_hit & head_end_calo_hit,
                                                    const snemo::datamodel::calibrated_calorimeter_hit & tail_begin_calo_hit)
    {
      const geomtools::geom_id & head_gid = head_end_calo_hit.get_geom_id();
      const geomtools::geom_id & tail_gid = tail_begin_calo_hit.get_geom_id();

      std::string head_label;
      std::string tail_label;

      int head_side = -1;
      int tail_side = -1;
      int head_wall = -1;
      int tail_wall = -1;

      gid_list_type the_head_second_neighbours;
      uint8_t second_mask = snemo::geometry::utils::NEIGHBOUR_SECOND; // NOT WORKING

      const snemo::geometry::calo_locator & calo_locator   = _locator_plugin_->get_calo_locator();
      const snemo::geometry::xcalo_locator & xcalo_locator = _locator_plugin_->get_xcalo_locator();
      const snemo::geometry::gveto_locator & gveto_locator = _locator_plugin_->get_gveto_locator();

      if (calo_locator.is_calo_block_in_current_module(head_gid)) {
        calo_locator.get_neighbours_ids(head_gid, the_head_second_neighbours, second_mask);
        head_side = calo_locator.extract_side(head_gid);
        head_label = snemo::datamodel::particle_track::vertex_on_main_calorimeter_label();
      }
      else if (xcalo_locator.is_calo_block_in_current_module(head_gid)) {
        xcalo_locator.get_neighbours_ids(head_gid, the_head_second_neighbours, second_mask);
        head_wall = xcalo_locator.extract_wall(head_gid);
        head_side = xcalo_locator.extract_side(head_gid);
        head_label = snemo::datamodel::particle_track::vertex_on_x_calorimeter_label();
      }
      else if (gveto_locator.is_calo_block_in_current_module(head_gid)) {
        gveto_locator.get_neighbours_ids(head_gid, the_head_second_neighbours, second_mask);
        head_wall = gveto_locator.extract_wall(head_gid);
        head_side = gveto_locator.extract_side(head_gid);
        head_label = snemo::datamodel::particle_track::vertex_on_gamma_veto_label();
      }
      else
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << head_gid << "' does not match any scintillator block !");

      if (calo_locator.is_calo_block_in_current_module(tail_gid)) {
        tail_side = calo_locator.extract_side(tail_gid);
        tail_label = snemo::datamodel::particle_track::vertex_on_main_calorimeter_label();
      }
      else if (xcalo_locator.is_calo_block_in_current_module(tail_gid)) {
        tail_wall = xcalo_locator.extract_wall(tail_gid);
        tail_side = xcalo_locator.extract_side(tail_gid);
        tail_label = snemo::datamodel::particle_track::vertex_on_x_calorimeter_label();
      }
      else if (gveto_locator.is_calo_block_in_current_module(tail_gid)) {
        tail_wall = gveto_locator.extract_wall(tail_gid);
        tail_side = gveto_locator.extract_side(tail_gid);
        tail_label = snemo::datamodel::particle_track::vertex_on_gamma_veto_label();
      }
      else
        DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << tail_gid << "' does not match any scintillator block !");

      /*** Tolerance for second neighbours  ***/
      std::cout << "test entry "<< the_head_second_neighbours.size() << std::endl;

      for (gid_list_type::const_iterator ineighbour = the_head_second_neighbours.begin();
           ineighbour != the_head_second_neighbours.end(); ++ineighbour)
        if(*ineighbour == tail_gid) {
          return false;
          DT_THROW_IF(true, std::logic_error,
                    "Current geom id '" << tail_gid << "' does not match any scintillator block !");
          std::cout << "Tolerance for second neighbours " << std::endl;
        }
        else
          std::cout << "test" << std::endl;

      /*** End second neighbours  ***/

      if(head_label == tail_label)
        {
          if(head_label == snemo::datamodel::particle_track::vertex_on_main_calorimeter_label())
            return head_side == tail_side;
          else
            return head_wall == tail_wall;
        }

      return false;
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
