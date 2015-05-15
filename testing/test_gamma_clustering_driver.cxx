// test_gamma_clustering_driver.cxx

// Standard library:
#include <cstdlib>
#include <iostream>
#include <string>
#include <exception>

// - Bayeux/geomtools:
#include <bayeux/geomtools/manager.h>

// This project:
#include <falaise/falaise.h>
#include <falaise/resource.h>
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/reconstruction/gamma_clustering_driver.h>

int main()
{
  FALAISE_INIT();
  int error_code = EXIT_SUCCESS;
  try {
    std::clog << "Test program for the 'gamma_clustering_driver' class." << std::endl;

    namespace sdm = snemo::datamodel;
    namespace srt = snemo::reconstruction;

    // Particle track data bank :
    sdm::particle_track_data PTD;
    sdm::calibrated_calorimeter_hit::collection_type & the_calos
      = PTD.grab_non_associated_calorimeters();

    // Fake cluster
    {
      datatools::handle<sdm::calibrated_calorimeter_hit> hCCH;
      hCCH.reset(new sdm::calibrated_calorimeter_hit);
      geomtools::geom_id gid(1302, 0, 1, 1, 4, 0);
      gid.set_any(4);
      hCCH.grab().set_geom_id(gid);
      hCCH.grab().set_time(0.0*CLHEP::ns);
      the_calos.push_back(hCCH);
    }
    {
      datatools::handle<sdm::calibrated_calorimeter_hit> hCCH;
      hCCH.reset(new sdm::calibrated_calorimeter_hit);
      geomtools::geom_id gid(1302, 0, 1, 1, 5, 0);
      gid.set_any(4);
      hCCH.grab().set_geom_id(gid);
      hCCH.grab().set_time(1.0*CLHEP::ns);
      hCCH.grab().set_sigma_time(0.5*CLHEP::ns);
      the_calos.push_back(hCCH);
    }
    // Add another calorimeter hit internally in time
    {
      datatools::handle<sdm::calibrated_calorimeter_hit> hCCH;
      hCCH.reset(new sdm::calibrated_calorimeter_hit);
      geomtools::geom_id gid(1302, 0, 0, 1, 5, 0);
      gid.set_any(4);
      hCCH.grab().set_geom_id(gid);
      hCCH.grab().set_time(2.0*CLHEP::ns);
      hCCH.grab().set_sigma_time(0.5*CLHEP::ns);
      the_calos.push_back(hCCH);
    }
    // Finally, create a random time calorimeter hits
    {
      datatools::handle<sdm::calibrated_calorimeter_hit> hCCH;
      hCCH.reset(new sdm::calibrated_calorimeter_hit);
      geomtools::geom_id gid(1302, 0, 0, 4, 8, 0);
      gid.set_any(4);
      hCCH.grab().set_geom_id(gid);
      hCCH.grab().set_time(5.0*CLHEP::ns);
      hCCH.grab().set_sigma_time(0.5*CLHEP::ns);
      the_calos.push_back(hCCH);
    }

    // Dump particle track data
    PTD.tree_dump();

    // Need an effective geometry manager ot make use of locators
    const std::string gmanager_config_file = falaise::get_resource_dir(true) +
      "/config/snemo/demonstrator/geometry/3.0/manager.conf";
    // Load properties from the configuration file
    datatools::properties gmanager_config;
    datatools::properties::read_config(gmanager_config_file,
                                       gmanager_config);
    geomtools::manager geo_mgr;
    geo_mgr.initialize(gmanager_config);

    srt::gamma_clustering_driver GCD;
    GCD.set_geometry_manager(geo_mgr);
    datatools::properties GCD_config;
    GCD_config.store("BGB.logging.priority", "debug");
    GCD.initialize(GCD_config);
    GCD.process(PTD);

  } catch (std::exception & x) {
    std::cerr << "error: " << x.what() << std::endl;
    error_code = EXIT_FAILURE;
  } catch (...) {
    std::cerr << "error: " << "unexpected error !" << std::endl;
    error_code = EXIT_FAILURE;
  }
  FALAISE_FINI();
  return error_code;
}
