#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "wgat.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_wgat_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <BAM|FASTA[.GZ]>")
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL")
      ("sam", "Output in SAM-format. Requres '-R' parameter.", cxxopts::value<bool>()->default_value("false"))
      ("R, read-group", "Sample name (required when using '--sam').", cxxopts::value<std::string>()->default_value(""))
      ("o, offset", "Extend coords by this amount", cxxopts::value<int>()->default_value("31"))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("1"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      const std::string bed = result["bed"].as<std::string>();
      params.is_sam = result["sam"].as<bool>();
      if(params.is_sam){
        params.read_group = result["read-group"].as<std::string>();
        if(params.read_group.empty()) {
          std::cout << "ERROR: '--sam' requires '-R' parameter\n" << options.help() << std::endl;
          exit(1);
        }
      }
      params.init_offset(result["offset"].as<int>());
      params.init_threads(result["threads"].as<int>());
      wgat(inputs.front(), bed, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_wgat(int argc, char **argv){
  //parse CLI arguments
  command_wgat_parser(argc, argv);
  return 0;
}
