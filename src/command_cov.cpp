#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "ocov.hpp"
#include <string>
#include <vector>

/**
 * Method to parse CLI arguments using typical main parameters.
 */

void command_cov_parser(int argc, char** argv){
  try{
    cxxopts::Options options(argv[0]);
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <BAM>")
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL")
      ("o, offset", "Extend coords by this amount", cxxopts::value<int>()->default_value("31"))
      ("m, mapq", "Minimum mapping quality.", cxxopts::value<int>()->default_value("0"))
      ("p, non-primary", "Use non-primmary read-alignments.", cxxopts::value<bool>()->default_value("false"))
      ("R, read-group", "Prepend stdout with this read-group/sample string.", cxxopts::value<std::string>()->default_value(""))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("1"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      const std::string bed = result["bed"].as<std::string>();
      params.init_mapq(result["mapq"].as<int>());
      params.init_offset(result["offset"].as<int>());
      params.init_threads(result["threads"].as<int>());
      params.read_group = result["read-group"].as<std::string>();
      ocov(params, bed, inputs.front());
    };
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}


int command_cov(int argc, char **argv){
  //parse CLI arguments
  command_cov_parser(argc, argv);
  return 0;
}
