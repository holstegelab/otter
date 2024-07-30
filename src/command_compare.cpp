#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "compare.hpp"
#include <vector>
#include <iostream>

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_compare_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <BAM> <BAM>")
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL")
      ("R, sample-name", "Sample name.", cxxopts::value<std::string>()->default_value(""))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("1"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty() || inputs.size() < 2) std::cout << options.help();
    else {
      OtterOpts params;
      std::string bed;
      if(result.count("bed")) bed = result["bed"].as<std::string>();
      else {
        std::cerr << "[ERROR] '--bed' parameter required\n" << options.help() << std::endl;
        exit(1);
      }
      params.read_group = result["sample-name"].as<std::string>();
      params.init_threads(result["threads"].as<int>());
      compare(params, bed, inputs[0], inputs[1]);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::exceptions::exception& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_compare(int argc, char **argv){
  //parse CLI arguments
  command_compare_parser(argc, argv);
  return 0;
}
