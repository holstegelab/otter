#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "wgat.hpp"
#include <vector>
#include <iostream>

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
      .custom_help("[parameters] <BAM|FASTA[.GZ]> ...")
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>())
      ("R, sample-name", "Sample name.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL")
      ("fasta", "Output in FASTA-fromat.", cxxopts::value<bool>()->default_value("false"))
      ("o, offset", "Extend start/end by this amount 'INT', or extend separately by these amounts 'INT,INT'", cxxopts::value<std::string>()->default_value("1,0"))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("1"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      std::string bed;
      if(result.count("bed")) bed = result["bed"].as<std::string>();
      else {
        std::cerr << "[ERROR] '--bed' parameter required\n" << options.help() << std::endl;
        exit(1);
      }
      if(result.count("sample-name")) params.read_group = result["sample-name"].as<std::string>();
      else{
        std::cerr << "[ERROR] '--sample-name' parameter required\n" << options.help() << std::endl;
        exit(1);
      }
      params.is_fa = result["fasta"].as<bool>();
      params.init_offset(result["offset"].as<std::string>());
      params.init_threads(result["threads"].as<int>());
      wgat(params, inputs.front(), bed);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::exceptions::exception e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_wgat(int argc, char **argv){
  //parse CLI arguments
  command_wgat_parser(argc, argv);
  return 0;
}
