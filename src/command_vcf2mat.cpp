#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "vcf2mat.hpp"
#include <vector>
#include <iostream>

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_vcf2mat_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <VCF[.GZ]>")
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL")
      ("k, kmer-size", "Kmer-size to use.", cxxopts::value<int>()->default_value("3"))
      ("t, threads", "Total threads to use.", cxxopts::value<int>()->default_value("1"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      const std::string bed = result["bed"].as<std::string>();
      int k_l = result["kmer-size"].as<int>();
      if(k_l < 1 || k_l > 32){
        std::cerr << "[ERROR] invalid '--kmer-size' (" << k_l << "). Needs to be 1 <= x <= 32." << options.help() << std::endl;
        exit(1);
      }
      params.init_threads(result["threads"].as<int>());
      vcf2mat(params, bed, inputs.front(), k_l);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::exceptions::exception& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_vcf2mat(int argc, char **argv){
  //parse CLI arguments
  command_vcf2mat_parser(argc, argv);
  return 0;
}
