#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "dist.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_dist_parser(int argc, char** argv){
  try{
    cxxopts::Options options(argv[0]);
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <FASTA>")
      .add_options()
      ("a, max-alleles", "Maximum alleles allowed.", cxxopts::value<int>()->default_value("2"))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("4"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      params.init_max_alleles(result["max-alleles"].as<int>());
      params.init_threads(result["threads"].as<int>());
      dist(inputs, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

int command_dist(int argc, char **argv){
  //parse CLI arguments
  command_dist_parser(argc, argv);
  return 0;
}
