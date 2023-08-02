#include "commands.hpp"
#include "cxxopts.hpp"
#include "length.hpp"
#include <string>
#include <vector>

/**
 * Method to parse CLI arguments using typical main parameters.
 */

void command_length_parser(int argc, char** argv){
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
      ("m, total-mincov", "Minimum total coverage per region", cxxopts::value<int>()->default_value("0"))
      ("c, allele-mincov", "Minimum coverage per allele sequence", cxxopts::value<int>()->default_value("0"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    if(inputs.empty()) std::cout << options.help();
    else {
      const std::string bed = result["bed"].as<std::string>();
      const int ac_mincov = result["allele-mincov"].as<int>();
      const int tc_mincov = result["total-mincov"].as<int>();
      length(inputs.front(), bed, ac_mincov, tc_mincov);
    };
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}


int command_length(int argc, char **argv){
  //parse CLI arguments
  command_length_parser(argc, argv);
  return 0;
}
