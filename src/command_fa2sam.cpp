#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "fa2sam.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_fa2sam_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <ref.fa|gz> <assembly.fa|gz>")
      .add_options(" REQUIRED")
      ("r, reference", "Path to reference genome.", cxxopts::value<std::string>()->default_value(""));
    options
      .add_options("OPTIONAL")
      ("R, read-group", "Output with this read-group tag", cxxopts::value<std::string>()->default_value(""));
    
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    if(inputs.empty()) std::cout << options.help();
    else {
      const std::string reference = result["reference"].as<std::string>();
      const std::string read_group = result["read-group"].as<std::string>();
      fa2sam(reference, read_group, inputs);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_fa2sam(int argc, char **argv){
  //parse CLI arguments
  command_fa2sam_parser(argc, argv);
  return 0;
}
