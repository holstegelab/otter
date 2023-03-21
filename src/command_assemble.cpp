#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "assemble.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_assemble_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
    options
      .allow_unrecognised_options()
      .show_positional_help()
      .set_width(120)
      .set_tab_expansion()
      .custom_help("[parameters] <BAM>")
      .add_options()
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>())
      ("r, reference", "Path to reference genome.", cxxopts::value<std::string>()->default_value(""))
      ("reads-only", "Output only (partial)spanning reads.", cxxopts::value<bool>()->default_value("false"))
      ("p, non-primary", "Use non-primmary read-alignments.", cxxopts::value<bool>()->default_value("false"))
      ("o, offset", "Extend coords by this amount", cxxopts::value<int>()->default_value("50"))
      ("a, max-alleles", "Maximum alleles allowed.", cxxopts::value<int>()->default_value("2"))
      ("m, mapq", "Minimum mapping quality.", cxxopts::value<int>()->default_value("0"))
      ("c, max-cov", "Ignore regions with coverage above this value.", cxxopts::value<int>()->default_value("200"))
      ("e, max-error", "Maximum tolerable error.", cxxopts::value<double>()->default_value("0.025"))
      ("h, bandwidth", "KDE bandwidth.", cxxopts::value<double>()->default_value("0.01"))
      ("f, flank-size", "Length of flanking seq re-alignment.", cxxopts::value<int>()->default_value("100"))
      ("s, min-sim", "Minimum similarity during re-alignment.", cxxopts::value<double>()->default_value("0.9"))
      ("t, threads", "Total number of threads.", cxxopts::value<int>()->default_value("4"));
    //parse CLI arguments
    auto result = options.parse(argc, argv);
    std::vector<std::string> inputs;
    for(auto & i : result.unmatched()) inputs.push_back(i);
    //no BAM files provided, output help message
    if(inputs.empty()) std::cout << options.help();
    else {
      OtterOpts params;
      const std::string bed = result["bed"].as<std::string>();
      const std::string reference = result["reference"].as<std::string>();
      const bool reads_only = result["reads-only"].as<bool>();
      params.nonprimary = result["non-primary"].as<bool>();
      params.init_offset(result["offset"].as<int>());
      params.init_max_alleles(result["max-alleles"].as<int>());
      params.init_mapq(result["mapq"].as<int>());
      params.init_max_cov(result["max-cov"].as<int>());
      params.init_threads(result["threads"].as<int>());
      params.init_max_error(result["max-error"].as<double>());
      params.init_bandwidth(result["bandwidth"].as<double>());
      params.init_flank(result["flank-size"].as<int>());
      params.init_min_sim(result["min-sim"].as<double>());
      assemble(inputs, bed, reference, reads_only, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_assemble(int argc, char **argv){
  //parse CLI arguments
  command_assemble_parser(argc, argv);
  return 0;
}
