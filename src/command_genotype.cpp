#include <vector>
#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "genotype.hpp"

/**
 * Method to parse CLI arguments using typical main parameters.
 */
void command_genotype_parser(int argc, char** argv){
  cxxopts::Options options(argv[0]);
  try{
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
      ("r, reference", "Use this reference genome to output in VCF format.", cxxopts::value<std::string>()->default_value(""))
      ("v, vcf", "Output in VCF format by using this sample name corresponding to reference genome already merged in input BAM file.", cxxopts::value<std::string>()->default_value(""))
      ("L, dist-length", "Trigger length-distance function for Hamiltonian cycle given comma-separated parameters: <length>,<fraction>,<se>.", cxxopts::value<std::string>()->default_value("1000,0.5,500"))
      ("m, total-mincov", "Minimum total coverage per region", cxxopts::value<int>()->default_value("0"))
      ("c, allele-mincov", "Minimum coverage per allele sequence", cxxopts::value<int>()->default_value("0"))
      ("e, max-error", "Minimum similarity during re-alignment.", cxxopts::value<double>()->default_value("0.05"))
      ("l, length", "Output multi-column TSV representing join, min, and max sequence lengths as genotypes", cxxopts::value<bool>()->default_value("false"))
      ("scaled", "Compute scaled genotypes", cxxopts::value<bool>()->default_value("false"))
      ("s, summary", "Output only summary information.", cxxopts::value<bool>()->default_value("false"))
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
      params.init_max_error(result["max-error"].as<double>());
      params.init_threads(result["threads"].as<int>());
      const int ac_mincov = result["allele-mincov"].as<int>();
      const int tc_mincov = result["total-mincov"].as<int>();
      const bool is_summary = result["summary"].as<bool>();
      const bool is_length = result["length"].as<bool>();
      const bool is_scaled = result["scaled"].as<bool>();
      const std::string reference = result["reference"].as<std::string>();
      const std::string vcf = result["vcf"].as<std::string>();
      params.init_dist_length(result["dist-length"].as<std::string>());
      genotype(inputs.front(), bed, params, ac_mincov, tc_mincov, is_summary, is_length, is_scaled, reference, vcf);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::OptionException& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_genotype(int argc, char **argv){
  //parse CLI arguments
  command_genotype_parser(argc, argv);
  return 0;
}
