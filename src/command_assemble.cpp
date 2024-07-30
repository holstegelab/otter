#include "commands.hpp"
#include "cxxopts.hpp"
#include "otter_opts.hpp"
#include "assemble.hpp"
#include <iostream>
#include <vector>

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
      .add_options(" REQUIRED")
      ("b, bed", "BED-formatted file of target regions.", cxxopts::value<std::string>())
      ("R, sample-name", "Sample name.", cxxopts::value<std::string>());
    options
      .add_options("OPTIONAL (INPUT/OUTPUT)")
      ("r, reference", "Path to reference genome (used for local realignments).", cxxopts::value<std::string>()->default_value(""))
      ("fasta", "Output in FASTA-format.", cxxopts::value<bool>()->default_value("false"))
      ("haps", "Usa haplotag information.", cxxopts::value<bool>()->default_value("false"))
      ("reads-only", "Output only (partial)spanning reads.", cxxopts::value<bool>()->default_value("false"))
      ("p, non-primary", "Use non-primmary read-alignments.", cxxopts::value<bool>()->default_value("false"))
      ("l, omit-nonspanning", "Omit non-spanning reads.", cxxopts::value<bool>()->default_value("false"))
      ("debug", "Turn on debug mode.", cxxopts::value<bool>()->default_value("false"));
    options
      .add_options("OPTIONAL (HEURISTICS)")
      ("o, offset", "Extend start/end by this amount 'INT', or extend separately by these amounts 'INT,INT'", cxxopts::value<std::string>()->default_value("1,0"))
      ("a, max-alleles", "Maximum alleles allowed.", cxxopts::value<int>()->default_value("2"))
      ("m, mapq", "Minimum mapping quality.", cxxopts::value<int>()->default_value("0"))
      ("q, read-quality", "Minimium (PacBio) read-quality.", cxxopts::value<double>()->default_value("0"))
      ("c, max-cov", "Ignore regions with coverage above this value.", cxxopts::value<int>()->default_value("200"))
      ("F, cov-fraction", "Minimum coverage fraction per sequence.", cxxopts::value<double>()->default_value("0.2"))
      ("A, cov-fraction-large", "Alternative minimum coverage fraction (INT,DOUBLE).", cxxopts::value<std::string>()->default_value("500,0.1"))
      ("e, max-error", "Maximum tolerable error.", cxxopts::value<double>()->default_value("0.01"))
      ("h, bandwidth", "KDE bandwidth.", cxxopts::value<double>()->default_value("0.01"))
      ("f, flank-size", "Length of flanking seq re-alignment.", cxxopts::value<int>()->default_value("100"))
      ("s, min-sim", "Minimum similarity during re-alignment.", cxxopts::value<double>()->default_value("0.9"))
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
      const std::string reference = result["reference"].as<std::string>();
      bool reads_only = result["reads-only"].as<bool>();
      params.nonprimary = result["non-primary"].as<bool>();
      params.omitnonspanning = result["omit-nonspanning"].as<bool>();
      params.is_fa = result["fasta"].as<bool>();
      params.ignore_haps = !result["haps"].as<bool>();
      params.init_offset(result["offset"].as<std::string>());
      params.init_max_alleles(result["max-alleles"].as<int>());
      params.init_mapq(result["mapq"].as<int>());
      params.init_read_quality(result["read-quality"].as<double>());
      params.init_max_cov(result["max-cov"].as<int>());
      params.init_min_cov_fraction(result["cov-fraction"].as<double>());
      params.init_threads(result["threads"].as<int>());
      params.init_max_error(result["max-error"].as<double>());
      params.init_bandwidth(result["bandwidth"].as<double>());
      params.init_flank(result["flank-size"].as<int>());
      params.init_min_sim(result["min-sim"].as<double>());
      params.init_min_cov_fraction2(result["cov-fraction-large"].as<std::string>());
      params.is_debug = result["debug"].as<bool>();
      assemble(inputs.front(), bed, reference, reads_only, params);
    }
  }

  //unable to make sense of CLI
  catch (const cxxopts::exceptions::exception& e) {
    std::cout << "Error parsing options: " << e.what() << '\n' << options.help() << std::endl;
    exit(1);
  }
}

int command_assemble(int argc, char **argv){
  //parse CLI arguments
  command_assemble_parser(argc, argv);
  return 0;
}
