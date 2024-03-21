#include <iostream>
#include <unistd.h>
#include <string>
#include "commands.hpp"


void printHelp(){
   std::cout << ("Usage:\n otter [command]") << std::endl;
   std::cout << ("      assemble      Run local assembly across a given set of target regions.") << std::endl;
   std::cout << ("      genotype      Generate genotypes for each local assembly.") << std::endl;
   std::cout << ("      wgat          Type whole-genome assemblies/alignments.") << std::endl;

   std::cout << ("DEPRECATED") << std::endl;  
   std::cout << ("      length        Output local assembly lengths per sample (superseded in 'genotype' command).") << std::endl;
   std::cout << ("      fa2sam        Convert (otter) fasta file to SAM-format.\n") << std::endl;
   std::cout << ("version: 0.6.0\n") << std::endl;
 }

int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "assemble") command_assemble(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "wgat") command_wgat(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "fa2sam") command_fa2sam(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "genotype") command_genotype(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "length") command_length(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "cov") command_cov(argc - 1, &argv[1]);
    else printHelp();
  }

  return 0;
}
