#include "anbamfilehelper.hpp"
#include "sam.h"
#include "commands.hpp"
#include <iostream>
#include <unistd.h>
#include <string>


const std::string OTTER_VERSION = "v1.0";

void printHelp(){
   std::cout << ("Usage:\n otter [command]") << std::endl;
   std::cout << ("      assemble      Locally assembly a given set of target regions.") << std::endl;
   std::cout << ("      genotype      Genotype target regions across one or more samples.") << std::endl;
   std::cout << ("      wgat          Genotype target regions in a whole-genome aligned assembly.") << std::endl;
   std::cout << ("      version       Output current version.\n") << std::endl;
 }

int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "assemble") command_assemble(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "wgat") command_wgat(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "genotype") command_genotype(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "vcf2mat") command_vcf2mat(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "compare") command_compare(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "version") std::cout << OTTER_VERSION << '\n';
    else printHelp();
  }

  return 0;
}
