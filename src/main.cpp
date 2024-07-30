#include "anbamfilehelper.hpp"
#include "sam.h"
#include "commands.hpp"
#include <iostream>
#include <unistd.h>
#include <string>


void printHelp(){
   std::cout << ("Usage:\n otter [command]") << std::endl;
   
 }

int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "assemble") command_assemble(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "wgat") command_wgat(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "genotype") command_genotype(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "vcf2mat") command_vcf2mat(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "compare") command_compare(argc - 1, &argv[1]);
    else printHelp();
  }

  return 0;
}
