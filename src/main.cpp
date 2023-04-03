#include <iostream>
#include <unistd.h>
#include <string>
#include "commands.hpp"


/**
 *  Function to pring out general (help) message of gvt
 */
void printHelp(){
   std::cout << ("Usage:\n otter [command]") << std::endl;
   std::cout << ("      assemble      Run local assembly across a given set of target regions.") << std::endl;
 }

/**
 * Main method. Parse input commands
 */
int main(int argc, char **argv){
  if(argc == 1) printHelp();
  else{
    if(std::string(argv[1]) == "assemble") command_assemble(argc - 1, &argv[1]);
    else if(std::string(argv[1]) == "fa2sam") command_fa2sam(argc - 1, &argv[1]);
    else printHelp();
  }

  return 0;
}
