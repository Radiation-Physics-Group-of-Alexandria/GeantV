#include <iostream>

#include "RngBenchmarker.h"

using namespace vecphys;

int main(int argc, char *argv[])
{
  // default run
  int nsample = 1000000;
  int nrepetition =  20;

  if (argc >= 2) nsample = atoi(argv[1]);
  if (argc >= 3) nrepetition = atoi(argv[2]);

  RngBenchmarker tester;
  tester.SetNSample(nsample);
  tester.SetRepetition(nrepetition);

  int status = tester.RunBenchmark();

  if (status == 1)
    std::cout << "Run RngBenchmark Failed" << std::endl;
  return status;
}
