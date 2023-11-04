#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <stack>
#include <map>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <typeinfo>
#include <thread>
#include <mutex>
#include <cstdio>
#include <filesystem>

using namespace std;

#include <math/wide_decimal/decwide_t.h>
#include "Eigen\Dense"
#include "tap_approaches.h"

using dec101_t = math::wide_decimal::decwide_t<INT32_C(120), std::uint32_t, void>;

template <typename Method>
void TrafficAssignmentTest(string test) {
  auto start = chrono::high_resolution_clock::now();
  cout << "test: " << test << '\n';
  Method tap_test(test);
  tap_test.SolveFlow();
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end - start;
  cout << "Duration: " << duration.count() << " sec \n";
  //tap_test.ShowStatistics();
}
void TrafficAssignmentTestSolutions() {
  string line;
  ifstream in("tests.txt");
  while (getline(in, line)) {
    //TrafficAssignmentTest<traffic_assignment::SolutionCheck <dec101_t>>(line);
    //traffic_assignment::SolutionCheck <dec101_t> a(line);
    //TrafficAssignmentTest<traffic_assignment::RouteBasedApproach <dec101_t>>(line);
    TrafficAssignmentTest<traffic_assignment::NewRouteBasedApproach <dec101_t>>(line);
    //TrafficAssignmentTest<traffic_assignment::RouteBasedApproach <long double>>(line);
  }
}

int main() {
  //std::cout << "Current path is " << filesystem::current_path() << '\n'; // (1)
  TrafficAssignmentTestSolutions();
  return 0;
}
