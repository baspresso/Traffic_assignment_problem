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

using namespace std;

#include <math/wide_decimal/decwide_t.h>
#include <Eigen/Dense>
#include "tap_approaches.h"

using dec101_t = math::wide_decimal::decwide_t<INT32_C(101), std::uint32_t, void>;

template <typename Method>
void TrafficAssignmentTest(string test) {
  Method tap_test(test);
  tap_test.SolveFlow();
  //tap_test.ShowStatistics();
}
void TrafficAssignmentTest_solutions() {
  string line;
  ifstream in("tests.txt");
  while (getline(in, line)) {
    auto start = chrono::high_resolution_clock::now();
    cout << "test: " << line << '\n';
    //TrafficAssignmentTest<tap_tr_equilibrium>(line);
    TrafficAssignmentTest<traffic_assignment::RouteBasedApproach <dec101_t>>(line);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<float> duration = end - start;
    cout << "Duration: " << duration.count() << " sec \n";
  }
}

int main() {
  std::cout << std::setprecision(15);
  TrafficAssignmentTest_solutions();
  return 0;
}
