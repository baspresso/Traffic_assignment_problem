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

using namespace std;

#include <math/wide_decimal/decwide_t.h>

using dec101_t = math::wide_decimal::decwide_t<INT32_C(101), std::uint32_t, void>;

#include "class_matrix.h"
#include "tap_approaches.h"


template <typename method>
void test_tap(string test)
{
	method tap_test(test);
	tap_test.solve_flow();
	//tap_test.show_stat();
}
void test_tap_solutions()
{
	string line;
	ifstream in("tests.txt");
	while (getline(in, line))
	{
		auto start = chrono::high_resolution_clock::now();
		cout << "test: " << line << '\n';
		//test_tap<tap_tr_equilibrium>(line);
		test_tap<tap <dec101_t>>(line);
		auto end = chrono::high_resolution_clock::now();
		chrono::duration<float> duration = end - start;
		cout << "Duration: " << duration.count() << " sec \n";

	}
}

int main()
{
	std::cout << std::setprecision(15);
	test_tap_solutions();
	return 0;
}
