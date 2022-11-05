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

using namespace std;

#include <math/wide_decimal/decwide_t.h>

using dec101_t = math::wide_decimal::decwide_t<INT32_C(101), std::uint32_t, void>;
#define double dec101_t

#include "tap_approaches.h"
#include "class_matrix.h"


class tap_tr_equilibrium : public tap_od_equilibrium
{
protected:
	const double alpha = 1e-6, beta = 1e-6;
	vector <vector <double>> od_route_flow;
public:
	tap_tr_equilibrium(string tr) : tap_od_equilibrium(tr) 
	{	
		od_route_flow.resize(number_of_od);
	}
	bool find_new_route(int t)
	{
		bool route_found = tap_od_equilibrium::find_new_route(t);
		if (route_found)
			if (od_routes[t].size() > 1)
				od_route_flow[t].push_back(0);
			else
				od_route_flow[t].push_back(get<2>(od_pairs[t]));
		return route_found;
	}
	bool check_zero_flow(double& flow_1, vector <int>& route_1, vector <int>& route_2)
	{
		set <int> edges_1_temp, edges_2_temp;
		vector <int> edges_1, edges_2;
		for (auto edge : route_1)
			edges_1_temp.insert(edge);
		for (auto edge : route_2)
			edges_2_temp.insert(edge);
		for (auto edge : route_1)
			if (edges_2_temp.find(edge) == edges_2_temp.end())
				edges_1.push_back(edge);
		for (auto edge : route_2)
			if (edges_1_temp.find(edge) == edges_1_temp.end())
				edges_2.push_back(edge);
		if (get_edges_delay(edges_1, -flow_1) >= get_edges_delay(edges_2, flow_1))
			return true;
		else
			return false;
	}
	void update_cur_flow(double flow, vector <int>& route)
	{
		for (auto& edge : route)
			current_flow[edge] += flow;
	}
	void balance_two_routes(double& flow_1, double& flow_2, vector <int>& route_1, vector <int>& route_2)
	{
		double l = 0, r = flow_1, m = l;
		set <int> edges_1_temp, edges_2_temp;
		vector <int> edges_1, edges_2;
		for (auto& edge : route_1)
			edges_1_temp.insert(edge);
		for (auto& edge : route_2)
			edges_2_temp.insert(edge);
		for (auto& edge : route_1)
			if (edges_2_temp.find(edge) == edges_2_temp.end())
				edges_1.push_back(edge);
		for (auto& edge : route_2)
			if (edges_1_temp.find(edge) == edges_1_temp.end())
				edges_2.push_back(edge);
		m = l + (r - l) / 2;
		while (abs(get_edges_delay(edges_1, -m) - get_edges_delay(edges_2, m)) > beta)
		{
			m = (1e5 * (r - l)) / 2;
			m += l * 1e5;
			m /= 1e5;
			if (get_edges_delay(edges_1, -m) - get_edges_delay(edges_2, m) >= 0)
				l = m;
			else
				r = m;
			//cout << m << ' ' << (r - l) << ' ' << (abs(get_edges_delay(edges_1, -m) - get_edges_delay(edges_2, m))) << '\n';
		}
		update_cur_flow(-flow_1, route_1);
		update_cur_flow(-flow_2, route_2);
		flow_1 -= m;
		flow_2 += m;
		update_cur_flow(flow_1, route_1);
		update_cur_flow(flow_2, route_2);

	}
	void balance_od_pair(int t)
	{
		int r_max, r_min;
		double temp, r_max_del, r_min_del;
		while (get_od_delta(t) > alpha)
		{
			r_max = -1; r_max_del = 0;
			r_min = -1; r_min_del = -1;
			for (int i = 0; i < od_routes[t].size(); i++)
			{
				temp = get_edges_delay(od_routes[t][i]);
				if (temp > r_max_del)
				{
					r_max = i;
					r_max_del = temp;
				}
				if (temp < r_min_del || r_min_del == -1)
				{
					r_min = i;
					r_min_del = temp;
				}
			}
			if (!check_zero_flow(od_route_flow[t][r_max], od_routes[t][r_max], od_routes[t][r_min]))
				balance_two_routes(od_route_flow[t][r_max], od_route_flow[t][r_min], od_routes[t][r_max], od_routes[t][r_min]);
			else
			{
				update_cur_flow(-od_route_flow[t][r_max], od_routes[t][r_max]);
				update_cur_flow(od_route_flow[t][r_max], od_routes[t][r_min]);
				od_route_flow[t][r_min] += od_route_flow[t][r_max];
				od_routes[t].erase(od_routes[t].begin() + r_max);
				od_route_flow[t].erase(od_route_flow[t].begin() + r_max);
			}
		}
	}
	void balance()
	{
		double delta = alpha + 1;
		while (delta > alpha)
		{
			for (int t = 0; t < number_of_od; t++)
				balance_od_pair(t);
			delta = get_delta();
			update_od_delays();
		}
	}
	void solve_flow() 
	{
		bool route_found = true;
		while (route_found)
		{
			route_found = false;
			for (int t = 0; t < number_of_od; t++)
				route_found = find_new_route(t) || route_found;
			balance();
			//cout << objective_function() << '\n';
		}
	}
};

class tap_matrix_equilibrium : public tap_od_equilibrium
{
protected:
	matrix <double> At;
	vector <vector <double>> od_flow;
	const double eps = 1e-10, alpha = 1e-10;
public:
	void generate_At()
	{
		At = matrix <double>(number_of_nodes, number_of_edges);
		for (int i = 0; i < number_of_edges; i++)
		{
			At[init_node[i]][i] = 1;
			At[term_node[i]][i] = -1;
		}
	}
	tap_matrix_equilibrium(string tr) : tap_od_equilibrium(tr)
	{
		generate_At();
		od_flow.resize(number_of_od, vector <double>(number_of_edges, 0));
	}
	double delay_der(int i, double extra_flow = 0)
	{
		double flow_on_edge = current_flow[i] + extra_flow;
		return free_flow_time[i] * 0.15 * 4 * pow(flow_on_edge / capacity[i], 3) / capacity[i];
	}
	bool find_new_route(int t)
	{
		bool route_found = tap_od_equilibrium::find_new_route(t);
		if (route_found && (od_routes[t].size() == 1))
			for (auto edge : od_routes[t][0])
				od_flow[t][edge] += get<2>(od_pairs[t]);
		return route_found;
	}
	void flow_change(int t, bool fl)
	{
		int temp = -1;
		if (!fl)
			temp = 1;
		for (int i = 0; i < number_of_edges; i++)
			current_flow[i] += od_flow[t][i] * temp;
	}
	tuple < map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <double> > compression(int t)
	{
		int origin = get<0>(od_pairs[t]), count = 0;
		set <int> au_edges, au_nodes;
		map <int, int> M_nodes, M_nodes_inv, M_edges, M_edges_inv;
		matrix <double> A;
		vector <vector <int>>& routes = od_routes[t];
		for (int i = 0; i < routes.size(); i++)
		{
			for (auto edge : routes[i])
				if (au_edges.find(edge) == au_edges.end())
				{
					au_edges.insert(edge);
					if (au_nodes.find(init_node[edge]) == au_nodes.end())
					{
						au_nodes.insert(init_node[edge]);
						M_nodes[init_node[edge]] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = init_node[edge];
					}
					if (au_nodes.find(term_node[edge]) == au_nodes.end())
					{
						au_nodes.insert(term_node[edge]);
						M_nodes[term_node[edge]] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = term_node[edge];
					}
				}
		}
		int n = au_nodes.size() - 1, m = au_edges.size();
		for (auto now : au_nodes)
			if (M_nodes[now] == n)
				swap(M_nodes[origin], M_nodes[now]);
		A = matrix <double>(n, m);
		count = 0;
		for (auto& edge : au_edges)
		{
			M_edges[edge] = count;
			M_edges_inv[count] = edge;
			if (init_node[edge] != origin)
				A[M_nodes[init_node[edge]]][count] = 1;
			if (term_node[edge] != origin)
				A[M_nodes[term_node[edge]]][count] = -1;
			count++;
		}
		return { M_nodes, M_nodes_inv, M_edges, M_edges_inv, A };
	}
	matrix <double> create_G(matrix <double>& flow, map <int, int>& M_edges_inv)
	{
		int n = flow.size().first;
		vector <vector <double>> Temp(n, vector <double>(n, 0));
		for (int i = 0; i < n; i++)
			Temp[i][i] = 1 / delay_der(M_edges_inv[i], flow[i][0]);
		return matrix <double>(n, n, Temp);
	}
	matrix <double> create_g(matrix <double>& flow, map <int, int>& M_edges_inv)
	{
		matrix <double> g = flow;
		for (int i = 0; i < flow.size().first; i++)
			g[i][0] = delay(M_edges_inv[i], flow[i][0]);
		return g;
	}
	void update_od_del(int t)
	{
		double delay = 0;
		if (od_routes[t].size() > 0)
			delay = get_edges_delay(od_routes[t][0]);
		for (auto& route : od_routes[t])
			delay = min(delay, get_edges_delay(route));
		od_delays[t] = delay;
	}
	void balance_od_pair(int t)
	{
		flow_change(t, true);
		for (auto& now : od_flow[t])
			now = 0;
		vector <vector <int>> pos_routes = od_routes[t];
		od_routes[t] = {};
		while (pos_routes.size() > 0)
		{
			int best_route = 0;
			double del_best = get_edges_delay(pos_routes[0]), del_now;
			for (int i = 1; i < pos_routes.size(); i++)
			{
				del_now = get_edges_delay(pos_routes[i]);
				if (del_now < del_best)
				{
					best_route = i;
					del_best = del_now;
				}
			}
			if (od_routes.size() > 0)
			{
				if (od_delays[t] > del_best)
				{
					od_routes[t].push_back(pos_routes[best_route]);
					balance_odp_routes(t);
				}
			}
			else
			{
				od_routes[t].push_back(pos_routes[best_route]);
				for (auto edge : od_routes[t][0])
					od_flow[t][edge] += get<2>(od_pairs[t]);
			}
			swap(pos_routes[best_route], pos_routes[pos_routes.size() - 1]);
			pos_routes.pop_back();
			update_od_del(t);
		}
	}
	// changed criteria
	void balance_odp_routes(int t)
	{
		flow_change(t, true);
		tuple <map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <double> > comp = compression(t);
		map <int, int> M_nodes = get<0>(comp), M_nodes_inv = get<1>(comp), M_edges = get<2>(comp), M_edges_inv = get<3>(comp);
		matrix <double> A = get<4>(comp), G, A_trans, E, g, Temp;
		int n = A.size().first, m = A.size().second;
		matrix <double>  flow_now(m, 1), flow_prev(m, 1);
		for (auto& route : od_routes[t])
			for (auto& edge : route)
				flow_now[M_edges[edge]][0] += get<2>(od_pairs[t]) / od_routes[t].size();
		double od_delta = eps + 1;
		while (od_delta > eps)
		{
			flow_prev = flow_now;
			G = create_G(flow_now, M_edges_inv);
			A_trans = A.trans();
			E = matrix<double>::create_e(m, m);
			g = create_g(flow_now, M_edges_inv);
			Temp = (A * G * A_trans).inv();
			Temp = A_trans * Temp * A * G;
			Temp = E - Temp;
			Temp = G * Temp * g;
			flow_now = flow_now - Temp;
			od_flow[t] = vector <double>(number_of_edges, 0);
			for (int i = 0; i < m; i++)
				od_flow[t][M_edges_inv[i]] = flow_now[i][0];
			flow_change(t, false);
			od_delta = get_od_delta(t);
			flow_change(t, true);
			//cout << od_delta << ' ' << (flow_now - flow_prev).norm2() << '\n';
		}
		od_flow[t] = vector <double>(number_of_edges, 0);
		for (int i = 0; i < m; i++)
			od_flow[t][M_edges_inv[i]] = flow_now[i][0];
		flow_change(t, false);
	}
	double get_delta()
	{
		double delta = 0;
		for (int t = 0; t < number_of_od; t++)
			delta = max(delta, get_od_delta(t));
		return delta;
	}
	void update_od_delays()
	{
		for (int t = 0; t < number_of_od; t++)
			update_od_del(t);
	}
	void balance()
	{
		double delta = alpha + 1;
		while (delta > alpha)
		{
			for (int t = 0; t < number_of_od; t++)
				if (od_routes[t].size() > 1)
					balance_od_pair(t);
			delta = get_delta();
			update_od_delays();
		}
	}
	void solve_flow()
	{
		bool route_found = true;
		while (route_found)
		{
			route_found = false;
			for (int t = 0; t < number_of_od; t++)
				route_found |= find_new_route(t);
			balance();
		}
	}
};
template <typename method>
void test_tap(string test)
{
	method tap_test(test);
	tap_test.solve_flow();
	tap_test.show_stat();
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
		test_tap<tap_matrix_equilibrium>(line);
		auto end = chrono::high_resolution_clock::now();
		chrono::duration<float> duration = end - start;
		cout << "Duration: " << duration.count() << " sec \n";

	}
}

int main()
{
	std::cout << std::setprecision(10);
	test_tap_solutions();
	return 0;
}
