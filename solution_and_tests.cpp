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
#include <typeinfo>

using namespace std;

#include <math/wide_decimal/decwide_t.h>

using dec101_t = math::wide_decimal::decwide_t<INT32_C(101), std::uint32_t, void>;

#include "tap_approaches.h"
#include "class_matrix.h"


template <typename v_type>
class tap_matrix_equilibrium : public tap <v_type>
{
protected:
	matrix <v_type> At;
	vector <vector <v_type>> od_flow;
	const v_type eps = 1e-7, alpha = 1e-7;
public:
	void generate_At()
	{
		At = matrix <v_type>(this->number_of_nodes, this->number_of_links);
		for (int i = 0; i < this->number_of_links; i++)
		{
			At[this->lnks[i].init][i] = 1;
			At[this->lnks[i].term][i] = -1;
		}
	}
	tap_matrix_equilibrium(string tr) : tap<v_type>(tr)
	{
		generate_At();
		od_flow.resize(this->number_of_od, vector <v_type>(this->number_of_links, 0));
	}
	bool find_new_route(int t)
	{
		bool route_found = tap<v_type>::find_new_route(t);
		if (route_found && (this->od_routes[t].size() == 1))
			for (auto now : this->od_routes[t][0])
				od_flow[t][now] += get<2>(this->od_pairs[t]);
		return route_found;
	}
	void flow_change(int t, bool fl)
	{
		int temp = -1;
		if (!fl)
			temp = 1;
		for (int i = 0; i < this->number_of_links; i++)
			this->lnks[i].flow += od_flow[t][i] * temp;
	}
	tuple < map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <v_type> > compression(int t)
	{
		int origin = get<0>(this->od_pairs[t]), cnt = 0;
		set <int> au_edges, au_nodes;
		map <int, int> M_nodes, M_nodes_inv, M_edges, M_edges_inv;
		matrix <v_type> A;
		vector <vector <int>>& routes = this->od_routes[t];
		for (int i = 0; i < routes.size(); i++)
		{
			for (auto now : routes[i])
				if (au_edges.find(now) == au_edges.end())
				{
					au_edges.insert(now);
					if (au_nodes.find(this->lnks[now].init) == au_nodes.end())
					{
						au_nodes.insert(this->lnks[now].init);
						M_nodes[this->lnks[now].init] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = this->lnks[now].init;
					}
					if (au_nodes.find(this->lnks[now].term) == au_nodes.end())
					{
						au_nodes.insert(this->lnks[now].term);
						M_nodes[this->lnks[now].term] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = this->lnks[now].term;
					}
				}
		}
		int n = au_nodes.size() - 1, m = au_edges.size();
		for (auto now : au_nodes)
			if (M_nodes[now] == n)
				swap(M_nodes[origin], M_nodes[now]);
		A = matrix <v_type>(n, m);
		cnt = 0;
		for (auto now : au_edges)
		{
			M_edges[now] = cnt;
			M_edges_inv[cnt] = now;
			if (this->lnks[now].init != origin)
				A[M_nodes[this->lnks[now].init]][cnt] = 1;
			if (this->lnks[now].term != origin)
				A[M_nodes[this->lnks[now].term]][cnt] = -1;
			cnt++;
		}
		return { M_nodes, M_nodes_inv, M_edges, M_edges_inv, A };
	}
	matrix <v_type> create_G(matrix <v_type>& flow, map <int, int>& M_edges_inv)
	{
		int n = flow.size().first;
		vector <vector <v_type>> Temp(n, vector <v_type>(n, 0));
		for (int i = 0; i < n; i++)
			Temp[i][i] = 1 / this->lnks[M_edges_inv[i]].delay_der(flow[i][0]);
		return matrix <v_type>(n, n, Temp);
	}
	matrix <v_type> create_g(matrix <v_type>& flow, map <int, int>& M_edges_inv)
	{
		matrix <v_type> g = flow;
		for (int i = 0; i < flow.size().first; i++)
			g[i][0] = this->lnks[M_edges_inv[i]].delay(flow[i][0]);
		return g;
	}
	void update_od_del(int t)
	{
		v_type delay = 0;
		if (this->od_routes[t].size() > 0)
			delay = this->get_links_delay(this->od_routes[t][0]);
		for (auto& route : this->od_routes[t])
			delay = min(delay, this->get_links_delay(route));
		this->od_delays[t] = delay;
	}
	void balance_od_pair(int t)
	{
		flow_change(t, true);
		for (auto& now : od_flow[t])
			now = 0;
		vector <vector <int>> pos_routes = this->od_routes[t];
		this->od_routes[t] = {};
		while (pos_routes.size() > 0)
		{
			int best_route = 0;
			v_type del_best = this->get_links_delay(pos_routes[0]), del_now;
			for (int i = 1; i < pos_routes.size(); i++)
			{
				del_now = this->get_links_delay(pos_routes[i]);
				if (del_now < del_best)
				{
					best_route = i;
					del_best = del_now;
				}
			}
			if (this->od_routes.size() > 0)
			{
				if (this->od_delays[t] > del_best)
				{
					this->od_routes[t].push_back(pos_routes[best_route]);
					balance_odp_routes(t);
				}
			}
			else
			{
				this->od_routes[t].push_back(pos_routes[best_route]);
				for (auto now : this->od_routes[t][0])
					od_flow[t][now] += get<2>(this->od_pairs[t]);
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
		tuple <map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <v_type> > comp = compression(t);
		map <int, int> M_nodes = get<0>(comp), M_nodes_inv = get<1>(comp), M_edges = get<2>(comp), M_edges_inv = get<3>(comp);
		matrix <v_type> A = get<4>(comp), G, A_trans, E, g, Temp;
		int n = A.size().first, m = A.size().second;
		matrix <v_type>  flow_now(m, 1), flow_prev(m, 1);
		for (auto& route : this->od_routes[t])
			for (auto& edge : route)
				flow_now[M_edges[edge]][0] += get<2>(this->od_pairs[t]) / this->od_routes[t].size();
		v_type od_delta = eps + 1;
		while (od_delta > eps)
		{
			flow_prev = flow_now;
			G = create_G(flow_now, M_edges_inv);
			A_trans = A.trans();
			E = matrix<v_type>::create_e(m, m);
			g = create_g(flow_now, M_edges_inv);
			Temp = (A * G * A_trans).inv();
			Temp = A_trans * Temp * A * G;
			Temp = E - Temp;
			Temp = G * Temp * g;
			flow_now = flow_now - Temp;
			od_flow[t] = vector <v_type>(this->number_of_links, 0);
			for (int i = 0; i < m; i++)
				od_flow[t][M_edges_inv[i]] = flow_now[i][0];
			flow_change(t, false);
			od_delta = this->get_od_delta(t);
			//this->show_stat();
			flow_change(t, true);
			//flow_now.show();
			//  << "!   " << get<2>(this->od_pairs[t]) << '\n';
			//cout << od_delta << ' ' << (flow_now - flow_prev).norm2() << '\n';
		}
		od_flow[t] = vector <v_type>(this->number_of_links, 0);
		for (int i = 0; i < m; i++)
			od_flow[t][M_edges_inv[i]] = flow_now[i][0];
		flow_change(t, false);
	}
	v_type get_delta()
	{
		v_type delta = 0;
		for (int t = 0; t < this->number_of_od; t++)
			delta = max(delta, this->get_od_delta(t));
		return delta;
	}
	void update_od_delays()
	{
		for (int t = 0; t < this->number_of_od; t++)
			update_od_del(t);
	}
	void balance()
	{
		v_type delta = alpha + 1;
		while (delta > alpha)
		{
			for (int t = 0; t < this->number_of_od; t++)
				if (this->od_routes[t].size() > 1)
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
			for (int t = 0; t < this->number_of_od; t++)
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
		test_tap<tap_matrix_equilibrium <dec101_t>>(line);
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
