template <class v_type>										// value type
struct link
{
public:
	int init, term, type;
	v_type capacity, length, free_flow_time, b, power, speed, toll, flow;
	link(int init, int term, int type, v_type capacity, v_type length,
		v_type free_flow_time, v_type b, v_type power, v_type speed, v_type toll) :
		init(init), term(term), type(type), capacity(capacity), length(length),
		free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll) { };
	v_type delay(v_type extra_flow = 0)
	{
		return free_flow_time + 0.15 * free_flow_time * pow(flow + extra_flow, 4) / pow(capacity, 4);
	}
	v_type delay_integ(v_type flow_on_link = -1.454555555)
	{
		if (flow_on_link == -1.454555555)
			flow_on_link = flow;
		return free_flow_time * (flow_on_link + 0.15 * capacity * pow(flow_on_link / capacity, 5) / 5);
	}
	v_type delay_der(v_type extra_flow = 0)
	{
		return free_flow_time * 0.15 * 4 * pow((flow + extra_flow) / capacity, 3) / capacity;
	}
	static v_type get_links_delay(vector <link <v_type>>& lnks, vector <int>& lnk_list, v_type extra_flow = 0)
	{
		v_type ans = 0;
		for (auto now : lnk_list)
			ans += lnks[now].delay(extra_flow);
		return ans;
	}
};

template <class v_type>
class od_pair
{
private:
	int origin, dest;
	v_type demand, cur_delay, next_delay, next_delta;
	vector <link <v_type>> &lnks;
	const vector <vector <int>> &adj_list;
	vector <vector <int>> cur_routes, next_routes;
	unordered_map <int, v_type> cur_lnks_flow, next_lnks_flow;
	const v_type eps = 1e-6;
	matrix <v_type> create_G(matrix <v_type>& temp_flow, map <int, int>& M_lnks_inv)
	{
		int n = temp_flow.size().first;
		vector <vector <v_type>> Temp(n, vector <v_type>(n, 0));
		for (int i = 0; i < n; i++)
			Temp[i][i] = 1 / lnks[M_lnks_inv[i]].delay_der(temp_flow[i][0] - cur_lnks_flow[M_lnks_inv[i]]);
		return matrix <v_type>(n, n, Temp);
	}
	matrix <v_type> create_g(matrix <v_type>& temp_flow, map <int, int>& M_lnks_inv)
	{
		matrix <v_type> g = temp_flow;
		for (int i = 0; i < temp_flow.size().first; i++)
			g[i][0] = lnks[M_lnks_inv[i]].delay(temp_flow[i][0] - cur_lnks_flow[M_lnks_inv[i]]);
		return g;
	}
	tuple < map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <v_type> > compression(vector <vector <int>>& routes)
	{
		int cnt = 0;
		set <int> au_edges, au_nodes;
		map <int, int> M_nodes, M_nodes_inv, M_edges, M_edges_inv;
		matrix <v_type> A;
		for (int i = 0; i < routes.size(); i++)
		{
			for (auto now : routes[i])
				if (au_edges.find(now) == au_edges.end())
				{
					au_edges.insert(now);
					if (au_nodes.find(lnks[now].init) == au_nodes.end())
					{
						au_nodes.insert(lnks[now].init);
						M_nodes[lnks[now].init] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = lnks[now].init;
					}
					if (au_nodes.find(lnks[now].term) == au_nodes.end())
					{
						au_nodes.insert(lnks[now].term);
						M_nodes[lnks[now].term] = au_nodes.size() - 1;
						M_nodes_inv[au_nodes.size() - 1] = lnks[now].term;
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
			if (lnks[now].init != origin)
				A[M_nodes[lnks[now].init]][cnt] = 1;
			if (lnks[now].term != origin)
				A[M_nodes[lnks[now].term]][cnt] = -1;
			cnt++;
		}
		return { M_nodes, M_nodes_inv, M_edges, M_edges_inv, A };
	}
	void update_next_data()
	{
		v_type del_min = 1e9, del_max = 0, del_route;
		if (next_routes.size() > 0)
		{
			del_route = 0;
			for (auto now : next_routes[0])
				del_route += lnks[now].delay(next_lnks_flow[now] - cur_lnks_flow[now]);
			del_min = del_route;
			del_max = del_route;
		}
		for (int i = 1; i < next_routes.size(); i++)
		{
			del_route = 0;
			for (auto now : next_routes[i])
				del_route += lnks[now].delay(next_lnks_flow[now] - cur_lnks_flow[now]);
			del_min = min(del_min, del_route);
			del_max = max(del_max, del_route);
		}
		next_delay = del_min;
		next_delta = del_max - del_min;
	}
	void balance_next_routes()
	{
		tuple <map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <v_type> > comp = compression(next_routes);
		map <int, int> M_nodes = get<0>(comp), M_nodes_inv = get<1>(comp), M_edges = get<2>(comp), M_edges_inv = get<3>(comp);
		matrix <v_type> A = get<4>(comp), G, A_trans, E, g, T;
		int n = A.size().first, m = A.size().second;

		matrix <v_type>  temp_flow_now(m, 1), temp_flow_prev(m, 1);
		for (int i = 0; i < next_routes.size(); i++)
			for (auto& now : next_routes[i])
				temp_flow_now[M_edges[now]][0] += demand / next_routes.size();
		next_delta = eps + 1;
		next_lnks_flow.clear();
		for (int i = 0; i < m; i++)
			next_lnks_flow[M_edges_inv[i]] = temp_flow_now[i][0];
		update_next_data();
		while (next_delta > eps)
		{
			temp_flow_prev = temp_flow_now;
			G = create_G(temp_flow_now, M_edges_inv);
			A_trans = A.trans();
			E = matrix<v_type>::create_e(m, m);
			g = create_g(temp_flow_now, M_edges_inv);
			T = (A * G * A_trans).inv();
			T = A_trans * T * A * G;
			T = E - T;
			T = G * T * g;
			temp_flow_now = temp_flow_now - T;
			next_lnks_flow.clear();
			for (int i = 0; i < m; i++)
				next_lnks_flow[M_edges_inv[i]] = temp_flow_now[i][0];
			update_next_data();
		}
	}
	void update_cur_delay()
	{
		v_type del_min = 1e9, del_route;
		if (cur_routes.size() > 0)
		{
			del_route = link<v_type>::get_links_delay(lnks, cur_routes[0]);
			del_min = del_route;
		}
		for (int i = 1; i < cur_routes.size(); i++)
		{
			del_route = link<v_type>::get_links_delay(lnks, cur_routes[i]);
			del_min = min(del_min, del_route);
		}
		cur_delay = del_min;
	}
public: 
	od_pair(int &origin, int &dest, v_type &demand, vector <link <v_type>> &lnks, vector <vector <int>> &adj_list) :
		origin(origin), dest(dest), demand(demand), lnks(lnks), adj_list(adj_list), cur_delay(0) { };
	bool find_new_route()
	{
		update_cur_delay();
		//cout << cur_delay << '\n';
		priority_queue <pair <v_type, int>, vector <pair <v_type, int>>, greater <pair <v_type, int>>> q;
		int u;
		q.push({ origin, -1 });
		unordered_set <int> processed;
		unordered_map <int, int> used_link;
		v_type temp;
		while (q.top().second == -1 || lnks[q.top().second].term != dest)
		{
			if (q.top().second != -1)
				u = lnks[q.top().second].term;
			else
				u = origin;
			if (processed.find(u) != processed.end())
			{
				q.pop();
				continue;
			}
			temp = q.top().first;
			if (processed.find(u) != processed.end())
				continue;
			processed.insert(u);
			used_link[u] = q.top().second;
			q.pop();
			for (auto now : adj_list[u])
				if (processed.find(lnks[now].term) == processed.end())
					q.push({ temp + lnks[now].delay(), now });
		}
		used_link[dest] = q.top().second;
		v_type delay_cur_route = q.top().first;
		if (delay_cur_route < cur_delay || cur_routes.size() == 0)
		{
			cur_delay = delay_cur_route;
			int now = dest;
			vector <int> new_route;
			while (now != origin)
			{
				new_route.push_back(used_link[now]);
				now = lnks[used_link[now]].init;
			}
			reverse(new_route.begin(), new_route.end());
			cur_routes.push_back(new_route);
			if (cur_routes.size() == 1)
			{
				for (auto now : cur_routes[0])
				{
					lnks[now].flow += demand;
					cur_lnks_flow[now] = demand;
				}
			}
			return true;
		}
		return false;
	}
	v_type get_cur_delta()
	{
		v_type del_min = 1e9, del_max = 0, del_route;
		if (cur_routes.size() > 0)
		{
			del_route = link<v_type>::get_links_delay(lnks, cur_routes[0]);
			del_min = del_route;
			del_max = del_route;
		}
		for (int i = 1; i < cur_routes.size(); i++)
		{
			del_route = link<v_type>::get_links_delay(lnks, cur_routes[i]);
			del_min = min(del_min, del_route);
			del_max = max(del_max, del_route);
		}
		return abs(del_max - del_min);
	}
	void calc_next_update()
	{
		vector <vector <int>> pos_routes = cur_routes;
		next_routes = {};
		next_delay = -1;
		next_lnks_flow = cur_lnks_flow;
		while (!pos_routes.empty())
		{
			int best_route = 0;
			v_type del_best = -1, del_now;
			for (int i = 0; i < pos_routes.size(); i++)
			{
				del_now = 0;
				for (auto now : pos_routes[i])
					del_now += lnks[now].delay(next_lnks_flow[now] - cur_lnks_flow[now]);
				if ((del_now < del_best) || (del_best == -1))
				{
					best_route = i;
					del_best = del_now;
				}
			}
			if ((del_best < next_delay) || (next_delay == -1))
			{
				next_routes.push_back(pos_routes[best_route]);
				balance_next_routes();
			}
			swap(pos_routes[best_route], pos_routes[pos_routes.size() - 1]);
			pos_routes.pop_back();
		}
	}
	void implement_next_update()
	{
		unordered_set <int> cur_lnks, next_lnks;
		for (int i = 0; i < cur_routes.size(); i++)
			for (auto now : cur_routes[i])
				cur_lnks.insert(now);
		for (auto now : cur_lnks)
			lnks[now].flow -= cur_lnks_flow[now];
		for (int i = 0; i < next_routes.size(); i++)
			for (auto now : next_routes[i])
				next_lnks.insert(now);
		for (auto now : next_lnks)
			lnks[now].flow += next_lnks_flow[now];
		cur_lnks_flow = next_lnks_flow;
		cur_routes = next_routes;
		update_next_data();
		cur_delay = next_delay;
	}
};

template <typename v_type>										// value type
class tap
{
private:
	vector <link <v_type>> lnks;
	vector <vector <int>> adj_list;							// links for each node
	vector <od_pair <v_type>> od_prs;
	int number_of_od, number_of_links, number_of_nodes;
	string test;
	v_type alpha = 1e-4;
public:
	bool isNumber(char c)
	{
		return c >= '0' && c <= '9';
	}
	v_type get_value(int& i, string& line)
	{
		while (!isNumber(line[i]))
			i++;
		string temp1 = "";
		while (isNumber(line[i]))
			temp1 += line[i++];
		if (line[i++] != '.')
			return stoi(temp1);
		string temp2 = "";
		while (isNumber(line[i]))
			temp2 += line[i++];
		return stold(temp1) + stold(temp2) * (pow(0.1, temp2.size()));
	}
	int get_value_int(int& i, string& line)
	{
		while (!isNumber(line[i]))
			i++;
		string temp1 = "";
		while (isNumber(line[i]))
			temp1 += line[i++];
		if (line[i++] != '.')
			return stoi(temp1);
		string temp2 = "";
		while (isNumber(line[i]))
			temp2 += line[i++];
		return stold(temp1);
	}
	void get_link(string& line)
	{
		int i = 0;
		number_of_links++;
		int init, term, type;
		v_type capacity, length, free_flow_time, b, power, speed, toll, flow;
		init = get_value_int(i, line) - 1; term = get_value_int(i, line) - 1;
		capacity = get_value(i, line);	length = get_value(i, line);
		free_flow_time = get_value(i, line); b = get_value(i, line);
		power = get_value(i, line);	speed = get_value(i, line);
		toll = get_value(i, line); type = get_value_int(i, line);
		lnks.push_back(link<v_type>(init, term, type, capacity, length, free_flow_time, b, power, speed, toll));
	}
	void get_od_pair(ifstream& in)
	{
		int i = 0, origin, dest;
		v_type demand;
		string line;
		getline(in, line);
		origin = get_value_int(i, line) - 1;
		getline(in, line);
		i = 0;
		for (int h = 0; h < number_of_nodes; h++)
		{
			dest = get_value_int(i, line) - 1;
			demand = get_value(i, line);
			i++;
			if (demand > 0)
			{
				number_of_od++;
				//cout << demand << '\n';
				od_prs.push_back(od_pair <v_type>(origin, dest, demand, lnks, adj_list));
			}
			if (i == line.size() - 1)
			{
				getline(in, line);
				i = 0;
			}
		}
	}
	void get_data()
	{
		string line;
		ifstream in(test + "_net.txt");
		while (getline(in, line))
			get_link(line);
		in = ifstream(test + "_trips.txt");
		getline(in, line);
		int i = 0;
		number_of_nodes = get_value_int(i, line);
		for (int j = 0; j < number_of_nodes; j++)
			get_od_pair(in);
		number_of_od = od_prs.size();
	}
	tap(string tn) : number_of_links(0), number_of_od(0), test(tn)
	{
		get_data();
		adj_list.resize(number_of_nodes);
		for (int i = 0; i < number_of_links; i++)
			adj_list[lnks[i].init].push_back(i);
	}
	v_type objective_function()
	{
		v_type ans = 0;
		for (int i = 0; i < number_of_links; i++)
			ans += lnks[i].delay_integ();
		return ans;
	}
	v_type get_best_answer()
	{
		v_type ans = 0, flow_on_link;
		string line;
		ifstream in(test + "_flow.txt");
		int j = 0;
		for (int i = 0; i < number_of_links; i++)
		{
			getline(in, line);
			j = 0;
			get_value(j, line); get_value(j, line);
			flow_on_link = get_value(j, line);
			ans += lnks[i].delay_integ(flow_on_link);
		}
		return ans;
	}
	void show_stat()
	{
		cout << "best ans: " << get_best_answer() << '\n';
		cout << "obj function: " << objective_function() << '\n';
		for (int i = 0; i < number_of_links; i++)
			cout << "from: " << lnks[i].init + 1 << " to " << lnks[i].term + 1 << " volume " << lnks[i].flow << " cost " <<
			lnks[i].delay() << '\n';
	}
	v_type get_delta()
	{
		v_type delta = 0;
		for (int t = 0; t < number_of_od; t++)
			delta = max(delta, od_prs[t].get_cur_delta());
		return delta;
	}
	bool fnd_routes()
	{
		bool fl = false;
		for (auto& now : od_prs)
			fl = (now.find_new_route() || fl);
		return fl;
	}
	void solve_flow()
	{
		while (fnd_routes())
			while (get_delta() > alpha)
				for (auto& now : od_prs)
				{
					now.calc_next_update();
					now.implement_next_update();
				}
		show_stat();
	}
};
