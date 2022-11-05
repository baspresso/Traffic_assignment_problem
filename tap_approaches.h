class tap
{
protected:
	vector <int> init_node, term_node, link_type;
	vector <double> capacity, length, free_flow_time, b, power, speed, toll;
	vector <double> current_flow;
	vector <tuple <int, int, double>> od_pairs;
	vector <vector <int>> adj_list;
	int number_of_od, number_of_edges, number_of_nodes;
	string test;
public:
	bool isNumber(char c)
	{
		return c >= '0' && c <= '9';
	}
	double get_value(int& i, string& line)
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
	void get_edge(string& line)
	{
		int i = 0;
		number_of_edges++;
		init_node.push_back(get_value_int(i, line));
		init_node[number_of_edges - 1]--;
		term_node.push_back(get_value_int(i, line));
		term_node[number_of_edges - 1]--;
		capacity.push_back(get_value(i, line));
		length.push_back(get_value(i, line));
		free_flow_time.push_back(get_value(i, line));
		b.push_back(get_value(i, line));
		power.push_back(get_value(i, line));
		speed.push_back(get_value(i, line));
		toll.push_back(get_value(i, line));
		link_type.push_back(get_value_int(i, line));
	}
	void get_od_pair(ifstream& in)
	{

		int i = 0, origin, dest;
		double demand;
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
				od_pairs.push_back({ origin, dest, demand });
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
			get_edge(line);
		in = ifstream(test + "_trips.txt");
		getline(in, line);
		int i = 0;
		number_of_nodes = get_value_int(i, line);
		for (int j = 0; j < number_of_nodes; j++)
			get_od_pair(in);
	}
	tap(string tn) : number_of_edges(0), number_of_od(0), test(tn)
	{
		get_data();
		current_flow.resize(number_of_edges, 0);
		adj_list.resize(number_of_nodes);
		for (int i = 0; i < number_of_edges; i++)
			adj_list[init_node[i]].push_back(i);
	}
	double delay(int i, double extra_flow = 0)
	{
		double flow_on_edge = current_flow[i] + extra_flow;
		return free_flow_time[i] + 0.15 * free_flow_time[i] * pow(flow_on_edge, 4) / pow(capacity[i], 4);
	}
	double get_edges_delay(vector <int>& edges, double extra_flow = 0)
	{
		double ans = 0;
		for (auto& edge : edges)
			ans += delay(edge, extra_flow);
		return ans;
	}
	double delay_integ(int i, double flow_on_edge)
	{
		return free_flow_time[i] * (flow_on_edge + 0.15 * capacity[i] * pow(flow_on_edge / capacity[i], 5) / 5);
	}
	double objective_function()
	{
		double ans = 0, flow_on_edge;
		for (int i = 0; i < number_of_edges; i++)
			ans += delay_integ(i, current_flow[i]);
		return ans;
	}
	double get_best_answer()
	{
		double ans = 0, flow_on_edge;
		string line;
		ifstream in(test + "_flow.txt");
		int j = 0;
		for (int i = 0; i < number_of_edges; i++)
		{
			getline(in, line);
			j = 0;
			get_value(j, line); get_value(j, line);
			flow_on_edge = get_value(j, line);
			ans += delay_integ(i, flow_on_edge);
		}
		return ans;
	}
	void show_stat()
	{
		cout << "best ans: " << get_best_answer() << '\n';
		cout << "obj function: " << objective_function() << '\n';
		for (int i = 0; i < current_flow.size(); i++)
			cout << "from: " << init_node[i] + 1 << " to " << term_node[i] + 1 << " volume " << current_flow[i] << " cost " <<
			delay(i) << '\n';
	}
	virtual void solve_flow() = 0;
};

class tap_od_equilibrium : public tap
{
protected:
	vector <vector <vector <int>>> od_routes;
	vector <double> od_delays;
public:
	tap_od_equilibrium(string tn) : tap(tn)
	{
		od_routes.resize(number_of_od);
		od_delays.resize(number_of_od);
	}
	double get_od_delta(int t)
	{
		double del_min = 1e9, del_max = 0, del_route;
		if (od_routes[t].size() > 0)
			del_min = get_edges_delay(od_routes[t][0]);
		for (int i = 0; i < od_routes[t].size(); i++)
		{
			del_route = get_edges_delay(od_routes[t][i]);
			del_min = min(del_min, del_route);
			del_max = max(del_max, del_route);
		}
		return del_max - del_min;
	}
	virtual void balance_od_pair(int t) = 0;
	double get_delta()
	{
		double delta = 0;
		for (int t = 0; t < number_of_od; t++)
			delta = max(delta, get_od_delta(t));
		return delta;
	}
	void update_od_delays()
	{
		double delay;
		for (int t = 0; t < number_of_od; t++)
		{
			if (od_routes[t].size() > 0)
				delay = get_edges_delay(od_routes[t][0]);
			for (auto& route : od_routes[t])
				delay = min(delay, get_edges_delay(route));
			od_delays[t] = delay;
		}
	}
	bool find_new_route(int t)
	{
		priority_queue <pair <double, int>, vector <pair <double, int>>, greater <pair <double, int>>> q;
		int origin = get<0>(od_pairs[t]), dest = get<1>(od_pairs[t]), u;
		q.push({ origin, -1 });
		set <int> processed;
		map <int, int> used_edge;
		double temp;
		while (q.top().second == -1 || term_node[q.top().second] != dest)
		{
			if (q.top().second != -1)
				u = term_node[q.top().second];
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
			used_edge[u] = q.top().second;
			q.pop();
			for (auto edge : adj_list[u])
				if (processed.find(term_node[edge]) == processed.end())
					q.push({ temp + delay(edge), edge });
		}
		used_edge[dest] = q.top().second;
		double delay_cur_route = q.top().first;
		if (delay_cur_route < od_delays[t] || od_routes[t].size() == 0)
		{
			int now = dest;
			vector <int> new_route;
			while (now != origin)
			{
				new_route.push_back(used_edge[now]);
				now = init_node[used_edge[now]];
			}
			reverse(new_route.begin(), new_route.end());
			
			od_routes[t].push_back(new_route);
			if (od_routes[t].size() == 1)
				for (auto edge : od_routes[t][0])
					current_flow[edge] += get<2>(od_pairs[t]);
			return true;
		}
		else
			return false;
	}
	virtual void solve_flow() = 0;
};
