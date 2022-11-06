// link struct implemented

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
};

template <typename v_type>										// value type
class tap
{
protected:
	vector <link <v_type>> lnks;
	vector <tuple <int, int, v_type>> od_pairs;
	vector <vector <int>> adj_list;							// links for each node
	int number_of_od, number_of_links, number_of_nodes;
	string test;
	vector <vector <vector <int>>> od_routes;
	vector <v_type> od_delays;
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
			get_link(line);
		in = ifstream(test + "_trips.txt");
		getline(in, line);
		int i = 0;
		number_of_nodes = get_value_int(i, line);
		for (int j = 0; j < number_of_nodes; j++)
			get_od_pair(in);
	}
	tap(string tn) : number_of_links(0), number_of_od(0), test(tn)
	{
		get_data();
		//lnks.resize(number_of_links);
		adj_list.resize(number_of_nodes);
		for (int i = 0; i < number_of_links; i++)
			adj_list[lnks[i].init].push_back(i);
		od_routes.resize(this->number_of_od);
		od_delays.resize(this->number_of_od);
	}
	v_type get_links_delay(vector <int>& lnk_list, v_type extra_flow = 0)
	{
		v_type ans = 0;
		for (auto now : lnk_list)
			ans += lnks[now].delay(extra_flow);
		return ans;
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
	v_type get_od_delta(int t)
	{
		v_type del_min = 1e9, del_max = 0, del_route;
		if (od_routes[t].size() > 0)
			del_min = this->get_links_delay(od_routes[t][0]);
		for (int i = 0; i < od_routes[t].size(); i++)
		{
			del_route = get_links_delay(od_routes[t][i]);
			//cout << del_route << '\n';
			del_min = min(del_min, del_route);
			del_max = max(del_max, del_route);
		}
		return del_max - del_min;
	}
	virtual void balance_od_pair(int t) = 0;
	v_type get_delta()
	{
		v_type delta = 0;
		for (int t = 0; t < number_of_od; t++)
			delta = max(delta, get_od_delta(t));
		return delta;
	}
	void update_od_delays()
	{
		v_type delay;
		for (int t = 0; t < number_of_od; t++)
		{
			if (od_routes[t].size() > 0)
				delay = get_links_delay(od_routes[t][0]);
			for (auto& route : od_routes[t])
				delay = min(delay, get_links_delay(route));
			od_delays[t] = delay;
		}
	}
	bool find_new_route(int t)
	{
		priority_queue <pair <v_type, int>, vector <pair <v_type, int>>, greater <pair <v_type, int>>> q;
		int origin = get<0>(od_pairs[t]), dest = get<1>(od_pairs[t]), u;
		q.push({ origin, -1 });
		unordered_set <int> processed;
		map <int, int> used_link;
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
					q.push({ temp + lnks[now].delay(), now});
		}
		used_link[dest] = q.top().second;
		v_type delay_cur_route = q.top().first;
		if (delay_cur_route < od_delays[t] || od_routes[t].size() == 0)
		{
			int now = dest;
			vector <int> new_route;
			while (now != origin)
			{
				new_route.push_back(used_link[now]);
				now = lnks[used_link[now]].init;
			}
			reverse(new_route.begin(), new_route.end());

			od_routes[t].push_back(new_route);
			if (od_routes[t].size() == 1)
				for (auto now : od_routes[t][0])
					lnks[now].flow += get<2>(od_pairs[t]);
			return true;
		}
		else
			return false;
	}
	virtual void solve_flow() = 0;
};


/*template <class v_type>
class tap_od_equilibrium : public tap <v_type>
{
protected:
public:
	tap_od_equilibrium(string tn) : tap(tn)
	{
		od_routes.resize(this->number_of_od);
		od_delays.resize(this->number_of_od);
	}
	v_type get_od_delta(int t)
	{
		v_type del_min = 1e9, del_max = 0, del_route;
		if (od_routes[t].size() > 0)
			del_min = this->get_links_delay(od_routes[t][0]);
		for (int i = 0; i < od_routes[t].size(); i++)
		{
			del_route = get_links_delay(od_routes[t][i]);
			del_min = min(del_min, del_route);
			del_max = max(del_max, del_route);
		}
		return del_max - del_min;
	}
	virtual void balance_od_pair(int t) = 0;
	v_type get_delta()
	{
		v_type delta = 0;
		for (int t = 0; t < number_of_od; t++)
			delta = max(delta, get_od_delta(t));
		return delta;
	}
	void update_od_delays()
	{
		v_type delay;
		for (int t = 0; t < number_of_od; t++)
		{
			if (od_routes[t].size() > 0)
				delay = get_links_delay(od_routes[t][0]);
			for (auto& route : od_routes[t])
				delay = min(delay, get_links_delay(route));
			od_delays[t] = delay;
		}
	}
	bool find_new_route(int t)
	{
		priority_queue <pair <v_type, int>, vector <pair <v_type, int>>, greater <pair <v_type, int>>> q;
		int origin = get<0>(od_pairs[t]), dest = get<1>(od_pairs[t]), u;
		q.push({ origin, -1 });
		unordered_set <int> processed;
		map <int, int> used_link;
		v_type temp;
		while (q.top().second == -1 || lnks[q.top().second].term_node != dest)
		{
			if (q.top().second != -1)
				u = lnks[q.top().second].term_node;
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
				if (processed.find(lnks[now].term_node) == processed.end())
					q.push({ temp + lnks[now].delay, now });
		}
		used_link[dest] = q.top().second;
		v_type delay_cur_route = q.top().first;
		if (delay_cur_route < od_delays[t] || od_routes[t].size() == 0)
		{
			int now = dest;
			vector <int> new_route;
			while (now != origin)
			{
				new_route.push_back(used_link[now]);
				now = lnks[used_link[now]].init;
			}
			reverse(new_route.begin(), new_route.end());
			
			od_routes[t].push_back(new_route);
			if (od_routes[t].size() == 1)
				for (auto now : od_routes[t][0])
					lnks[now].current_flow += get<2>(od_pairs[t]);
			return true;
		}
		else
			return false;
	}
	virtual void solve_flow() = 0;
};*/
