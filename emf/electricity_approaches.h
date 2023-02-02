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
	v_type delay(v_type temp_flow = -1)
	{
		if (temp_flow == -1)
			temp_flow = flow;
		return free_flow_time + 0.15 * free_flow_time * pow(temp_flow, 4) / pow(capacity, 4);
	}
	v_type delay_integ(v_type temp_flow = -1)
	{
		if (temp_flow == -1)
			temp_flow = flow;
		return free_flow_time * (temp_flow + 0.15 * capacity * pow(temp_flow / capacity, 5) / 5);
	}
	v_type delay_der(v_type temp_flow = -1)
	{
		if (temp_flow == -1)
			temp_flow = flow;
		return free_flow_time * 0.15 * 4 * pow(temp_flow / capacity, 3) / capacity;
	}
	static v_type get_links_delay(vector <link <v_type>>& lnks, vector <int>& lnk_list)
	{
		v_type ans = 0;
		for (auto now : lnk_list)
			ans += lnks[now].delay();
		return ans;
	}
};

template <class v_type>
class emf_problem
{
private:
	int emf, number_of_links, number_of_nodes, number_of_zones;
	const v_type eps = 1e-6;
	v_type demand;
	vector <link <v_type>> lnks;
	string test;
	matrix <v_type> create_G()
	{
		vector <vector <v_type>> Temp(number_of_links, vector <v_type>(number_of_links, 0));
		for (int i = 0; i < number_of_links; i++)
			Temp[i][i] = 1 / lnks[i].delay_der();
		return matrix <v_type>(number_of_links, number_of_links, Temp);
	}
	matrix <v_type> create_g()
	{
		matrix <v_type> g(number_of_links, 1);
		for (int i = 0; i < number_of_links; i++)
			g[i][0] = lnks[i].delay();
		return g;
	}
	matrix <v_type> create_A(int v)
	{
		matrix <v_type> A;
		A = matrix <v_type>(number_of_nodes - 1, number_of_links);
		for (int i = 0; i < number_of_links; i++)
		{
			int temp;
			if (lnks[i].init != v)
			{
				if (lnks[i].init > v)
					temp = -1;
				else
					temp = 0;
				A[lnks[i].init + temp][i] = 1;
			}
			if (lnks[i].term != v)
			{
				if (lnks[i].term > v)
					temp = -1;
				else
					temp = 0;
				A[lnks[i].term + temp][i] = -1;
			}
		}
		return A;
	}
	bool isNumber(char c)
	{
		return c >= '0' && c <= '9';
	}
	v_type get_value(int& i, const string& line)
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
	int get_value_int(int& i, const string& line)
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
	void get_link(const string& line)
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
	void get_data()
	{
		string line;
		ifstream in(test + "_net.txt");
		while (getline(in, line))
			get_link(line);
		in = ifstream(test + "_info.txt");
		getline(in, line);
		int i = 0;
		number_of_zones = get_value_int(i, line);
		getline(in, line);
		i = 0;
		number_of_nodes = get_value_int(i, line);
		i = 0;
		getline(in, line);
		emf = get_value_int(i, line) - 1;
		i = 0;
		getline(in, line);
		demand = get_value_int(i, line);
	}

public:
	emf_problem(string test) : number_of_links(0), test(test)
	{
		get_data();
	}
	void solve_flow()
	{
		int v;
		for (auto &now : lnks)
			if ((now.init == emf) || (now.term == emf))
			{
				if (now.init == emf)
				{
					now.flow = demand;
					v = now.term;
				}
				else
				{
					now.flow = -demand;
					v = now.init;
				}
				break;
			}
		matrix <v_type> A = create_A(v);
		matrix <v_type> G, A_trans, E, g, T;
		matrix <v_type> flow_now(number_of_links, 1);
		v_type delta = eps + 1;
		while (delta > eps)
		{
			for (int i = 0; i < number_of_links; i++)
				flow_now[i][0] = lnks[i].flow;
			G = create_G();
			A_trans = A.trans();
			E = matrix<v_type>::create_e(number_of_links, number_of_links);
			g = create_g();
			T = (A * G * A_trans).inv();
			T = A_trans * T * A * G;
			T = E - T;
			T = G * T * g;
			flow_now = flow_now - T;
			for (int i = 0; i < number_of_links; i++)
				lnks[i].flow = flow_now[i][0];
			delta = T.norm2();
		}
		flow_now.show();
	}
};