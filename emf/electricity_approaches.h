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
	v_type delay(v_type temp_flow)
	{
		if (temp_flow == -1)
			temp_flow = flow;
		return free_flow_time + 0.15 * free_flow_time * pow(temp_flow, 4) / pow(capacity, 4);
	}
	v_type delay()
	{
		v_type temp_flow = flow;
		return free_flow_time + 0.15 * free_flow_time * pow(temp_flow, 4) / pow(capacity, 4);
	}
	v_type delay_integ(v_type temp_flow = -1)
	{
		if (temp_flow == -1)
			temp_flow = flow;
		return free_flow_time * (temp_flow + 0.15 * capacity * pow(temp_flow / capacity, 5) / 5);
	}
	v_type delay_der(v_type temp_flow)
	{
		return free_flow_time * 0.15 * 4 * pow(temp_flow / capacity, 3) / capacity;
	}
	v_type delay_der()
	{
		v_type temp_flow = flow;	
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
class lin_equation
{
private:
	matrix <v_type> A, b;
	matrix <v_type> gaussian_elimination(matrix <v_type> A_mod, matrix <v_type> b_mod)
	{
		int n = A_mod.size().first;
		for (int i = 0; i < n; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				v_type t = A_mod[j][i] / A_mod[i][i];
				for (int k = i; k < n; k++)
					A_mod[j][k] -= A_mod[i][k] * t;
				b_mod[j][0] -= b_mod[i][0] * t;
			}
		}
		for (int i = n - 1; i >= 0; i--)
		{
			b_mod[i][0] /= A_mod[i][i];
			for (int j = 0; j < i; j++)
				b_mod[j][0] -= b_mod[i][0] * A_mod[j][i];
		}
		return b_mod;
	}
public:
	lin_equation(const matrix <v_type>& A, const matrix <v_type>& b) : A(A), b(b)
	{	};
	matrix <v_type> rnd_non_zero_solution()
	{
		int n = A.size().first, m = A.size().second;
		//A.show(); b.show();
		vector <matrix <v_type>> basis(n);
		vector <bool> engaged(n, false);
		vector <int> indx(n);
		set <int> S;
		for (int i = 0; i < m; i++)
		{
			matrix <v_type> now = A.get_column(i);
			for (int j = 0; j < n; j++)
				if (now[j][0] != 0)
					if (!engaged[j])
					{
						engaged[j] = true;
						basis[j] = now;
						indx[j] = i;
						break;
					}
					else
					{
						//(basis[j] / basis[j][j][0] * now[j][0]).show();
						now -= basis[j] / basis[j][j][0] * now[j][0];
						//now 
					}
		}
		for (auto now : indx)
			S.insert(now);
		int attempt_cnt = 1e2;
		while (attempt_cnt-- > 0)
		{
			matrix <v_type> A_mod(basis[indx[0]]), b_mod(b);
			for (int i = 1; i < n; i++)
				A_mod.concatenate(A.get_column(indx[i]));
			matrix <v_type> solution(m, 1);
			mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
			v_type limitation = 0;
			for (int i = 0; i < b.size().first; i++)
				limitation += abs(b[i][0]);
			for (int i = 0; i < m; i++)
				if (S.find(i) == S.end())
				{
					v_type a = rng(), b = rng();
					if (a > b)
						swap(a, b);
					solution[i][0] = a / b * limitation;
				}
			b_mod -= A * solution;
			matrix <v_type> temp = gaussian_elimination(A_mod, b_mod);
			bool fl = true;
			for (int i = 0; i < n; i++)
			{
				solution[indx[i]][0] = temp[i][0];
				if (temp[i][0] == 0)
					fl = false;
			}
			if (fl)
				return solution;

		}
		// need exception
		cout << '!';
		return matrix <v_type> (m, 1);

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
		matrix <v_type> A(number_of_nodes - 1, number_of_links);
		for (int i = 0; i < number_of_links; i++)
		{
			int temp;
			if (lnks[i].init != v)
			{
				temp = (lnks[i].init > v) ? -1 : 0;
				A[lnks[i].init + temp][i] = 1;
			}
			if (lnks[i].term != v)
			{
				temp = (lnks[i].term > v) ? -1 : 0;
				A[lnks[i].term + temp][i] = -1;
			}
		}
		//A.show();
		return A;
	}
	matrix <v_type> create_b(int v)
	{
		matrix <v_type> b(number_of_nodes - 1, 1);
		if (v < emf)
			b[emf - 1][0] = demand;
		else
			b[emf][0] = demand;
		return b;
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
	emf_problem(const string& test) : number_of_links(0), test(test)
	{
		get_data();
	}
	void solve_flow()
	{
		mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
		int v = rng() % number_of_nodes;
		while (v == emf)
			v = rng() % number_of_nodes;
		matrix <v_type> A = create_A(v);
		matrix <v_type> b = create_b(v);
		lin_equation <v_type> lin_sys(A, b);
		matrix <v_type> G, A_trans, E, g, T;
		matrix <v_type> flow_now = lin_sys.rnd_non_zero_solution();
		for (int i = 0; i < number_of_links; i++)
			lnks[i].flow = flow_now[i][0];
		v_type delta = eps + 1;
		while (delta > eps)
		{
			//flow_now.show();
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
		cout <<	"v = " << v + 1 << '\n';
		flow_now.show();
		for (int i = 0; i < number_of_links; i++)
			cout << lnks[i].delay() << '\n';
	}
};
