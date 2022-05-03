template <class type>
class matrix
{
private:
	int n, m;
	vector <vector <type>> values;
public:
	matrix(int k1, int k2, vector <vector <type>>& A): n(k1), m(k2), values(A)
	{	}
	matrix(int k1, int k2): n(k1), m(k2)
	{
		values = vector <vector <type>>(n, vector <type>(m, 0));
	}
	matrix() : n(0), m(0)
	{ }
	pair <int, int> size()
	{
		return { n, m };
	}
	vector <type>& operator[](int i)
	{
		return values[i];
	}
	matrix operator* (matrix& A)
	{
		vector <vector <type>> B(n);
		int d = A.size().second;
		for (int i = 0; i < n; i++)
		{
			B[i].resize(d, 0);
			for (int j = 0; j < d; j++)
				for (int k = 0; k < m; k++)
					B[i][j] += values[i][k] * A[k][j];
		}
		return matrix(n, d, B);
	}
	matrix operator* (double t)
	{
		vector <vector <type>> B(n);
		for (int i = 0; i < n; i++)
		{
			B[i].resize(m, 0);
			for (int j = 0; j < m; j++)
				B[i][j] = t * values[i][j];
		}
		return matrix(n, m, B);
	}
	matrix operator/ (double t)
	{
		vector <vector <type>> B(n);
		for (int i = 0; i < n; i++)
		{
			B[i].resize(m, 0);
			for (int j = 0; j < m; j++)
				B[i][j] = values[i][j] / t;
		}
		return matrix(n, m, B);
	}
	matrix operator+ (matrix A)
	{
		vector <vector <type>> B(n);
		for (int i = 0; i < n; i++)
		{
			B[i].resize(m, 0);
			for (int j = 0; j < m; j++)
				B[i][j] = values[i][j] + A[i][j];
		}
		return matrix(n, m, B);
	}
	matrix operator- (matrix A)
	{
		vector <vector <type>> B(n);
		for (int i = 0; i < n; i++)
		{
			B[i].resize(m, 0);
			for (int j = 0; j < m; j++)
				B[i][j] = values[i][j] - A[i][j];
		}
		return matrix(n, m, B);
	}
	matrix trans()
	{
		vector <vector <type>> B(m);
		for (int i = 0; i < m; i++)
		{
			B[i].resize(n);
			for (int j = 0; j < n; j++)
				B[i][j] = values[j][i];
		}
		return matrix(m, n, B);
	}
	type norm2()
	{
		type temp = 0;
		for (int i = 0; i < n; i++)
			temp += values[i][0] * values[i][0];
		return sqrt(temp);
	}
	matrix <type> get_first_column()
	{
		vector <vector <type>> T(n, { 0 });
		for (int i = 0; i < n; i++)
			T[i][0] = values[i][0];
		return matrix(n, 1, T);
	}
	double det(vector <vector <double>> A = {})
	{
		int a_size = A.size();
		if (a_size == 0)
		{
			A = values;
			a_size = n;
		}
		double ans = 1;
		for (int i = 0; i < a_size; i++)
		{
			int t = i;
			while (t < a_size && A[i][t] == 0)
				t++;
			if (t == a_size)
				return 0;
			if (i != t)
				ans *= -1;
			swap(A[i], A[t]);
			ans *= A[i][i];
			for (int j = i + 1; j < a_size; j++)
			{
				double temp = A[j][i] / A[i][i];
				for (int k = i; k < a_size; k++)
					A[j][k] -= A[i][k] * temp;
			}
		}
		return ans;
	}
	double alg_comp(int ti, int tj)
	{
		vector <vector <double>> B = values;
		int t = B.size();
		t = values.size();
		vector <vector <double>> C(t - 1, vector <double>(t - 1));
		for (int i = 0; i < t - 1; i++)
			for (int j = 0; j < t - 1; j++)
				C[i][j] = B[i + (i >= ti)][j + (j >= tj)];
		return det(C) * pow(-1, ti + tj);
	}
	void show()
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				cout << values[i][j] << ' ';
			cout << '\n';
		}
		cout << '\n' << '\n';
	}
	static matrix <double> create_e(int n, int m = 1)
	{
		vector <vector <double>> A(n);
		for (int i = 0; i < n; i++)
		{
			A[i].resize(m, 0);
			if (i < m)
				A[i][i] = 1;
		}
		return matrix <double>(n, m, A);
	}
	matrix <double> inv()
	{
		matrix <double> A(n, m, values);
		matrix <double> B = matrix<double>::create_e(n, n);
		for (int i = 0; i < n; i++)
		{
			int t = i;
			while (A[t][i] == 0)
				t++;
			swap(A[i], A[t]);
			swap(B[i], B[t]);
			double temp;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;
				temp = A[j][i] / A[i][i];
				for (int k = 0; k < n; k++)
				{
					A[j][k] -= A[i][k] * temp;
					B[j][k] -= B[i][k] * temp;
				}
			}
			temp = A[i][i];
			for (int j = 0; j < n; j++)
			{
				A[i][j] /= temp;
				B[i][j] /= temp;
			}
		}
		return B;
	}
};