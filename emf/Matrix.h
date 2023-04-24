namespace Matrix {
template <class T>
class Matrix {
private:
  int first_dimension_, second_dimension_;
  std::vector <std::vector <T>> values_;
public:
  Matrix(int first_dimension, int second_dimension, const std::vector <std::vector <T>>& values) :
    first_dimension_(first_dimension), second_dimension_(second_dimension), values_(values) {  }
  Matrix(int first_dimension, int second_dimension) : first_dimension_(first_dimension), second_dimension_(second_dimension) {
    values_ = std::vector <std::vector <T>>(first_dimension_, std::vector <T>(second_dimension_, 0));
  }
  Matrix() : first_dimension_(0), second_dimension_(0) { }
  std::pair <int, int> Size() {
    return { first_dimension_, second_dimension_ };
  }
  std::vector <T>& operator[](int i) {
    return values_[i];
  }
  const std::vector <T>& operator[](int i) const {
    return values_[i];
  }
  Matrix <T> operator* (const Matrix <T>& a) {
    int a_second_dimension = a.Size().second;
    Matrix <T> result(first_dimension_, a_second_dimension);
    for (int i = 0; i < first_dimension_; i++)
      for (int j = 0; j < a_second_dimension; j++)
        for (int k = 0; k < second_dimension_; k++)
          result[i][j] += values_[i][k] * a[k][j];
    return result;
  }
  Matrix <T> operator* (const T& k) {
    Matrix <T> result(*this);
    for (int i = 0; i < first_dimension_; i++)
      for (int j = 0; j < second_dimension_; j++)
        result[i][j] *= k;
    return result;
  }
  Matrix <T> operator/ (const T& k) {
    Matrix <T> result(*this);
    for (int i = 0; i < first_dimension_; i++)
      for (int j = 0; j < second_dimension_; j++)
        result[i][j] /= k;
    return result;
  }
  Matrix <T> operator+ (const Matrix <T>& a) {
    Matrix <T> result(*this);
    for (int i = 0; i < first_dimension_; i++)
      for (int j = 0; j < second_dimension_; j++)
        result[i][j] += a[i][j];
    return result;
  }
  Matrix <T>& operator += (const Matrix <T>& a) {
    return *this = *this + a;
  }
  Matrix operator- (const Matrix <T>& a) {
    Matrix <T> result(*this);
    for (int i = 0; i < first_dimension_; i++)
      for (int j = 0; j < second_dimension_; j++)
        result[i][j] -= a[i][j];
    return result;
  }
  Matrix <T> operator-= (const Matrix <T>& a) {
    return *this = *this - a;
  }
  Matrix <T> Transposed() const {
    Matrix <T> result(second_dimension_, first_dimension_);
    for (int i = 0; i < second_dimension_; i++)
      for (int j = 0; j < first_dimension_; j++)
        result[i][j] = values_[j][i];
    return result;
  }
  T NormOneDimension() const {
    T norm_value = 0;
    for (int i = 0; i < first_dimension_; i++)
      norm_value += abs(values_[i][0]);
    return norm_value;
  }
  Matrix <T> GetColumn(int j) const {
    Matrix <T> result(first_dimension_, 1);
    for (int i = 0; i < first_dimension_; i++)
      result[i][0] = values_[i][j];
    return result;
  }
  T Det() const {
    std::vector <std::vector <T>> values_copy(values_);
    T ans = 1;
    for (int i = 0; i < first_dimension_; i++) {
      int t = i;
      while (t < first_dimension_ && values_copy[i][t] == 0)
        t++;
      if (t == first_dimension_)
        return 0;
      if (i != t)
        ans *= -1;
      swap(values_copy[i], values_copy[t]);
      ans *= values_copy[i][i];
      for (int j = i + 1; j < first_dimension_; j++) {
        T temp = values_copy[j][i] / values_copy[i][i];
        for (int k = i; k < first_dimension_; k++)
          values_copy[j][k] -= values_copy[i][k] * temp;
      }
    }
    return ans;
  }
  void Show() const {
    for (int i = 0; i < first_dimension_; i++) {
      for (int j = 0; j < second_dimension_; j++)
        std::cout << values_[i][j] << ' ';
      std::cout << '\n';
    }
  }
  Matrix <T> Inversed() const;
};

template <class T>
Matrix <T> CreateIdentityMatrix(const int dimension) {
  Matrix <T> result(dimension, dimension);
  for (int i = 0; i < dimension; i++)
      result[i][i] = 1;
  return result;
}

template <class T>
Matrix <T> Concatenate(const Matrix <T>& a, const Matrix <T>& b) {
  int first_dimension = a.Size().first, second_dimension = a.Size().second + b.Size().second(),
      a_second_dimension = a.Size().second, b_second_dimension = b.Size().second;
  Matrix <T> result(first_dimension, second_dimension);
  for (int i = 0; i < first_dimension; i++) {
    for (int j = 0; j < a_second_dimension; j++)
      result[i][j] = a[i][j];
    for (int j = 0; j < b_second_dimension; j++)
      result[i][a_second_dimension + j] = b[i][j];
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::Inversed() const {
  Matrix <T> a(*this);
  Matrix <T> b = CreateIdentityMatrix<T>(first_dimension_);
  for (int i = 0; i < first_dimension_; i++) {
    int t = i;
    while (a[t][i] == 0)
      t++;
    swap(a[i], a[t]);
    swap(a[i], a[t]);
    T temp;
    for (int j = 0; j < first_dimension_; j++) {
      if (j == i)
        continue;
      temp = a[j][i] / a[i][i];
      for (int k = 0; k < first_dimension_; k++) {
        a[j][k] -= a[i][k] * temp;
        b[j][k] -= b[i][k] * temp;
      }
    }
    temp = a[i][i];
    for (int j = 0; j < first_dimension_; j++) {
      a[i][j] /= temp;
      b[i][j] /= temp;
    }
  }
  return b;
}
}

template <class T>
Matrix::Matrix<T> operator* (const T k, const Matrix::Matrix <T>& a) {
  return a * k;
}
