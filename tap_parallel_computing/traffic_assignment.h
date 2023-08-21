namespace traffic_assignment {
  template <class T>
  struct Link {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  public:
    int init, term, type;
    T capacity, length, free_flow_time, b, power, speed, toll, flow;
    Link(int init, int term, int type, T capacity, T length,
      T free_flow_time, T b, T power, T speed, T toll) :
      init(init), term(term), type(type), capacity(capacity), length(length),
      free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll) { };
    T Delay(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      return free_flow_time + 0.15 * free_flow_time * pow(temp_flow, 4) / pow(capacity, 4);
    }
    T DelayInteg(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      return free_flow_time * (temp_flow + 0.15 * capacity * pow(temp_flow / capacity, 5) / 5);
    }
    T DelayDer(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      return free_flow_time * 0.15 * 4 * pow(temp_flow / capacity, 3) / capacity;
    }
    T DelaySecondDer(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      return free_flow_time * 0.15 * 12 * pow(temp_flow, 2) / pow(capacity, 4);
    }
    static T GetLinksDelay(vector <Link <T>>& links, const vector <int>& links_list) {
      T ans = 0;
      for (auto now : links_list)
        ans += links[now].Delay();
      return ans;
    }
  };

  template <class T>
  class OriginDestinationPair {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  private:
    int origin_, dest_;
    T demand_, current_delay_, next_delay_, next_delta_;
    vector <Link <T>>& links_;
    const vector <vector <int>>& adjacency_list_;
    vector <vector <int>> current_routes_, next_routes_;
    unordered_map <int, T> current_links_flow_, next_links_flow_;
    unordered_map <int, T> fixed_links_flow_;
    const T eps_ = 1e-6;
    MatrixXd CreateG(MatrixXd& temp_flow, unordered_map <int, int>& mapping_links_inv) {
      int n = temp_flow.rows();
      MatrixXd g(n, n);
      g = MatrixXd::Identity(n, n);
      for (int i = 0; i < n; i++)
        g(i, i) = 1 / links_[mapping_links_inv[i]].DelayDer(fixed_links_flow_[mapping_links_inv[i]] + temp_flow(i, 0));
      return g;
    }
    MatrixXd CreateDelayColumn(MatrixXd& temp_flow, unordered_map <int, int>& M_links_inv) {
      MatrixXd delay_column(temp_flow.rows(), 1);
      for (int i = 0; i < temp_flow.rows(); i++)
        delay_column(i, 0) = links_[M_links_inv[i]].Delay(fixed_links_flow_[M_links_inv[i]] + temp_flow(i, 0));
      return delay_column;
    }
    tuple < unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, MatrixXd > Compression(const vector <vector <int>>& routes) {
      int count = 0;
      set <int> actually_used_links, actually_used_nodes;
      unordered_map <int, int> mapping_nodes, mapping_nodes_inv, mapping_edges, mapping_edges_inv;
      for (int i = 0; i < routes.size(); i++) {
        for (auto now : routes[i])
          if (actually_used_links.find(now) == actually_used_links.end()) {
            actually_used_links.insert(now);
            if (actually_used_nodes.find(links_[now].init) == actually_used_nodes.end()) {
              actually_used_nodes.insert(links_[now].init);
              mapping_nodes[links_[now].init] = actually_used_nodes.size() - 1;
              mapping_nodes_inv[actually_used_nodes.size() - 1] = links_[now].init;
            }
            if (actually_used_nodes.find(links_[now].term) == actually_used_nodes.end()) {
              actually_used_nodes.insert(links_[now].term);
              mapping_nodes[links_[now].term] = actually_used_nodes.size() - 1;
              mapping_nodes_inv[actually_used_nodes.size() - 1] = links_[now].term;
            }
          }
      }
      int n = actually_used_nodes.size() - 1, m = actually_used_links.size();
      for (auto now : actually_used_nodes)
        if (mapping_nodes[now] == n)
          swap(mapping_nodes[origin_], mapping_nodes[now]);
      MatrixXd compressed_adjacency_matrix = MatrixXd::Zero(n, m);
      count = 0;
      for (auto now : actually_used_links) {
        mapping_edges[now] = count;
        mapping_edges_inv[count] = now;
        if (links_[now].init != origin_)
          compressed_adjacency_matrix(mapping_nodes[links_[now].init], count) = 1;
        if (links_[now].term != origin_)
          compressed_adjacency_matrix(mapping_nodes[links_[now].term], count) = -1;
        count++;
      }
      return { mapping_nodes, mapping_nodes_inv, mapping_edges, mapping_edges_inv, compressed_adjacency_matrix };
    }
    void UpdateNextData() {
      T delay_min = 1e9, delay_max = 0, delay_route;
      if (next_routes_.size() > 0) {
        delay_route = 0;
        for (auto now : next_routes_[0])
          delay_route += links_[now].Delay(fixed_links_flow_[now] + next_links_flow_[now]);
        delay_min = delay_route;
        delay_max = delay_route;
      }
      for (int i = 1; i < next_routes_.size(); i++) {
        delay_route = 0;
        for (auto now : next_routes_[i])
          delay_route += links_[now].Delay(fixed_links_flow_[now] + next_links_flow_[now]);
        delay_min = min(delay_min, delay_route);
        delay_max = max(delay_max, delay_route);
      }
      next_delay_ = delay_min;
      next_delta_ = delay_max - delay_min;
    }
    void BalanceNextRoutes() {
      tuple <unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, MatrixXd > comp = Compression(next_routes_);
      unordered_map <int, int> mapping_nodes = get<0>(comp), mapping_nodes_inv = get<1>(comp), mapping_edges = get<2>(comp), mapping_edges_inv = get<3>(comp);
      MatrixXd compressed_adjacency_matrix = get<4>(comp);
      int n = compressed_adjacency_matrix.rows(), m = compressed_adjacency_matrix.cols();
      MatrixXd temp_flow_now(m, 1), temp_flow_prev(m, 1);
      for (int i = 0; i < next_routes_.size(); i++)
        for (auto& now : next_routes_[i])
          temp_flow_now(mapping_edges[now], 0) += demand_ / next_routes_.size();
      next_delta_ = eps_ + 1;
      next_links_flow_.clear();
      for (int i = 0; i < m; i++)
        next_links_flow_[mapping_edges_inv[i]] = temp_flow_now(i, 0);
      UpdateNextData();
      while (next_delta_ > eps_) {
        temp_flow_prev = temp_flow_now;
        MatrixXd g = CreateG(temp_flow_now, mapping_edges_inv);
        temp_flow_now -= g * (MatrixXd::Identity(m, m) - compressed_adjacency_matrix.transpose() *
          (compressed_adjacency_matrix * g * compressed_adjacency_matrix.transpose()).inverse() *
          compressed_adjacency_matrix * g) * CreateDelayColumn(temp_flow_now, mapping_edges_inv);
        next_links_flow_.clear();
        for (int i = 0; i < m; i++)
          next_links_flow_[mapping_edges_inv[i]] = temp_flow_now(i, 0);
        UpdateNextData();
      }
    }
    void UpdateCurrentDelay() {
      T delay_min = 1e9, delay_route;
      if (current_routes_.size() > 0) {
        delay_route = Link<T>::GetLinksDelay(links_, current_routes_[0]);
        delay_min = delay_route;
      }
      for (int i = 1; i < current_routes_.size(); i++) {
        delay_route = Link<T>::GetLinksDelay(links_, current_routes_[i]);
        delay_min = min(delay_min, delay_route);
      }
      current_delay_ = delay_min;
    }
    void GetFixedFlow(mutex& mtx) {
      fixed_links_flow_.clear();
      mtx.lock();
      //auto start = chrono::high_resolution_clock::now();
      for (int i = 0; i < current_routes_.size(); i++)
        for (auto now : current_routes_[i])
          fixed_links_flow_[now] = links_[now].flow - current_links_flow_[now];
      //auto end = chrono::high_resolution_clock::now();
      //chrono::duration<float> duration = end - start;
      //cout << "Duration of calc_next " << duration.count() << " sec \n";
      mtx.unlock();
    }
  public:
    OriginDestinationPair(int& origin, int& dest, T& demand, vector <Link <T>>& links, vector <vector <int>>& adjacency_list) :
      origin_(origin), dest_(dest), demand_(demand), links_(links), adjacency_list_(adjacency_list), current_delay_(0) { };
    vector <int> BestRoute() {
      UpdateCurrentDelay();
      //cout << current_delay_ << '\n';
      priority_queue <pair <T, int>, vector <pair <T, int>>, greater <pair <T, int>>> q;
      int u;
      q.push({ origin_, -1 });
      unordered_set <int> processed;
      unordered_map <int, int> used_link;
      T temp;
      while (q.top().second == -1 || links_[q.top().second].term != dest_) {
        if (q.top().second != -1)
          u = links_[q.top().second].term;
        else
          u = origin_;
        if (processed.find(u) != processed.end()) {
          q.pop();
          continue;
        }
        temp = q.top().first;
        if (processed.find(u) != processed.end())
          continue;
        processed.insert(u);
        used_link[u] = q.top().second;
        q.pop();
        for (auto now : adjacency_list_[u])
          if (processed.find(links_[now].term) == processed.end())
            q.push({ temp + links_[now].Delay(), now });
      }
      used_link[dest_] = q.top().second;
      int now = dest_;
      vector <int> new_route;
      while (now != origin_) {
        new_route.push_back(used_link[now]);
        now = links_[used_link[now]].init;
      }
      reverse(new_route.begin(), new_route.end());
      return new_route;
    }
    bool FindNewRoute() {
      vector <int> new_route = BestRoute();
      if (Link<T>::GetLinksDelay(links_, new_route) < current_delay_ || current_routes_.size() == 0) {
        current_delay_ = Link<T>::GetLinksDelay(links_, new_route);
        current_routes_.push_back(new_route);
        if (current_routes_.size() == 1)
          for (auto now : current_routes_[0]) {
            links_[now].flow += demand_;
            current_links_flow_[now] = demand_;
          }
        return true;
      }
      return false;
    }
    T GetCurrentDelta() {
      T delay_min = 1e9, delay_max = 0, delay_route;
      if (current_routes_.size() > 0) {
        delay_route = Link<T>::GetLinksDelay(links_, current_routes_[0]);
        delay_min = delay_route;
        delay_max = delay_route;
      }
      for (int i = 1; i < current_routes_.size(); i++) {
        delay_route = Link<T>::GetLinksDelay(links_, current_routes_[i]);
        delay_min = min(delay_min, delay_route);
        delay_max = max(delay_max, delay_route);
      }
      return abs(delay_max - delay_min);
    }
    void CalculateNextUpdate(mutex& mtx) {
      //auto start = chrono::high_resolution_clock::now();
      vector <vector <int>> possible_routes = current_routes_;
      next_routes_ = {};
      next_delay_ = -1;
      next_links_flow_ = current_links_flow_;
      GetFixedFlow(mtx);
      while (!possible_routes.empty()) {
        int best_route = 0;
        T delay_best = -1, delay_now;
        for (int i = 0; i < possible_routes.size(); i++) {
          delay_now = 0;
          for (auto now : possible_routes[i])
            delay_now += links_[now].Delay(fixed_links_flow_[now] + next_links_flow_[now]);
          if ((delay_now < delay_best) || (delay_best == -1)) {
            best_route = i;
            delay_best = delay_now;
          }
        }
        if ((delay_best < next_delay_) || (next_delay_ == -1)) {
          next_routes_.push_back(possible_routes[best_route]);
          BalanceNextRoutes();
        }
        swap(possible_routes[best_route], possible_routes[possible_routes.size() - 1]);
        possible_routes.pop_back();
        //auto end = chrono::high_resolution_clock::now();
        //chrono::duration<float> duration = end - start;
        //cout << "Duration of calc_next " << duration.count() << " sec \n";
      }
    }
    bool ImplementNextUpdate(mutex& mtx) {
      //auto start = chrono::high_resolution_clock::now();
      unordered_set <int> all_links, cur_links, next_links;
      for (int i = 0; i < current_routes_.size(); i++)
        for (auto now : current_routes_[i]) {
          cur_links.insert(now);
          all_links.insert(now);
        }
      for (int i = 0; i < next_routes_.size(); i++)
        for (auto now : next_routes_[i]) {
          next_links.insert(now);
          all_links.insert(now);
        }
      T current_result = 0, next_result = 0;
      bool result;
      mtx.lock();
      for (auto now : all_links)
        current_result += links_[now].DelayInteg();
      for (auto now : all_links)
        next_result += links_[now].DelayInteg(fixed_links_flow_[now] + next_links_flow_[now]);
      //cout << "diff: " << current_result - next_result << '\n';
      if (current_result >= next_result) {
        for (auto now : cur_links)
          links_[now].flow -= current_links_flow_[now];
        for (auto now : next_links)
          links_[now].flow += next_links_flow_[now];
        current_links_flow_ = next_links_flow_;
        current_routes_ = next_routes_;
        current_delay_ = next_delay_;
        result = true;
      }
      else
        result = false;
      mtx.unlock();
      return result;
    }
    int RoutesCount() {
      return current_routes_.size();
    }
    T GetDemand() {
      return demand_;
    }
  };

  template <typename T>
  class TrafficAssignmentApproach {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  protected:
    T best_ans_;
    vector <Link <T>> links_;
    vector <vector <int>> adjacency_list_;                    // Links for each node
    vector <OriginDestinationPair <T>> origin_destination_pairs_;
    int number_of_origin_destination_pairs_, number_of_links_, number_of_nodes_, number_of_zones_;
    string test_name_;
    ofstream time_out, rgap_out, func_out;
    std::chrono::time_point<struct std::chrono::steady_clock, class std::chrono::duration<__int64, struct std::ratio<1, 1000000000> > > start_;
    std::chrono::duration<float> time_on_statistics_;
    bool IsNumber(char c) {
      return c >= '0' && c <= '9';
    }
    T GetValue(int& i, const string& line) {
      while (!IsNumber(line[i]))
        i++;
      string temp1 = "";
      while (IsNumber(line[i]))
        temp1 += line[i++];
      if (line[i++] != '.')
        return stoi(temp1);
      string temp2 = "";
      while (IsNumber(line[i]))
        temp2 += line[i++];
      return stold(temp1) + stold(temp2) * (pow(0.1, temp2.size()));
    }
    int GetValueInt(int& i, const string& line) {
      while (!IsNumber(line[i]))
        i++;
      string temp1 = "";
      while (IsNumber(line[i]))
        temp1 += line[i++];
      if (line[i++] != '.')
        return stoi(temp1);
      string temp2 = "";
      while (IsNumber(line[i]))
        temp2 += line[i++];
      return stold(temp1);
    }
    void GetLink(const string& line) {
      int i = 0;
      number_of_links_++;
      int init, term, type;
      T capacity, length, free_flow_time, b, power, speed, toll, flow;
      init = GetValueInt(i, line) - 1; term = GetValueInt(i, line) - 1;
      capacity = GetValue(i, line);  length = GetValue(i, line);
      free_flow_time = GetValue(i, line); b = GetValue(i, line);
      power = GetValue(i, line);  speed = GetValue(i, line);
      toll = GetValue(i, line); type = GetValueInt(i, line);
      links_.push_back(Link<T>(init, term, type, capacity, length, free_flow_time, b, power, speed, toll));
    }
    void GetOriginDestinationPair(ifstream& in) {
      int i = 0, origin, dest;
      T demand;
      string line;
      getline(in, line);
      origin = GetValueInt(i, line) - 1;
      getline(in, line);
      i = 0;
      int count = 0;
      while (count < number_of_zones_ - 1) {
        dest = GetValueInt(i, line) - 1;
        if (dest != origin)
          count++;
        demand = GetValue(i, line);
        i++;
        if (demand > 0) {
          number_of_origin_destination_pairs_++;
          //cout << demand_ << '\n';
          origin_destination_pairs_.push_back(OriginDestinationPair <T>(origin, dest, demand, links_, adjacency_list_));
        }
        if (i >= line.size() - 1) {
          getline(in, line);
          i = 0;
        }
      }
    }
    void GetData() {
      string line;
      ifstream in(test_name_ + "_net.txt");
      while (getline(in, line))
        GetLink(line);
      in = ifstream(test_name_ + "_trips.txt");
      getline(in, line);
      int i = 0;
      number_of_zones_ = GetValueInt(i, line);
      getline(in, line);
      i = 0;
      number_of_nodes_ = GetValueInt(i, line);
      for (int j = 0; j < number_of_zones_; j++)
        GetOriginDestinationPair(in);
      number_of_origin_destination_pairs_ = origin_destination_pairs_.size();
    }
    T RelativeGap() {
      T result = 1, numerator = 0, denominator = 0;
      for (auto& od_pair : origin_destination_pairs_) {
        vector <int> best_route = od_pair.BestRoute();
        numerator += od_pair.GetDemand() * Link<T>::GetLinksDelay(links_, best_route);
      }
      for (auto& link : links_)
        denominator += link.flow * link.Delay();
      result -= numerator / denominator;
      return result;
    }
    T GetBestAnswer() {
      T ans = 0, flow_on_link;
      string line;
      ifstream in(test_name_ + "_flow.txt");
      int j = 0;
      for (int i = 0; i < number_of_links_; i++) {
        getline(in, line);
        j = 0;
        GetValue(j, line); GetValue(j, line);
        flow_on_link = GetValue(j, line);
        ans += links_[i].DelayInteg(flow_on_link);
      }
      return ans;
    }
  public:
    TrafficAssignmentApproach(string test_name) : number_of_links_(0), number_of_origin_destination_pairs_(0), test_name_(test_name), time_on_statistics_(0) {
      start_ = std::chrono::high_resolution_clock::now();
      time_out.open("time_out.txt");
      func_out.open("func_out.txt");
      rgap_out.open("rgap_out.txt");
      time_out.close();
      func_out.close();
      rgap_out.close();
      GetData();
      adjacency_list_.resize(number_of_nodes_);
      for (int i = 0; i < number_of_links_; i++)
        adjacency_list_[links_[i].init].push_back(i);
      best_ans_ = GetBestAnswer();
    }
    virtual ~TrafficAssignmentApproach() {}
    T ObjectiveFunction() {
      T ans = 0;
      for (int i = 0; i < number_of_links_; i++)
        ans += links_[i].DelayInteg();
      return ans;
    }
    void GetStatistics() {
      auto statistics_start = std::chrono::high_resolution_clock::now();
      time_out.open("time_out.txt", std::ios::app);
      func_out.open("func_out.txt", std::ios::app);
      rgap_out.open("rgap_out.txt", std::ios::app);
      rgap_out << setprecision(15) << RelativeGap() << '\n';
      func_out << setprecision(15) << ObjectiveFunction() << '\n';
      time_on_statistics_ += std::chrono::high_resolution_clock::now() - statistics_start;
      time_out << (chrono::duration_cast<chrono::seconds>(std::chrono::high_resolution_clock::now() - start_ - time_on_statistics_)).count() << '\n';
      time_out.close();
      func_out.close();
      rgap_out.close();
      //for (int i = 0; i < number_of_links_; i++)
      //  cout << "from: " << links[i].init + 1 << " to " << links[i].term + 1 << " volume " << links[i].flow << " cost " <<
      //  links[i].Delay() << '\n';
    }
    virtual void SolveFlow() {};
  };

  template <typename T>                    // value type
  class RouteBasedApproach : public TrafficAssignmentApproach <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  private:
    queue <int> ready_to_implement_, tasks_;
    string test_name_;
    T alpha_ = 1e-4;
    mutex mtx_, mtx_tasks_;
    T GetDelta() {
      T delta = 0;
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++)
        delta = max(delta, TrafficAssignmentApproach<T>::origin_destination_pairs_[t].GetCurrentDelta());
      //cout << delta << '\n';
      return delta;
    }
    bool NewRouteFound() {
      bool fl = false;
      for (auto& now : TrafficAssignmentApproach<T>::origin_destination_pairs_)
        fl = (now.FindNewRoute() || fl);
      return fl;
    }
    void ThreadCalculation(const bool& finished) {
      while (!finished) {
        if (tasks_.empty())
          this_thread::sleep_for(100ms);
        int j = -1;
        mtx_tasks_.lock();
        if (!tasks_.empty()) {
          j = tasks_.front();
          tasks_.pop();
        }
        mtx_tasks_.unlock();
        if (j == -1)
          continue;
        TrafficAssignmentApproach<T>::origin_destination_pairs_[j].CalculateNextUpdate(mtx_);
        TrafficAssignmentApproach<T>::origin_destination_pairs_[j].ImplementNextUpdate(mtx_);
      }
    }
  public:
    RouteBasedApproach(string test_name) : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(test_name) {
    }
    void SolveFlow() override {
      int num_threads = thread::hardware_concurrency() - 1;
      num_threads = 1;
      vector <thread> threads_list(num_threads);
      vector <int> tasks_order(TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_);
      vector <int> now_in_process(num_threads, 0);
      for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++)
        tasks_order[i] = i;
      int count = 0;
      bool finished = false;
      for (int i = 0; i < num_threads; i++)
        threads_list[i] = thread(&RouteBasedApproach<T>::ThreadCalculation, this, ref(finished));
      while (count++ < 20) {
        NewRouteFound();
        auto start = chrono::high_resolution_clock::now();
        while (GetDelta() > alpha_) {
          random_shuffle(tasks_order.begin(), tasks_order.end());
          mtx_tasks_.lock();
          for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++)
            if (TrafficAssignmentApproach<T>::origin_destination_pairs_[tasks_order[i]].RoutesCount() > 1)
              tasks_.push(tasks_order[i]);
          mtx_tasks_.unlock();
          while (!tasks_.empty())
            this_thread::sleep_for(200ms);
        }
        TrafficAssignmentApproach<T>::GetStatistics();
        //ShowStatistics();
      }
      finished = true;
      for_each(threads_list.begin(), threads_list.end(), mem_fn(&thread::join));
      //TrafficAssignmentApproach<T>::ShowStatistics();
    }
  };

  template <typename T>                    // value type
  class LinkBasedApproach : public TrafficAssignmentApproach <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  private:
    MatrixXd flow_;
    const int number_of_iterations_ = 1e2;
    const T delta = 0.2;
    void UpdateFlowInfo() {
      for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_links_; i++)
        TrafficAssignmentApproach<T>::links_[i].flow = flow_(i, 0);
    }
    // use only ChangeFlow in order to change flow_
    void ChangeFlow(MatrixXd temp) {
      flow_ += temp;
      UpdateFlowInfo();
    }
    MatrixXd AllOrNothingSolution() {
      MatrixXd temp(TrafficAssignmentApproach<T>::number_of_links_, 1);
      for (auto& od_pair : TrafficAssignmentApproach<T>::origin_destination_pairs_) {
        vector <int> route = od_pair.BestRoute();
        for (auto now : route)
          temp(now, 0) += od_pair.GetDemand();
      }
      return temp;
    }
    MatrixXd FrankWolfDirection() {
      return AllOrNothingSolution() - flow_;
    }
    MatrixXd GenerateHessian() {
      MatrixXd result(TrafficAssignmentApproach<T>::number_of_links_, TrafficAssignmentApproach<T>::number_of_links_);
      for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_links_; i++)
        result(i, i) = TrafficAssignmentApproach<T>::links_[i].DelaySecondDer();
      return result;
    }
    T CalculateN(const MatrixXd& direction_dash, const MatrixXd& hessian_matrix, const MatrixXd& all_or_nothing) {
      return (direction_dash.transpose() * hessian_matrix * all_or_nothing)(0, 0);
    }
    T CalculateD(const MatrixXd& direction_dash, const MatrixXd& hessian_matrix, const MatrixXd& all_or_nothing) {
      return (direction_dash.transpose() * hessian_matrix * (all_or_nothing - direction_dash))(0, 0);
    }
    T CalculateAlpha(const MatrixXd& direction_dash, const MatrixXd& hessian_matrix, const MatrixXd& all_or_nothing) {
      T n = CalculateN(direction_dash, hessian_matrix, all_or_nothing), d = CalculateD(direction_dash, hessian_matrix, all_or_nothing);
      if (d > 0) {
        if (n / d > 1 - delta)
          return 1 - delta;
        else
          return n / d;
      }
      else
        return 0;
      return (direction_dash.transpose() * hessian_matrix * (all_or_nothing - direction_dash))(0, 0);
    }
    MatrixXd GenerateS(const MatrixXd& previous_s, const MatrixXd& direction_dash, const MatrixXd& hessian_matrix, const MatrixXd& all_or_nothing) {
      T alpha = CalculateAlpha(direction_dash, hessian_matrix, all_or_nothing);
      return alpha * previous_s + (1 - alpha) * all_or_nothing;
    }
    MatrixXd ConjugateFrankWolfDirection(const MatrixXd& s) {
      return s - flow_;
    }
    bool TryToImplementUpdate(MatrixXd flow_update) {
      T cur_func = TrafficAssignmentApproach<T>::ObjectiveFunction();
      ChangeFlow(flow_update);
      T updated_func = TrafficAssignmentApproach<T>::ObjectiveFunction();
      if (cur_func < updated_func) {
        ChangeFlow(-flow_update);
        return false;
      }
      return true;
    }
  public:
    LinkBasedApproach(string test_name) : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(test_name) {
      flow_ = MatrixXd(TrafficAssignmentApproach<T>::number_of_links_, 1);
    }
    void SolveFlow() override {
      ChangeFlow(AllOrNothingSolution());
      T tau = 1;
      MatrixXd descent_direction = flow_, s = flow_;
      for (int cnt_iterations = 0; cnt_iterations < number_of_iterations_; cnt_iterations++) {
        s = GenerateS(s, (1 - tau) * descent_direction, GenerateHessian(), AllOrNothingSolution());
        descent_direction = ConjugateFrankWolfDirection(s);
        tau = 0;
        for (T k = 0.5, i = 0; i < 10; i++, k /= 2)
          if (TryToImplementUpdate(k * descent_direction))
            tau += k;
        TrafficAssignmentApproach<T>::ShowStatistics();
      }
      TrafficAssignmentApproach<T>::ShowStatistics();
    }
  };
}
