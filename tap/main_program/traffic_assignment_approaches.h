namespace traffic_assignment {
  #define MAX_ROUTE_DELAY 1e9

  template <class T>
  struct Link {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  public:
    int init, term, type;
    long double power, capacity, length, free_flow_time, b, speed, toll;
    T flow;
    Link(int init, int term, int type, T capacity, T length,
      T free_flow_time, T b, T power, T speed, T toll) :
      init(init), term(term), type(type), capacity(capacity), length(length),
      free_flow_time(free_flow_time), b(b), power(power), speed(speed), toll(toll), flow(0) { };
    T Delay(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      //cout << power << endl;
      return free_flow_time + b * free_flow_time * pow(temp_flow, power) / pow(capacity, power);
    }
    T DelayInteg(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      //cout << power << endl;
      return free_flow_time * (temp_flow + b * capacity * pow(temp_flow / capacity, power + 1) / (power + 1));
    }
    T DelayDer(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      //cout << power << endl;
      return free_flow_time * b * power * pow(temp_flow / capacity, power - 1) / capacity;
    }
    T DelaySecondDer(T temp_flow = -1) {
      if (temp_flow == -1)
        temp_flow = flow;
      //cout << power << endl;
      return free_flow_time * b * power * (power - 1) * pow(temp_flow, power - 2) / pow(capacity, power);
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
  protected:
    int origin_, dest_;
    T demand_, current_delay_, next_delay_, next_delta_;
    vector <Link <T>>& links_;
    const vector <vector <int>>& adjacency_list_;
    vector <vector <int>> routes_;
    unordered_map <int, T> links_flow_;
    const T eps_ = 1e-32;
    void ClearFlow() {
      unordered_set <int> used_links;
      for (const auto& route : routes_)
        for (const auto& link : route)
          used_links.insert(link);
      for (const auto& link : used_links) {
        links_[link].flow -= links_flow_[link];
        links_flow_[link] = 0;
      }
    }
    void UpdateFlow(const MatrixXd& temp_flow, const unordered_map <int, int>& mapping_links_inv) {
      int n = temp_flow.rows();
      for (int i = 0; i < n; i++) {
        links_flow_[mapping_links_inv.at(i)] = temp_flow(i, 0);
        links_[mapping_links_inv.at(i)].flow += temp_flow(i, 0);
      }
    }
    void UpdateRoutesFlow(const MatrixXd& routes_flow) {
      for (int i = 0; i < routes_.size(); i++)
        for (auto now : routes_[i]) {
          links_flow_[now] += routes_flow(i, 0);
          links_[now].flow += routes_flow(i, 0);
        }
    }
    T GetDelay() {
      T delay = 0;
      for (const auto& route : routes_)
        delay = max(delay, Link<T>::GetLinksDelay(links_, route));
      return delay;
    }
    MatrixXd CreateG(const MatrixXd& temp_flow, const unordered_map <int, int>& mapping_links_inv) {
      int n = temp_flow.rows();
      MatrixXd g(n, n);
      g = MatrixXd::Identity(n, n);
      for (int i = 0; i < n; i++)
        g(i, i) = 1 / links_[mapping_links_inv.at(i)].DelayDer(links_[mapping_links_inv.at(i)].flow + temp_flow(i, 0));
      return g;
    }
    MatrixXd CreateDelayColumn(const MatrixXd& temp_flow, const unordered_map <int, int>& mapping_links_inv) {
      MatrixXd delay_column(temp_flow.rows(), 1);
      for (int i = 0; i < temp_flow.rows(); i++)
        delay_column(i, 0) = links_[mapping_links_inv.at(i)].Delay(links_[mapping_links_inv.at(i)].flow + temp_flow(i, 0));
      return delay_column;
    }
    MatrixXd CreateRoutesDelayColumn() {
      MatrixXd delay_column(routes_.size(), 1);
      for (int i = 0; i < routes_.size(); i++)
        delay_column(i, 0) = Link<T>::GetLinksDelay(links_, routes_[i]);
      return delay_column;
    }
    MatrixXd CreateEColumn() {
      MatrixXd e(routes_.size(), 1);
      for (int i = 0; i < routes_.size(); i++)
        e(i, 0) = 1;
      return e;
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
    void BalanceRoutes();
    MatrixXd CreateRoutesJacobiMatrix() {
      int m = routes_.size();
      MatrixXd res(m, m);
      for (int i = 0; i < m; i++) {
        unordered_set <int> links_cur_route;
        for (auto now : routes_[i])
          links_cur_route.insert(now);
        for (int j = 0; j < m; j++) {
          res(i, j) = 0;
          for (auto now : routes_[j])
            if (links_cur_route.count(now))
              res(i, j) += links_[now].DelayDer();
        }
      }
      return res;
    }
    T ConstRouteDelay() {
      T delay_const_route = -1;
      for (const auto& route : routes_) {
        bool fl = true;
        for (const auto& now : route) {
          if (links_[now].power != 0)
            fl = false;
        }
        if (fl)
          delay_const_route = Link<T>::GetLinksDelay(links_, route);
      }
      return delay_const_route;
    }
    MatrixXd CreateT() {
      T delay_const_route = ConstRouteDelay();
      if (delay_const_route == -1) {
        MatrixXd e = CreateEColumn(), jacobi_inv = CreateRoutesJacobiMatrix().inverse();
        return (e.transpose() * jacobi_inv * CreateRoutesDelayColumn()) / (e.transpose() * jacobi_inv * e)(0, 0);
      }
      else {
        MatrixXd res(1, 1);
        res(0, 0) = delay_const_route;
        return res;
      }
    }
    T CalculateExcess(const MatrixXd& routes_flow) {
      T total_routes_flow = 0;
      for (int i = 0; i < routes_flow.rows(); i++)
        total_routes_flow += routes_flow(i, 0);
      return total_routes_flow - demand_;
    }
    void BalanceRoutes2() {
      ClearFlow();
      int m = routes_.size();
      T excess = 0;
      MatrixXd routes_flow(m, 1);
      for (int i = 0; i < m; i++)
        routes_flow(i, 0) = demand_ / m;
      UpdateRoutesFlow(routes_flow);
      while (GetCurrentDelta() > eps_) {
        routes_flow -= CreateRoutesJacobiMatrix().inverse() * (CreateRoutesDelayColumn() - CreateEColumn() * CreateT());
        excess = CalculateExcess(routes_flow);
        //cout << excess << '\n';
        for (int i = 0; i < m; i++)
          routes_flow(i, 0) -= excess / m;
        ClearFlow();
        UpdateRoutesFlow(routes_flow);
        //cout << routes_flow << '\n';
      }
    }
  public:
    OriginDestinationPair(int& origin, int& dest, T& demand, vector <Link <T>>& links, vector <vector <int>>& adjacency_list) :
      origin_(origin), dest_(dest), demand_(demand), links_(links), adjacency_list_(adjacency_list), current_delay_(0) { };
    vector <int> BestRoute() {
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
    bool AddNewRoute(vector <int> new_route) {
      for (auto route : routes_)
        if (new_route == route)
          return false;
      routes_.push_back(new_route);
      for (auto now : routes_[routes_.size() - 1])
        if (routes_.size() == 1) {
          links_[now].flow += demand_;
          links_flow_[now] = demand_;
        }
        else if (!links_flow_.count(now))
          links_flow_[now] = 0;
      return true;
    }
    T GetCurrentDelta() {
      T delay_min = MAX_ROUTE_DELAY, delay_max = 0, delay_route = 0;
      if (routes_.size() > 0) {
        delay_route = Link<T>::GetLinksDelay(links_, routes_[0]);
        delay_min = delay_route;
        delay_max = delay_route;
      }
      for (int i = 1; i < routes_.size(); i++) {
        delay_route = Link<T>::GetLinksDelay(links_, routes_[i]);
        delay_min = min(delay_min, delay_route);
        delay_max = max(delay_max, delay_route);
      }
      //if (abs(delay_max - delay_min) > 1)
      //cout << abs(delay_max - delay_min) << '\n';
      return abs(delay_max - delay_min);
    }
    void NextUpdate() {
      //auto start = chrono::high_resolution_clock::now();
      ClearFlow();
      vector <vector <int>> possible_routes = routes_;
      routes_.clear();
      T delay = MAX_ROUTE_DELAY;
      while (!possible_routes.empty()) {
        int best_route = 0;
        T delay_best = MAX_ROUTE_DELAY, delay_now = 0;
        for (int i = 0; i < possible_routes.size(); i++) {
          delay_now = Link<T>::GetLinksDelay(links_, possible_routes[i]);
          if (delay_now < delay_best) {
            best_route = i;
            delay_best = delay_now;
          }
        }
        if (delay_best < delay) {
          routes_.push_back(possible_routes[best_route]);
          //BalanceRoutes();
          BalanceRoutes2();
          delay = GetDelay();
        }
        swap(possible_routes[best_route], possible_routes[possible_routes.size() - 1]);
        possible_routes.pop_back();
        //auto end = chrono::high_resolution_clock::now();
        //chrono::duration<float> duration = end - start;
        //cout << "Duration of calc_next " << duration.count() << " sec \n";
      }
    }
    int RoutesCount() {
      return routes_.size();
    }
    T GetDemand() {
      return demand_;
    }
    pair <int, int> GetOriginDestination() {
      return { origin_, dest_ };
    }
    vector <T> RoutesDelays() {
      vector <T> routes_delays(routes_.size());
      for (int i = 0; i < routes_.size(); i++)
        routes_delays[i] = Link<T>::GetLinksDelay(links_, routes_[i]);
      return routes_delays;
    }
    vector <vector <int>> RoutesInfo() {
      return routes_;
    }
    void SetDefaultFlow() {
      ClearFlow();
      for (const auto& route : routes_) {
        for (auto now : route) {
          links_[now].flow += demand_ / routes_.size();
          links_flow_[now] += demand_ / routes_.size();
        }
      }
    }
    void SetNewRoutesInfo(const vector <vector <int>>& new_routes, const vector <T>& new_routes_flow) {
      ClearFlow();
      routes_ = new_routes;
      for (int i = 0; i < new_routes_flow.size(); i++) 
        for (auto now : new_routes[i]) {
          links_[now].flow += new_routes_flow[i];
          links_flow_[now] += demand_ / routes_.size();
        }
    }
  };

  template <class T>
  void OriginDestinationPair<T>::BalanceRoutes() {
    ClearFlow();
    tuple <unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, unordered_map <int, int>, MatrixXd > comp = Compression(routes_);
    unordered_map <int, int> mapping_nodes = get<0>(comp), mapping_nodes_inv = get<1>(comp), mapping_edges = get<2>(comp), mapping_edges_inv = get<3>(comp);
    MatrixXd compressed_adjacency_matrix = get<4>(comp);
    int n = compressed_adjacency_matrix.rows(), m = compressed_adjacency_matrix.cols();
    MatrixXd temp_flow(m, 1);
    for (const auto& route : routes_)
      for (const auto& now : route)
        temp_flow(mapping_edges[now], 0) += demand_ / routes_.size();
    UpdateFlow(temp_flow, mapping_edges_inv);
    while (GetCurrentDelta() > eps_) {
      ClearFlow();
      MatrixXd g = CreateG(temp_flow, mapping_edges_inv);
      temp_flow -= g * (MatrixXd::Identity(m, m) - compressed_adjacency_matrix.transpose() *
        (compressed_adjacency_matrix * g * compressed_adjacency_matrix.transpose()).inverse() *
        compressed_adjacency_matrix * g) * CreateDelayColumn(temp_flow, mapping_edges_inv);
      //cout << temp_flow << '\n';
      UpdateFlow(temp_flow, mapping_edges_inv);
    }
  }

  template <typename T>
  class TrafficAssignmentApproach {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  protected:
    T best_ans_;
    vector <Link <T>> links_;
    vector <vector <int>> adjacency_list_;                    // Links for each node
    vector <OriginDestinationPair <T>> origin_destination_pairs_;
    // stores information about all destinations for the single origin
    vector <unordered_map <int, int>> origin_info_;
    int number_of_origin_destination_pairs_, number_of_links_, number_of_nodes_, number_of_zones_;
    string test_name_;
    ofstream time_out, rgap_out, func_out;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    long double time_on_statistics_;
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
          origin_destination_pairs_.push_back(OriginDestinationPair <T>(origin, dest, demand, links_, adjacency_list_));
          origin_info_[origin][dest] = number_of_origin_destination_pairs_ - 1;
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
      origin_info_.resize(number_of_nodes_);
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
    vector <int> RestoreRoute(int origin, int dest, const unordered_map <int, int>& used_link) {
      int now = dest;
      vector <int> new_route;
      while (now != origin) {
        new_route.push_back(used_link.at(now));
        now = links_[used_link.at(now)].init;
      }
      reverse(new_route.begin(), new_route.end());
      return new_route;
    }
    // Dijkstra's algorithm for finding shortest paths to destinatins for the single origin
    void SingleOriginBestRoutes(int origin) {
      priority_queue <pair <T, int>, vector <pair <T, int>>, greater <pair <T, int>>> q;
      int u;
      q.push({ 0, -1 });
      unordered_set <int> processed;
      unordered_map <int, int> used_link;
      T cur_delay;
      int cnt_not_processed = origin_info_[origin].size();
      vector <int> destinations;
      while (cnt_not_processed > 0) {
        if (q.top().second != -1)
          u = links_[q.top().second].term;
        else
          u = origin;
        if (processed.count(u)) {
          q.pop();
          continue;
        }
        cur_delay = q.top().first;
        processed.insert(u);
        used_link[u] = q.top().second;
        q.pop();
        for (auto now : adjacency_list_[u])
          if (!processed.count(links_[now].term))
            q.push({ cur_delay + links_[now].Delay(), now });
        if (origin_info_[origin].count(u)) {
          cnt_not_processed--;
          destinations.push_back(u);
        }
      }
      for (auto dest : destinations)
        // adds new route for od-pair with destination dest
        origin_destination_pairs_[origin_info_[origin][dest]].AddNewRoute(RestoreRoute(origin, dest, used_link));
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
      //cout << RelativeGap() << '\n';
      rgap_out << setprecision(30) << RelativeGap() << '\n';
      //cout << setprecision(30) << ObjectiveFunction() << '\n';
      func_out << setprecision(30) << ObjectiveFunction() << '\n';
      //cout << (chrono::duration_cast<chrono::milliseconds>(std::chrono::high_resolution_clock::now() - statistics_start)).count() << '\n';
      time_on_statistics_ += 0.001 * (chrono::duration_cast<chrono::seconds>(std::chrono::high_resolution_clock::now() - statistics_start)).count();
      time_out << (chrono::duration_cast<chrono::seconds>(std::chrono::high_resolution_clock::now() - start_)).count() - time_on_statistics_ << '\n';
      time_out.close();
      func_out.close();
      rgap_out.close();
    }
    virtual void SolveFlow() {};
  };

  template <typename T>                    // value type
  class RouteBasedApproach : public TrafficAssignmentApproach <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  private:
    queue <int> ready_to_implement_, tasks_;
    string test_name_;
    T alpha_ = 1e-30;
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
  public:
    RouteBasedApproach(string test_name) : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(test_name) {
    }
    void SolveFlow() override {
      int count = 0;
      while (count++ < 20) {
        //NewRouteFound();
        for (int origin = 0; origin < TrafficAssignmentApproach<T>::number_of_nodes_; origin++)
          TrafficAssignmentApproach<T>::SingleOriginBestRoutes(origin);
        auto start = chrono::high_resolution_clock::now();
        while (GetDelta() > alpha_) {
          for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++) {
            TrafficAssignmentApproach<T>::origin_destination_pairs_[t].NextUpdate();
            //cout << TrafficAssignmentApproach<T>::origin_destination_pairs_[t].GetCurrentDelta() << '\n';
          }
        }
        TrafficAssignmentApproach<T>::GetStatistics();
      }
    }
  };

  template <typename T>                    // value type
  class NewRouteBasedApproach : public TrafficAssignmentApproach <T> {
    using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  private:
    queue <int> ready_to_implement_, tasks_;
    string test_name_;
    T alpha_ = 1e-10;
    mutex mtx_, mtx_tasks_;
    T GetDelta() {
      T delta = 0;
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++)
        delta = max(delta, TrafficAssignmentApproach<T>::origin_destination_pairs_[t].GetCurrentDelta());
      //cout << delta << '\n';
      return delta;
    }
    int TotalNumberOfRoutes() {
      int total_route_count = 0;
      for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++)
        total_route_count += TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesCount();
      return total_route_count;
    }
    MatrixXd EMatrix() {
      int total_route_count = TotalNumberOfRoutes();
      MatrixXd e(TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_, total_route_count);
      for (int i = 0, temp = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++) {
        for (int j = 0; j < TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesCount(); j++)
          e(i, temp + j) = 1;
        temp += TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesCount();
      }
      return e;
    }
    MatrixXd AllRoutesDelaysColumn() {
      int total_route_count = TotalNumberOfRoutes();
      MatrixXd all_routes_delays(total_route_count, 1);
      for (int i = 0, temp = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++) {
        vector <T> od_routes_delays = TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesDelays();
        for (int j = 0; j < od_routes_delays.size(); j++) 
          all_routes_delays(temp + j, 0) = od_routes_delays[j];
        temp += od_routes_delays.size();
      }
      return all_routes_delays;
    }
    MatrixXd AllRoutesJacobiMatrix() {
      int total_route_count = TotalNumberOfRoutes();
      vector <vector <int>> all_routes(total_route_count);
      for (int i = 0, temp = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++) {
        vector <vector <int>> od_routes = TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesInfo();
        for (int j = 0; j < od_routes.size(); j++)
          all_routes[temp + j] = od_routes[j];
        temp += od_routes.size();
      }
      MatrixXd res(total_route_count, total_route_count);
      for (int i = 0; i < total_route_count; i++) {
        unordered_set <int> links_cur_route;
        for (auto now : all_routes[i])
          links_cur_route.insert(now);
        for (int j = 0; j < total_route_count; j++) {
          res(i, j) = 0;
          for (auto now : all_routes[j])
            if (links_cur_route.count(now))
              res(i, j) += TrafficAssignmentApproach<T>::links_[now].DelayDer();
        }
      }
      return res;
    }
    MatrixXd TMatrix(const MatrixXd& e, const MatrixXd& rev_jacobi) {
      return (e * rev_jacobi * e.transpose()).inverse() * (e * rev_jacobi * AllRoutesDelaysColumn());
    }
    MatrixXd NextUpdate(const MatrixXd& all_routes_flows) {
      int total_route_count = TotalNumberOfRoutes();
      MatrixXd e = EMatrix(), rev_jacobi = AllRoutesJacobiMatrix().inverse();
      return all_routes_flows - rev_jacobi * (AllRoutesDelaysColumn() - e.transpose() * TMatrix(e, rev_jacobi));
    }
    vector <T> CorrectNextRouteFlowInfo(const int origin_dest_pair, const vector <T>& next_flow) {
      int cur_routes_cnt = TrafficAssignmentApproach<T>::origin_destination_pairs_[origin_dest_pair].RoutesCount();
      vector <vector <int>> next_pos_routes, prev_routes = TrafficAssignmentApproach<T>::origin_destination_pairs_[origin_dest_pair].RoutesInfo();
      vector <T> next_pos_flow;
      T temp = 0;
      for (int i = 0; i < cur_routes_cnt; i++) 
        if (next_flow[i] > 0) {
          next_pos_routes.push_back(prev_routes[i]);
          next_pos_flow.push_back(next_flow[i]);
          temp += next_flow[i];
        }
      for (int i = 0; i < next_pos_routes.size(); i++)
        next_pos_flow[i] *= TrafficAssignmentApproach<T>::origin_destination_pairs_[origin_dest_pair].GetDemand() / temp;
      TrafficAssignmentApproach<T>::origin_destination_pairs_[origin_dest_pair].SetNewRoutesInfo(next_pos_routes, next_pos_flow);
      return next_pos_flow;
    }
    MatrixXd CorrectAllRouteFlowInfo(const MatrixXd& all_routes_flow) {
      int cur_route_cnt = TotalNumberOfRoutes(), next_total_cnt = 0;
      for (int i = 0; i < cur_route_cnt; i++)
        if (all_routes_flow(i, 0) > 0)
          next_total_cnt++;
      MatrixXd next_all_routes_flow(next_total_cnt, 1);
      int k = 0;
      for (int i = 0, temp = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++) {
        int cur_routes_cnt = TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesCount();
        vector <T> next_flow(cur_routes_cnt);
        for (int j = 0; j < cur_routes_cnt; j++)
          next_flow[j] = all_routes_flow(temp + j, 0);
        temp += cur_routes_cnt;
        vector <T> next_pos_routes = CorrectNextRouteFlowInfo(i, next_flow);
        for (int j = 0; j < next_pos_routes.size(); j++)
          next_all_routes_flow(k + j, 0) = next_pos_routes[j];
        k += next_pos_routes.size();
      }
      return next_all_routes_flow;
    }
  public:
    NewRouteBasedApproach(string test_name) : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(test_name) {
    }
    void SolveFlow() override {
      int count = 0;
      while (count++ < 20) {
        for (int origin = 0; origin < TrafficAssignmentApproach<T>::number_of_nodes_; origin++)
          TrafficAssignmentApproach<T>::SingleOriginBestRoutes(origin);
        int total_route_count = TotalNumberOfRoutes();
        MatrixXd all_routes_flow(total_route_count, 1);
        for (int i = 0, temp = 0; i < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; i++) {
          TrafficAssignmentApproach<T>::origin_destination_pairs_[i].SetDefaultFlow();
          int cur_routes_cnt = TrafficAssignmentApproach<T>::origin_destination_pairs_[i].RoutesCount();
          T cur_demand = TrafficAssignmentApproach<T>::origin_destination_pairs_[i].GetDemand();
          for (int j = 0; j < cur_routes_cnt; j++)
            all_routes_flow(temp + j, 0) = cur_demand / cur_routes_cnt;
          temp += cur_routes_cnt;
        }
        auto start = chrono::high_resolution_clock::now();
        while (GetDelta() > alpha_)
          all_routes_flow = CorrectAllRouteFlowInfo(NextUpdate(all_routes_flow));
        TrafficAssignmentApproach<T>::GetStatistics();
      }
    }
  };

  template <typename T>
  class SolutionCheck : public TrafficAssignmentApproach <T> {
  private: 
    map <pair <int, int>, int> origin_destination_number_, link_number_;
    vector <vector <vector <int>>> od_routes_info_;
    vector <vector <T>> od_routes_flow_;
    //vector <T> solution_od_routes_flow, solution_od_routes_info;
    int GetOriginDestinationPair(ifstream& in) {
      int origin, dest;
      in >> origin >> dest;
      origin--; dest--;
      //cout << origin << ' ' << dest << ' ' << origin_destination_number_[{origin, dest}] << ' ';
      return origin_destination_number_[{origin, dest}];
    }
    vector <T> GetOriginDestinationFLows(ifstream& in) {
      char c = '#';
      vector <T> flows;
      T temp = 0, k = 1;
      while (c != '{')
        in >> c;
      //cout << '[';
      while (c != '}') {
        if (in.peek() == ' ') {
          in.ignore();
          continue;
        }
        if (isdigit(in.peek())) {
          in >> temp;
          //cout << temp << ' ';
          flows.push_back(temp);
        }
        else 
          in >> c;
      }
      //cout << "]\n";
      return flows;
    }
    vector <int> ConvertNodesRouteIntoLinks(const vector <int>& nodes) {
      vector <int> res;
      for (int i = 1; i < nodes.size(); i++) 
        res.push_back(link_number_[{nodes[i - 1], nodes[i]}]);
      return res;
    }
    vector <int> GetRoute(ifstream& in) {
      vector <int> nodes, route;
      int temp;
      char c; in >> c;
      while (c != '[')
        in >> c;
      while (c != ']') {
        if (in.peek() == ' ') {
          in.ignore();
          continue;
        }
        if (isdigit(in.peek())) {
          in >> temp;
          temp--;
          //cout << temp << ' ';
          nodes.push_back(temp);
        }
        else
          in >> c;
      }
      return ConvertNodesRouteIntoLinks(nodes);
    }
    vector <vector <int>> GetOriginDestinationRoutes(ifstream& in, int number_of_routes) {
      vector <vector <int>> routes;
      for (int i = 0; i < number_of_routes; i++)
        routes.push_back(GetRoute(in));
      return routes;
    }
    void GetSolutionInfo(string test_name) {
      string line;
      ifstream in(test_name + "_flows.txt");
      getline(in, line);
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++) {
        cout << setprecision(15);
        int i = GetOriginDestinationPair(in);
        //cout << "OD_pair: i = " << i << "\n     ";
        od_routes_flow_[i] = GetOriginDestinationFLows(in);
        od_routes_info_[i] = GetOriginDestinationRoutes(in, od_routes_flow_[i].size());
        /*for (int j = 0; j < od_routes_flow_[i].size(); j++) {
          cout << "ROUTE: ";
          for (auto now : od_routes_info_[i][j])
            cout << now << ' ';
          cout << "ROUTE FLOW: " << od_routes_flow_[i][j] << "\n     ";
        }
        cout << '\n';*/
        T temp; in >> temp;
      }
    }
    void ContributeFlowValues() {
      for (int i = 0; i < TrafficAssignmentApproach<T>::number_of_links_; i++)
        TrafficAssignmentApproach<T>::links_[i].flow = 0;
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++) {
        for (int i = 0; i < od_routes_info_[t].size(); i++)
          for (auto now : od_routes_info_[t][i])
            TrafficAssignmentApproach<T>::links_[now].flow += od_routes_flow_[t][i];
      }
    }
    T IncludeExcessFlow() {
      T total_excess = 0;
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++) {
        T excess = TrafficAssignmentApproach<T>::origin_destination_pairs_[t].GetDemand();
        for (int i = 0; i < od_routes_flow_[t].size(); i++)
          excess -= od_routes_flow_[t][i];
        for (int i = 0; i < od_routes_info_[t].size(); i++)
          od_routes_flow_[t][i] += excess / od_routes_info_[t].size();
        total_excess += excess;
      }
      ContributeFlowValues();
      return total_excess;
    }
  public:
    SolutionCheck(string test_name) : TrafficAssignmentApproach<T>::TrafficAssignmentApproach(test_name) {
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_; t++)
        origin_destination_number_[TrafficAssignmentApproach<T>::origin_destination_pairs_[t].GetOriginDestination()] = t;
      for (int t = 0; t < TrafficAssignmentApproach<T>::number_of_links_; t++)
        link_number_[{TrafficAssignmentApproach<T>::links_[t].init, TrafficAssignmentApproach<T>::links_[t].term}] = t;
      od_routes_info_.resize(TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_);
      od_routes_flow_.resize(TrafficAssignmentApproach<T>::number_of_origin_destination_pairs_);
      GetSolutionInfo(test_name);
      ContributeFlowValues();
      cout << setprecision(30) << TrafficAssignmentApproach<T>::ObjectiveFunction() << '\n';
      cout << IncludeExcessFlow() << '\n';
      cout << setprecision(30) << TrafficAssignmentApproach<T>::ObjectiveFunction() << '\n';

    }
    void SolveFlow() {
      cout << 1;
    }
    void Check() {
      T total_flow_loss = 0;

    }
  };
  /*
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
  */
}
