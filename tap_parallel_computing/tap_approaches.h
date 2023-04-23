namespace traffic_assignment {
template <class T>                    // value type
struct Link {
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
  static T GetLinksDelay(vector <Link <T>>& links, const vector <int>& links_list) {
    T ans = 0;
    for (auto now : links_list)
      ans += links[now].Delay();
    return ans;
  }
};

template <class T>
class OriginDestinationPair {
private:
  int origin_, dest_;
  T demand_, current_delay_, next_delay_, next_delta_;
  vector <Link <T>> &links_;
  const vector <vector <int>> &adjacency_list_;
  vector <vector <int>> current_routes_, next_routes_;
  unordered_map <int, T> current_links_flow_, next_links_flow_;
  unordered_map <int, T> fixed_links_flow_;
  const T eps_ = 1e-6;
  matrix <T> CreateG(matrix <T>& temp_flow, map <int, int>& M_links_inv) {
    int n = temp_flow.size().first;
    vector <vector <T>> Temp(n, vector <T>(n, 0));  
    for (int i = 0; i < n; i++)
      Temp[i][i] = 1 / links_[M_links_inv[i]].DelayDer(fixed_links_flow_[M_links_inv[i]] + temp_flow[i][0]);
    return matrix <T>(n, n, Temp);
  }
  matrix <T> CreateDelayColumn(matrix <T>& temp_flow, map <int, int>& M_links_inv) {
    matrix <T> g = temp_flow;
    for (int i = 0; i < temp_flow.size().first; i++)
      g[i][0] = links_[M_links_inv[i]].Delay(fixed_links_flow_[M_links_inv[i]] + temp_flow[i][0]);
    return g;
  }
  tuple < map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <T> > Compression(const vector <vector <int>>& routes) {
    int count = 0;
    set <int> actually_used_links, actually_used_nodes;
    map <int, int> mapping_nodes, mapping_nodes_inv, mapping_edges, mapping_edges_inv;
    matrix <T> A;
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
    A = matrix <T>(n, m);
    count = 0;
    for (auto now : actually_used_links) {
      mapping_edges[now] = count;
      mapping_edges_inv[count] = now;
      if (links_[now].init != origin_)
        A[mapping_nodes[links_[now].init]][count] = 1;
      if (links_[now].term != origin_)
        A[mapping_nodes[links_[now].term]][count] = -1;
      count++;
    }
    return { mapping_nodes, mapping_nodes_inv, mapping_edges, mapping_edges_inv, A };
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
    tuple <map <int, int>, map <int, int>, map <int, int>, map <int, int>, matrix <T> > comp = Compression(next_routes_);
    map <int, int> mapping_nodes = get<0>(comp), mapping_nodes_inv = get<1>(comp), mapping_edges = get<2>(comp), mapping_edges_inv = get<3>(comp);
    matrix <T> A = get<4>(comp), G, A_trans, E, g, t;
    int n = A.size().first, m = A.size().second;
    matrix <T> temp_flow_now(m, 1), temp_flow_prev(m, 1);
    for (int i = 0; i < next_routes_.size(); i++)
      for (auto& now : next_routes_[i])
        temp_flow_now[mapping_edges[now]][0] += demand_ / next_routes_.size();
    next_delta_ = eps_ + 1;
    next_links_flow_.clear();
    for (int i = 0; i < m; i++)
      next_links_flow_[mapping_edges_inv[i]] = temp_flow_now[i][0];
    UpdateNextData();
    while (next_delta_ > eps_) {
      temp_flow_prev = temp_flow_now;
      G = CreateG(temp_flow_now, mapping_edges_inv);
      A_trans = A.trans();
      E = matrix<T>::create_e(m, m);
      g = CreateDelayColumn(temp_flow_now, mapping_edges_inv);
      t = (A * G * A_trans).inv();
      t = A_trans * t * A * G;
      t = E - t;
      t = G * t * g;
      temp_flow_now = temp_flow_now - t;
      next_links_flow_.clear();
      for (int i = 0; i < m; i++)
        next_links_flow_[mapping_edges_inv[i]] = temp_flow_now[i][0];
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
  OriginDestinationPair(int &origin, int &dest, T &demand, vector <Link <T>> &links, vector <vector <int>> &adjacency_list) :
    origin_(origin), dest_(dest), demand_(demand), links_(links), adjacency_list_(adjacency_list), current_delay_(0) { };
  bool FindNewRoute() {
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
    T delay_cur_route = q.top().first;
    if (delay_cur_route < current_delay_ || current_routes_.size() == 0) {
      current_delay_ = delay_cur_route;
      int now = dest_;
      vector <int> new_route;
      while (now != origin_) {
        new_route.push_back(used_link[now]);
        now = links_[used_link[now]].init;
      }
      reverse(new_route.begin(), new_route.end());
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
  bool ImplementNextUpdate() {
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
      return true;
    }
    else
      return false;
  }
  int RoutesCount() {
    return current_routes_.size();
  }
};

template <typename T>                    // value type
class RouteBasedApproach {
private:
  vector <Link <T>> links_;
  vector <vector <int>> adjacency_list_;              // Links for each node
  vector <OriginDestinationPair <T>> origin_destination_pairs_;
  int number_of_origin_destination_pairs_, number_of_links_, number_of_nodes_, number_of_zones_;
  queue <int> ready_to_implement_, tasks_;
  int tasks_count;
  string test;
  T alpha = 1e-4;
  mutex mtx_, mtx_tasks_;
public:
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
    ifstream in(test + "_net.txt");
    while (getline(in, line))
      GetLink(line);
    in = ifstream(test + "_trips.txt");
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
  RouteBasedApproach(string tn) : number_of_links_(0), number_of_origin_destination_pairs_(0), test(tn) {
    GetData();
    adjacency_list_.resize(number_of_nodes_);
    for (int i = 0; i < number_of_links_; i++)
      adjacency_list_[links_[i].init].push_back(i);
  }
  T ObjectiveFunction() {
    T ans = 0;
    for (int i = 0; i < number_of_links_; i++)
      ans += links_[i].DelayInteg();
    return ans;
  }
  T GetBestAnswer() {
    T ans = 0, flow_on_link;
    string line;
    ifstream in(test + "_flow.txt");
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
  void ShowStatistics() {
    cout << "best ans: " << GetBestAnswer() << '\n';
    cout << "obj function: " << ObjectiveFunction() << '\n';
    //for (int i = 0; i < number_of_links_; i++)
    //  cout << "from: " << links[i].init + 1 << " to " << links[i].term + 1 << " volume " << links[i].flow << " cost " <<
    //  links[i].Delay() << '\n';
  }
  T GetDelta() {
    T delta = 0;
    for (int t = 0; t < number_of_origin_destination_pairs_; t++)
      delta = max(delta, origin_destination_pairs_[t].GetCurrentDelta());
    //cout << delta << '\n';
    return delta;
  }
  bool NewRouteFound() {
    bool fl = false;
    for (auto& now : origin_destination_pairs_)
      fl = (now.FindNewRoute() || fl);
    return fl;
  }
  void ThreadCalculation() {
    while (!tasks_.empty()) {
      int j = -1;
      mtx_tasks_.lock();
      if (!tasks_.empty()) {
        j = tasks_.front();
        tasks_.pop();
      }
      mtx_tasks_.unlock();
      if (j == -1)
        continue;
      origin_destination_pairs_[j].CalculateNextUpdate(mtx_);
      mtx_.lock();
      origin_destination_pairs_[j].ImplementNextUpdate();
      mtx_.unlock();
    }
  }
  void SolveFlow() {
    int num_threads = thread::hardware_concurrency() - 1;
    //num_threads = 1;
    vector <thread> threads_list(num_threads);
    vector <int> tasks_order(number_of_origin_destination_pairs_);
    vector <int> now_in_process(num_threads, 0);
    for (int i = 0; i < number_of_origin_destination_pairs_; i++)
      tasks_order[i] = i;
    int count = 0;
    while (NewRouteFound()) {
      auto start = chrono::high_resolution_clock::now();
      while (GetDelta() > alpha) {
        random_shuffle(tasks_order.begin(), tasks_order.end());
        for (int i = 0; i < number_of_origin_destination_pairs_; i++)
          if (origin_destination_pairs_[tasks_order[i]].RoutesCount() > 1)
            tasks_.push(tasks_order[i]);
        for (int i = 0; i < num_threads; i++)
          threads_list[i] = thread(&RouteBasedApproach<T>::ThreadCalculation, this);
        for_each(threads_list.begin(), threads_list.end(), mem_fn(&thread::join));
      }
      auto end = chrono::high_resolution_clock::now();
      chrono::duration<float> duration = end - start;
      cout << "Duration of balance " << duration.count() << " sec \n";
      //ShowStatistics();
    }
    ShowStatistics();
  }
};
}
