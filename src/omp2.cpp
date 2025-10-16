// omp2.cpp
// Compilar:
//   /opt/homebrew/bin/g++-15 -O3 -std=c++17 -fopenmp src/omp2.cpp -o src/omp2
// Ejecutar:
//   ./src/omp2 --threads 8 --in input/20000_data.csv --eps 0.03 --minpts 10 --out output/20000_results.csv

#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_map>
#include <omp.h>

using namespace std;

// ------------------------------- Tipos básicos ---------------------------------
struct Point { vector<double> x; };
struct Params { double eps; int minPts; };
struct Result {
  vector<int> label;        // -1 (interno), -2 noise, >=0 cluster
  vector<uint8_t> is_core;  // 1 si core, 0 no-core
  int n_clusters = 0;
  int n_core = 0, n_border = 0, n_noise = 0;
};

static inline double dist2(const Point& a, const Point& b){
  double s=0; for (size_t k=0;k<a.x.size();++k){ double d=a.x[k]-b.x[k]; s+=d*d; } return s;
}

static inline void compute_counts(Result& R){
  int core=0, border=0, noise=0;
  for (size_t i=0;i<R.label.size();++i){
    if (R.label[i] == -2) noise++;
    else if (R.is_core[i]) core++;
    else border++;
  }
  R.n_core=core; R.n_border=border; R.n_noise=noise;
}

// ---------------------------------- DSU ----------------------------------------
struct DSU {
  vector<int> p, r;
  DSU(int n=0){ init(n); }
  void init(int n){ p.resize(n); r.assign(n,0); for(int i=0;i<n;i++) p[i]=i; }
  int find(int x){ while(p[x]!=x){ p[x]=p[p[x]]; x=p[x]; } return x; }
  void unite(int a, int b){
    a=find(a); b=find(b); if(a==b) return;
    if (r[a]<r[b]) p[a]=b; else if (r[b]<r[a]) p[b]=a; else { p[b]=a; r[a]++; }
  }
};

// ------------------------- Grilla espacial (hash 2D) ---------------------------
struct CellKey { long long x, y; };
struct CellKeyHash{
  size_t operator()(const CellKey& k) const noexcept {
    return std::hash<long long>{}((k.x<<32) ^ (k.y & 0xffffffff));
  }
};
struct CellKeyEq{ bool operator()(const CellKey& a, const CellKey& b) const noexcept { return a.x==b.x && a.y==b.y; } };

// ------------------------------- DBSCAN OMP2 -----------------------------------
Result dbscan_omp2(const vector<Point>& P, const Params& param, int threads){
  omp_set_num_threads(threads);
  const int n = (int)P.size();
  const double eps = param.eps; const double eps2 = eps*eps;

  const int T = omp_get_max_threads();

  auto cell_of = [&](const Point& p){
    double x0 = p.x.size()>=1 ? p.x[0] : 0.0;
    double x1 = p.x.size()>=2 ? p.x[1] : 0.0;
    long long cx = (long long)floor(x0/eps);
    long long cy = (long long)floor(x1/eps);
    return CellKey{cx,cy};
  };

  // 1) Construcción de bins con buffers por hilo y fusión
  vector<vector<pair<CellKey,int>>> local_bins; local_bins.resize(T);
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_bins[tid].reserve(n/threads + 1);
    #pragma omp for schedule(static)
    for(int i=0;i<n;i++){
      local_bins[tid].push_back({cell_of(P[i]), i});
    }
  }
  unordered_map<CellKey, vector<int>, CellKeyHash, CellKeyEq> bins;
  bins.reserve(static_cast<size_t>(n * 1.3));
  for(auto& v: local_bins){ for(auto& kv : v) bins[kv.first].push_back(kv.second); }

  auto neighbors_in_bins = [&](int i, vector<int>& out){
    CellKey c = cell_of(P[i]);
    for(int dx=-1; dx<=1; ++dx){
      for(int dy=-1; dy<=1; ++dy){
        CellKey nb{c.x+dx, c.y+dy};
        auto it = bins.find(nb);
        if (it==bins.end()) continue;
        const auto& vec = it->second;
        for(int j: vec){ if (dist2(P[i],P[j]) <= eps2) out.push_back(j); }
      }
    }
  };

  auto count_neighbors_ge_min = [&](int i, int minReq){
    int cnt = 0;
    CellKey c = cell_of(P[i]);
    for(int dx=-1; dx<=1; ++dx){
      for(int dy=-1; dy<=1; ++dy){
        CellKey nb{c.x+dx, c.y+dy};
        auto it = bins.find(nb);
        if (it==bins.end()) continue;
        const auto& vec = it->second;
        for(int j: vec){ if (dist2(P[i],P[j]) <= eps2){ if(++cnt>=minReq) return true; } }
      }
    }
    return false;
  };

  // 2) CORE en paralelo y 3) recolección de aristas CORE-CORE
  vector<uint8_t> is_core(n,0);
  vector<vector<pair<int,int>>> edges_local(T);
  for(auto &v : edges_local) v.reserve(1024);

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    vector<int> buf; buf.reserve(64);

    #pragma omp for schedule(static)
    for(int i=0;i<n;i++){
      if (count_neighbors_ge_min(i, param.minPts)) is_core[i]=1;
    }

    #pragma omp for schedule(static)
    for(int i=0;i<n;i++){
      if (!is_core[i]) continue;
      buf.clear(); neighbors_in_bins(i, buf);
      for(int j : buf){ if (j<=i) continue; if (!is_core[j]) continue; 
        edges_local[tid].emplace_back(i,j);
      }
    }
  }

  DSU dsu(n);
  vector<pair<int,int>> edges; size_t total_edges=0; 
  for(auto &v: edges_local) total_edges += v.size();
  edges.reserve(total_edges);
  for(auto &v: edges_local){ edges.insert(edges.end(), v.begin(), v.end()); }
  for(const auto &e: edges){ dsu.unite(e.first, e.second); }

  // 4) Etiquetas para CORE por raíz del DSU
  Result R; R.is_core = is_core; R.label.assign(n, -2); // noise por defecto
  unordered_map<int,int> root2cid; root2cid.reserve((size_t)(n*0.5));
  int cid=0;
  for(int i=0;i<n;i++) if (R.is_core[i]){
    int r = dsu.find(i);
    auto it = root2cid.find(r);
    if (it==root2cid.end()) root2cid[r]=cid++;
    R.label[i] = root2cid[r];
  }

  // 5) BORDER: hereda la etiqueta de cualquier core vecino
  #pragma omp parallel for schedule(static)
  for(int i=0;i<n;i++){
    if (R.is_core[i]) continue;
    vector<int> buf; buf.reserve(64);
    neighbors_in_bins(i, buf);
    for(int j: buf){ if (R.is_core[j]) { R.label[i] = R.label[j]; break; } }
  }

  R.n_clusters = cid; compute_counts(R); return R;
}

// ------------------------------- I/O helpers -----------------------------------
static vector<Point> read_csv_points(const string& filename){
  ifstream in(filename);
  vector<Point> P;
  string line;
  while (getline(in, line)){
    stringstream ss(line);
    Point p; double val; char c;
    while (ss >> val){
      p.x.push_back(val);
      ss >> c;
    }
    if(!p.x.empty()) P.push_back(p);
  }
  return P;
}

static void write_csv_results(const string& filename, const vector<Point>& P, const Result& R){
  ofstream out(filename);
  for(size_t i=0;i<P.size();++i){
    int lbl = (R.label[i] == -2 ? 0 : 1);
    out << P[i].x[0] << "," << P[i].x[1] << "," << lbl << "\n";
  }
}

static vector<Point> make_synthetic(int n, int d, unsigned seed){
  mt19937 rng(seed);
  normal_distribution<double> g1(0.0, 1.0), g2(6.0, 1.0);
  vector<Point> P(n); for(auto& p:P) p.x.resize(d);
  for(int i=0;i<n;i++){
    bool c = (i % 2 == 0);
    for(int k=0;k<d;k++) P[i].x[k] = c ? g1(rng) : g2(rng);
  }
  return P;
}

// ------------------------------------ main -------------------------------------
int main(int argc, char** argv){
  int threads = 8;
  int n = 30000, d = 2;
  double eps = 1.5;
  int minPts = 8;
  unsigned seed = 42;
  string in_file = "", out_file = "";

  for (int i=1; i<argc; ++i){
    string a = argv[i];
    if (a=="--threads" && i+1<argc) threads = stoi(argv[++i]);
    else if (a=="--n" && i+1<argc) n = stoi(argv[++i]);
    else if (a=="--d" && i+1<argc) d = stoi(argv[++i]);
    else if (a=="--eps" && i+1<argc) eps = stod(argv[++i]);
    else if (a=="--minpts" && i+1<argc) minPts = stoi(argv[++i]);
    else if (a=="--seed" && i+1<argc) seed = (unsigned)stoul(argv[++i]);
    else if (a=="--in" && i+1<argc) in_file = argv[++i];
    else if (a=="--out" && i+1<argc) out_file = argv[++i];
  }

  vector<Point> P = in_file.empty() ? make_synthetic(n,d,seed) : read_csv_points(in_file);
  Params param{eps, minPts};

  double t0 = omp_get_wtime();
  Result R = dbscan_omp2(P, param, threads);
  double t1 = omp_get_wtime();

  cout << "impl=omp2 threads=" << threads
        << " n=" << n << " d=" << d
        << " eps=" << eps << " minPts=" << minPts << "\n";
  cout << "time_s=" << (t1 - t0) << "\n";
  cout << "clusters=" << R.n_clusters
        << " core=" << R.n_core
        << " border=" << R.n_border
        << " noise=" << R.n_noise << "\n";

  if(!out_file.empty()) write_csv_results(out_file, P, R);
  return 0;
}