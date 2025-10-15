// serial.cpp — Implementación serial de DBSCAN (actualizado con soporte CSV)
// Compilar:
//   /opt/homebrew/bin/g++-15 -O3 -std=c++17 -fopenmp serial.cpp -o serial
// Ejemplo:
//   ./serial --in 4000_data.csv --eps 0.03 --minpts 10 --out 4000_results.csv

#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
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

// ------------------------------- DBSCAN Serial ---------------------------------
Result dbscan_serial(const vector<Point>& P, const Params& param){
  const int n = (int)P.size();
  const double eps2 = param.eps * param.eps;
  Result R; R.label.assign(n, -1); R.is_core.assign(n, 0);

  vector<vector<int>> neigh(n);
  for(int i=0;i<n;i++){
    neigh[i].reserve(64);
    for(int j=0;j<n;j++){
      if (dist2(P[i],P[j]) <= eps2) neigh[i].push_back(j);
    }
    if ((int)neigh[i].size() >= param.minPts) R.is_core[i]=1;
  }

  int cid=0; queue<int> q;
  for(int i=0;i<n;i++){
    if (R.label[i] != -1) continue;
    if (!R.is_core[i]) { R.label[i] = -2; continue; }

    R.label[i]=cid; q.push(i);
    while(!q.empty()){
      int u=q.front(); q.pop();
      if (!R.is_core[u]) continue;
      for(int v : neigh[u]){
        if (R.label[v] == -2) R.label[v] = cid;
        if (R.label[v] == -1){
          R.label[v] = cid;
          if (R.is_core[v]) q.push(v);
        }
      }
    }
    cid++;
  }

  R.n_clusters = cid;
  compute_counts(R);
  return R;
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
  int n = 30000, d = 2;
  double eps = 1.5;
  int minPts = 8;
  unsigned seed = 42;
  string in_file = "", out_file = "";

  for (int i=1; i<argc; ++i){
    string a = argv[i];
    if (a=="--n" && i+1<argc) n = stoi(argv[++i]);
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
  Result R = dbscan_serial(P, param);
  double t1 = omp_get_wtime();

  cout << "impl=serial n=" << n << " d=" << d
       << " eps=" << eps << " minPts=" << minPts << "\n";
  cout << "time_s=" << (t1 - t0) << "\n";
  cout << "clusters=" << R.n_clusters
       << " core=" << R.n_core
       << " border=" << R.n_border
       << " noise=" << R.n_noise << "\n";

  if(!out_file.empty()) write_csv_results(out_file, P, R);
  return 0;
}
