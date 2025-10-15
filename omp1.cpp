// omp1.cpp — DBSCAN Paralelo (Versión 1: ingenua)
// Paraleliza SOLO el cálculo de vecindarios (O(n^2));
// la expansión de clusters se mantiene en serie para evitar condiciones de carrera.
// Compilación:
//   /opt/homebrew/bin/g++-14 -O3 -std=c++17 -fopenmp omp1.cpp -o omp1
// Ejemplo:
//   ./omp1 --threads 8 --n 50000 --d 2 --eps 1.5 --minpts 8 --seed 42

#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <iostream>
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

// -------------------------- DBSCAN OMP1 (ingenua) ------------------------------
Result dbscan_omp1(const vector<Point>& P, const Params& param, int threads){
  omp_set_num_threads(threads);
  const int n = (int)P.size();
  const double eps2 = param.eps * param.eps;

  vector<vector<int>> neigh(n);
  vector<uint8_t> is_core(n, 0);

  // 1) Vecindarios en paralelo (cada hilo trabaja sobre su i)
  #pragma omp parallel for schedule(static)
  for(int i=0;i<n;i++){
    vector<int> local; local.reserve(64);
    for(int j=0;j<n;j++){
      if (dist2(P[i],P[j]) <= eps2) local.push_back(j);
    }
    neigh[i] = std::move(local);
    if ((int)neigh[i].size() >= param.minPts) is_core[i]=1;
  }

  // 2) Expansión en serie (idéntica a la serial) para mantener correctitud
  Result R; R.label.assign(n, -1); R.is_core = std::move(is_core);
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

// ------------------------------- Dataset sintético ------------------------------
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
  int threads = 8;           // por defecto
  int n = 30000, d = 2;      // tamaño y dimensión del dataset
  double eps = 1.5;          // radio
  int minPts = 8;            // densidad mínima
  unsigned seed = 42;        // reproducibilidad

  for (int i=1; i<argc; ++i){
    if (string(argv[i])=="--threads" && i+1<argc) threads = stoi(argv[++i]);
    else if (string(argv[i])=="--n" && i+1<argc) n = stoi(argv[++i]);
    else if (string(argv[i])=="--d" && i+1<argc) d = stoi(argv[++i]);
    else if (string(argv[i])=="--eps" && i+1<argc) eps = stod(argv[++i]);
    else if (string(argv[i])=="--minpts" && i+1<argc) minPts = stoi(argv[++i]);
    else if (string(argv[i])=="--seed" && i+1<argc) seed = (unsigned)stoul(argv[++i]);
  }

  auto P = make_synthetic(n, d, seed);
  Params param{eps, minPts};

  double t0 = omp_get_wtime();
  Result R = dbscan_omp1(P, param, threads);
  double t1 = omp_get_wtime();

  cout << "impl=omp1 threads=" << threads
       << " n=" << n << " d=" << d
       << " eps=" << eps << " minPts=" << minPts << "\n";
  cout << "time_s=" << (t1 - t0) << "\n";
  cout << "clusters=" << R.n_clusters
       << " core=" << R.n_core
       << " border=" << R.n_border
       << " noise=" << R.n_noise << "\n";

  return 0;
}