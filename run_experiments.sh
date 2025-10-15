#!/usr/bin/env bash
# -----------------------------------------------------------
# Ejecuta todos los experimentos requeridos para el proyecto
# DBSCAN Serial / OMP1 / OMP2
# -----------------------------------------------------------

set -euo pipefail

# ---------------- CONFIGURACIÓN ----------------
DATA="4000_data.csv"     # Archivo generado por Python (x,y)
USE_DATA=1               # 1 = usar DATA, 0 = dataset sintético
EPS_LIST=(0.03)          # Puedes agregar más eps (ej: 0.03 0.05 0.08)
MINPTS=10
THREADS=(1 2 4 8 16)
REPS=3                   # Número de repeticiones por combinación
RESULTS_DIR="results"
RAW_CSV="$RESULTS_DIR/experiments_raw.csv"

mkdir -p "$RESULTS_DIR"

# Crear encabezado si no existe
if [[ ! -f "$RAW_CSV" ]]; then
  echo "impl,threads,eps,minpts,rep,time_s,clusters,core,border,noise" > "$RAW_CSV"
fi

# ---------------- FUNCIONES --------------------
append_row() {
  local impl="$1"; shift
  local threads="$1"; shift
  local eps="$1"; shift
  local minpts="$1"; shift
  local rep="$1"; shift
  local log_file="$1"; shift

  local time=$(awk -F'=' '/^time_s=/{print $2}' "$log_file" | tail -n1)
  local line=$(awk '/^clusters=/{print $0}' "$log_file" | tail -n1)
  local clusters=$(echo "$line" | awk -F'[ =]+' '{print $2}')
  local core=$(echo "$line" | awk -F'[ =]+' '{print $4}')
  local border=$(echo "$line" | awk -F'[ =]+' '{print $6}')
  local noise=$(echo "$line" | awk -F'[ =]+' '{print $8}')
  echo "$impl,$threads,$eps,$minpts,$rep,$time,$clusters,$core,$border,$noise" >> "$RAW_CSV"
}

# ---------------- CORRIDAS --------------------
for eps in "${EPS_LIST[@]}"; do
  echo "== EPS=$eps, MINPTS=$MINPTS =="

  # -------- SERIAL --------
  for rep in $(seq 1 $REPS); do
    out_file="$RESULTS_DIR/serial_eps${eps}_min${MINPTS}_rep${rep}.csv"
    if [[ "$USE_DATA" -eq 1 ]]; then
      cmd=(./serial --in "$DATA" --eps "$eps" --minpts "$MINPTS" --out "$out_file")
    else
      cmd=(./serial --eps "$eps" --minpts "$MINPTS" --out "$out_file")
    fi
    log="$RESULTS_DIR/serial_eps${eps}_min${MINPTS}_rep${rep}.log"
    echo "[serial] rep=$rep"
    "${cmd[@]}" | tee "$log"
    append_row "serial" 1 "$eps" "$MINPTS" "$rep" "$log"
  done

  # -------- OMP1 --------
  for t in "${THREADS[@]}"; do
    for rep in $(seq 1 $REPS); do
      out_file="$RESULTS_DIR/omp1_t${t}_eps${eps}_min${MINPTS}_rep${rep}.csv"
      if [[ "$USE_DATA" -eq 1 ]]; then
        cmd=(./omp1 --threads "$t" --in "$DATA" --eps "$eps" --minpts "$MINPTS" --out "$out_file")
      else
        cmd=(./omp1 --threads "$t" --eps "$eps" --minpts "$MINPTS" --out "$out_file")
      fi
      log="$RESULTS_DIR/omp1_t${t}_eps${eps}_min${MINPTS}_rep${rep}.log"
      echo "[omp1] t=$t rep=$rep"
      "${cmd[@]}" | tee "$log"
      append_row "omp1" "$t" "$eps" "$MINPTS" "$rep" "$log"
    done
  done

  # -------- OMP2 --------
  for t in "${THREADS[@]}"; do
    for rep in $(seq 1 $REPS); do
      out_file="$RESULTS_DIR/omp2_t${t}_eps${eps}_min${MINPTS}_rep${rep}.csv"
      if [[ "$USE_DATA" -eq 1 ]]; then
        cmd=(./omp2 --threads "$t" --in "$DATA" --eps "$eps" --minpts "$MINPTS" --out "$out_file")
      else
        cmd=(./omp2 --threads "$t" --eps "$eps" --minpts "$MINPTS" --out "$out_file")
      fi
      log="$RESULTS_DIR/omp2_t${t}_eps${eps}_min${MINPTS}_rep${rep}.log"
      echo "[omp2] t=$t rep=$rep"
      "${cmd[@]}" | tee "$log"
      append_row "omp2" "$t" "$eps" "$MINPTS" "$rep" "$log"
    done
  done

done

echo
echo "✅ Listo. Resultados guardados en:"
echo "   $RAW_CSV"
echo "   (usa analyze_results.py para procesarlos)"
