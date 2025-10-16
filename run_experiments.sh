#!/bin/bash
# Script: run_experiments.sh
# Ejecuta las tres versiones (serial, omp1, omp2) para medir rendimiento y genera un CSV con resultados robusto. Además de un archivo de clasificación para omp2.

EPS=0.03
MINPTS=10
REPS=10

CORES_VIRTUALES=$(sysctl -n hw.logicalcpu)
THREADS_LIST=(1 $(($CORES_VIRTUALES/2)) $CORES_VIRTUALES $(($CORES_VIRTUALES*2)))
POINTS_LIST=(20000 40000 80000 120000 140000 160000 180000 200000)

RENDIMIENTO_DIR="output/rendimiento"
CLASIFICACION_DIR="output/clasificacion"
mkdir -p "$RENDIMIENTO_DIR"
mkdir -p "$CLASIFICACION_DIR"

for N in "${POINTS_LIST[@]}"; do
    DATA="input/${N}_data.csv"
    OUT_CSV="$RENDIMIENTO_DIR/${N}_results.csv"
    
    echo ">>> Dataset: $DATA"
    echo "impl,threads,eps,minpts,rep,time_s,clusters,core,border,noise" > "$OUT_CSV"

    for REP in $(seq 1 $REPS); do
        # ---------------- Serial ----------------
        OUTPUT=$(./src/serial --in "$DATA" --eps "$EPS" --minpts "$MINPTS")
        TIME=$(echo "$OUTPUT" | grep "time_s=" | cut -d= -f2)
        CLUSTERS=$(echo "$OUTPUT" | awk '/clusters=/ {for(i=1;i<=NF;i++) if($i ~ /^clusters=/){split($i,a,"="); print a[2]}}')
        CORE=$(echo "$OUTPUT" | awk '/core=/ {for(i=1;i<=NF;i++) if($i ~ /^core=/){split($i,a,"="); print a[2]}}')
        BORDER=$(echo "$OUTPUT" | awk '/border=/ {for(i=1;i<=NF;i++) if($i ~ /^border=/){split($i,a,"="); print a[2]}}')
        NOISE=$(echo "$OUTPUT" | awk '/noise=/ {for(i=1;i<=NF;i++) if($i ~ /^noise=/){split($i,a,"="); print a[2]}}')
        echo "serial,1,$EPS,$MINPTS,$REP,$TIME,$CLUSTERS,$CORE,$BORDER,$NOISE" >> "$OUT_CSV"

        # ---------------- omp1 ----------------
        for T in "${THREADS_LIST[@]}"; do
            OUTPUT=$(./src/omp1 --threads $T --in "$DATA" --eps "$EPS" --minpts "$MINPTS")
            TIME=$(echo "$OUTPUT" | grep "time_s=" | cut -d= -f2)
            CLUSTERS=$(echo "$OUTPUT" | awk '/clusters=/ {for(i=1;i<=NF;i++) if($i ~ /^clusters=/){split($i,a,"="); print a[2]}}')
            CORE=$(echo "$OUTPUT" | awk '/core=/ {for(i=1;i<=NF;i++) if($i ~ /^core=/){split($i,a,"="); print a[2]}}')
            BORDER=$(echo "$OUTPUT" | awk '/border=/ {for(i=1;i<=NF;i++) if($i ~ /^border=/){split($i,a,"="); print a[2]}}')
            NOISE=$(echo "$OUTPUT" | awk '/noise=/ {for(i=1;i<=NF;i++) if($i ~ /^noise=/){split($i,a,"="); print a[2]}}')
            echo "omp1,$T,$EPS,$MINPTS,$REP,$TIME,$CLUSTERS,$CORE,$BORDER,$NOISE" >> "$OUT_CSV"
        done

        # ---------------- omp2 ----------------
        for T in "${THREADS_LIST[@]}"; do
            # SOLO guardar clasificación en primera repetición con máximo de hilos
            SAVE_CLASS=""
            if [ "$REP" -eq 1 ] && [ "$T" -eq "$CORES_VIRTUALES" ]; then
                SAVE_CLASS="--out $CLASIFICACION_DIR/omp2_${N}.csv"
            fi
            
            OUTPUT=$(./src/omp2 --threads $T --in "$DATA" --eps "$EPS" --minpts "$MINPTS" $SAVE_CLASS)
            TIME=$(echo "$OUTPUT" | grep "time_s=" | cut -d= -f2)
            CLUSTERS=$(echo "$OUTPUT" | awk '/clusters=/ {for(i=1;i<=NF;i++) if($i ~ /^clusters=/){split($i,a,"="); print a[2]}}')
            CORE=$(echo "$OUTPUT" | awk '/core=/ {for(i=1;i<=NF;i++) if($i ~ /^core=/){split($i,a,"="); print a[2]}}')
            BORDER=$(echo "$OUTPUT" | awk '/border=/ {for(i=1;i<=NF;i++) if($i ~ /^border=/){split($i,a,"="); print a[2]}}')
            NOISE=$(echo "$OUTPUT" | awk '/noise=/ {for(i=1;i<=NF;i++) if($i ~ /^noise=/){split($i,a,"="); print a[2]}}')
            echo "omp2,$T,$EPS,$MINPTS,$REP,$TIME,$CLUSTERS,$CORE,$BORDER,$NOISE" >> "$OUT_CSV"
        done
    done
done

echo "✅ Experimentos completados."
echo "✅ Resultados de rendimiento en: $RENDIMIENTO_DIR/"
echo "✅ Archivos de clasificación (omp2) en: $CLASIFICACION_DIR/"