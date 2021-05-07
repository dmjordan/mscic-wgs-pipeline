#!/bin/bash

#ml hail/0.2.60dev
SAVED_ENV=$CONDA_DEFAULT_ENV
ml spark/2.4.5 anaconda3
unset PYTHONPATH
ml -python
source $(conda info --base)/etc/profile.d/conda.sh
conda activate
conda activate $SAVED_ENV
export HAIL_HOME=$(python -c 'import sysconfig; print(sysconfig.get_path("purelib"))')/hail

bsub -n 48 \
    -W 12:00 \
    -q private \
    -m schade01-1 \
    -R rusage[mem=8G] \
    -P acc_mscic \
    -I <<EOF
export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/

lsf-spark-submit.sh \
    --jars $HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
    --executor-memory 5G \
    --driver-memory 5G \
    $@
EOF

