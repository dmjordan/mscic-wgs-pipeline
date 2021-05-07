#!/bin/bash

SAVED_ENV=$CONDA_DEFAULT_ENV
ml spark/2.4.5 anaconda3
unset PYTHONPATH
source $(conda info --base)/etc/profile.d/conda.sh
conda activate
conda activate $SAVED_ENV
export HAIL_HOME=$(python -c 'import sysconfig; print(sysconfig.get_path("purelib"))')/hail

export SPARK_MASTER_PORT=6311
export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/
#export SPARK_WORKER_INSTANCES=48
export SPARK_WORKER_CORES=48

$SPARK_HOME/sbin/start-master.sh -p $SPARK_MASTER_PORT
$SPARK_HOME/sbin/start-slave.sh spark://$(hostname):$SPARK_MASTER_PORT

lsf-spark-submit.sh \
    --jars $HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
    --executor-memory 512G \
    --driver-memory 512G \
    $@
