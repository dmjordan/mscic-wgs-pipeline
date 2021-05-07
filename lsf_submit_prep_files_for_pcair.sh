#BSUB -J jupyter
#BSUB -n 256
#BSUB -W 24:00
#BSUB -q premium
#BSUB -R rusage[mem=16G]

#BSUB -P acc_mscic

conda deactivate
ml hail/0.2.60dev
export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/

lsf-spark-submit.sh \
    --jars $HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
    --executor-memory 14G \
    --driver-memory 14G \
    /sc/arion/projects/mscic1/scripts/WGS/prep_files_for_pcair.py
