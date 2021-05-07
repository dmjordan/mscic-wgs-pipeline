#BSUB -J hello_hail
#BSUB -n 8 
#BSUB -W 1:00
#BSUB -oo hello_hail.%J.out

#BSUB -P acc_mscic

conda deactivate
ml hail/0.2.60dev

export SPARK_LOG_DIR="spark_log/"
export SPARK_WORKER_DIR="spark_work/"

lsf-spark-submit.sh \
    --jars $HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
    /sc/arion/projects/mscic1/scripts/WGS/hello_hail.py
