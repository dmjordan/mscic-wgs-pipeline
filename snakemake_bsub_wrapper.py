#!/usr/bin/env python
import shlex
import sys, subprocess, sysconfig
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

bsub_cmd = "bsub -q {queue} -P acc_{project} -W {time_min} -n {threads} -R rusage[mem={mem_mb}]".format_map(job_properties)

if job_properties.get('hail', 0) == 1:
    submit_script = shlex.quote(f"""ml spark/2.4.5 java
    ml -python
    export HAIL_HOME={sysconfig.get_path("purelib")}/hail
    export SPARK_LOCAL_DIRS=/local/tmp/
    export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
    export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/
    export LD_PRELOAD=/usr/lib64/libgslcblas.so

    lsf-spark-submit.sh \
    --jars $HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
    --conf spark.executor.extraJavaOptions='-Djava.io.tmpdir=/local/tmp' \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
    --executor-memory {job_properties['mem_mb'] - 4000} \
    --driver-memory {job_properties['mem_mb'] - 4} {jobscript}""")
else:
    submit_script = jobscript
subprocess.run(bsub_cmd + " " + submit_script)