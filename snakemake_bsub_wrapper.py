#!/usr/bin/env python
import shlex, re
import sys, subprocess, sysconfig
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

bsub_cmd = "bsub -q {resources[queue]} -P acc_{resources[project]} -W {resources[time_min]} -n {resources[cpus]} -R rusage[mem={resources[mem_mb]}] -oo {rule}.%J.log".format_map(job_properties)
if job_properties['resources'].get("single_host", 0) == 1:
    bsub_cmd += " -R span[hosts=1]"


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
output = subprocess.check_output(bsub_cmd + " " + submit_script, shell=True, text=True)
match = re.match(f"Job <(\\d+)> is submitted to queue <{job_properties['resources']['queue']}>.", output)
if match is not None:
    print(match.group(1))
