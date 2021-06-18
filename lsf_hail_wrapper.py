import subprocess
import sysconfig, os, sys
import typing

# declaring the snakemake object so my IDE stops yelling at me
if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

new_env_variables = {
    'HAIL_HOME': os.path.join(sysconfig.get_path("purelib"), "hail"),
    "SPARK_LOCAL_DIRS": "/local/tmp",
    "SPARK_LOG_DIR": "/sc/arion/projects/mscic1/scratch/hail/logs/",
    "SPARK_WORKER_DIR": "/sc/arion/projects/mscic1/scratch/hail/worker/",
    "LD_PRELOAD": "/usr/lib64/libgslcblas.so"
}

env = os.environ.copy()
env.update(new_env_variables)

hail_args = list(snakemake.input)
if snakemake.params.get("pass_output", False):
    hail_args += list(snakemake.output)
hail_args_formatted = " ".join(hail_args)

hail_script_path = os.path.join(snakemake.config["scriptsdir"], "hail_wgs.py")

spark_submit_script = \
f"""ml spark/2.4.5 java
ml -python
lsf-spark-submit.sh \
--jars $HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.executor.extraJavaOptions='-Djava.io.tmpdir=/local/tmp' \
--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
--executor-memory {snakemake.resources.mem_mb - 4000}M \
--driver-memory {snakemake.resources.mem_mb - 4000}M \
{hail_script_path} {snakemake.params.hail_cmd} {hail_args_formatted}"""

hail_process = subprocess.run(spark_submit_script, env=env, shell=True)
sys.exit(hail_process.returncode)
