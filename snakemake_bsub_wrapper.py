#!/usr/bin/env python
import shlex, re
import sys, subprocess, sysconfig
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

bsub_cmd = "bsub -q {resources[queue]} -P acc_{resources[project]} -W {resources[time_min]} -n {resources[cpus]} -R rusage[mem={resources[mem_mb]}] -oo {rule}.%J.log".format_map(job_properties)
if job_properties['resources'].get("single_host", 0) == 1:
    bsub_cmd += " -R span[hosts=1]"
output = subprocess.check_output(bsub_cmd + " " + jobscript, shell=True, text=True)
match = re.match(f"Job <(\\d+)> is submitted to queue <{job_properties['resources']['queue']}>.", output)
if match is not None:
    print(match.group(1))
