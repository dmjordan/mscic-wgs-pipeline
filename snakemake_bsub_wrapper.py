#!/usr/bin/env python
import shlex, re
import sys, subprocess, sysconfig
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

if job_properties['type'] == "group":
    job_properties['rulename'] = job_properties['groupid']
else:
    job_properties['rulename'] = job_properties['rule']

bsub_cmd = "bsub -q {resources[queue]} -P acc_{resources[project]} -W {resources[time_min]} -n {resources[cpus]} -R rusage[mem={resources[mem_mb]}] -oo {rulename}.%J.log".format_map(job_properties)
if job_properties['resources'].get("single_host", 0) == 1:
    bsub_cmd += " -R span[hosts=1]"
if job_properties['resources'].get('host') is not None:
    bsub_cmd += " -m {resources[host]}".format_map(job_properties)
output = subprocess.check_output(bsub_cmd + " " + jobscript, shell=True, text=True)
match = re.match(f"Job <(\\d+)> is submitted to queue <{job_properties['resources']['queue']}>.", output)
if match is not None:
    print(match.group(1))
