#!/usr/bin/env python

import sys, time
from subprocess import Popen, PIPE
import pandas as pd

jobid = sys.argv[1]

def get_bjobs(jobid):
    bjobs_process = Popen(["bjobs", jobid], text=True, stdout=PIPE)
    try:
        table = pd.read_fwf(bjobs_process.stdout)
    except pd.errors.EmptyDataError:
        return None
    else:
        return table.STAT[0]


for i in range(60):
    result = get_bjobs(jobid)
    if result is not None:
        break


if result in {"PEND", "RUN"}:
    print("running")
elif result == "DONE":
    print("success")
else:
    print("failed")
