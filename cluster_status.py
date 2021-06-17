#!/usr/bin/env python

import sys
from subprocess import Popen, PIPE
import pandas as pd

jobid = sys.argv[1]

bjobs_process = Popen(["bjobs", jobid], text=True, stdout=PIPE)
table = pd.read_fwf(bjobs_process.stdout)

if table.STAT[0] in {"PEND", "RUN"}:
    print("running")
elif table.STAT[0] == "DONE":
    print("success")
else:
    print("failed")
