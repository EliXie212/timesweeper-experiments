import os
import sys
from initializeVar import *

sys.path.insert(1, '/pine/scr/e/m/emae/timeSeriesSweeps')

import runCmdAsJobTEST

prefixLs = ['hard_soft_neut_ttv_ali', 'hard_soft_neut_ttv_haps', 'hard_soft_neut_ttv_sfs', 'hard_soft_neut_ttv_jsfs']

simTypeToScript = {"":"../../keras_CNN_loadNrun.py", "1Samp":"../../keras_DNN_loadNrun.py"}
for simType in ["", "1Samp"]:
    outDir="{}/classifiers{}".format(baseDir, simType)
    os.system("mkdir -p {}".format(outDir))

    for prefix in prefixLs:
        cmd = "python {0} -i {1}/npzs{2}/{3}.npz -c {4}/{3}.mod".format(simTypeToScript[simType], baseDir, simType, prefix, outDir)
        runCmdAsJobTEST.runCmdAsJobWithoutWaitingWithLog(cmd, "trainTS", "trainTS.txt", "12:00:00", "general", "32GB", outDir+"/{}.log".format(prefix))
