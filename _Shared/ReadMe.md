# ReadMe.md

**Table of Contents**

[TOCM]

[TOC]

# MultiRunModelUpdating.m
The script can perform multiple runs of model updating by starting the search from randomized starting points.
## Syntax
MultiRunModelUpdating
## Description
Before running this script, the workspace must have following variables before running the script. See documentation of StructModelUpdating.m for descriptions.

**structModel, expModes, updatingOpts, optimzOpts, filename**

If following variables do not exist, default actions will be assigned.

**numRuns**: number of starting points.  Default = 1
**randSeed**: random seed value for rng.  By default the random seed will NOT be fixed.

If multiple runs are requested through "numRuns" variable and the runs are disrupted half way, the script can pick up and continue the runs by loading previous run results from the .mat named by the "filename" variable. **IF A USER WANTS TO START ALL THE RUNS AFRESH, THE .MAT FILE MUST BE DELETED FIRST.**

