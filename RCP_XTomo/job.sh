#!/bin/ksh
# @ job_type = serial
# @ executable = a.out
# @ arguments = a b c
# @ output = $(executable).$(jobid).$(stepid).out
# @ error = $(executable).$(jobid).$(stepid).err
# @ initialdir = 
# @ notify_user = gan@imf.fzk.de
# @ notification = complete
# @ class = public
# @ environment = COPY_ALL
# @ requirements = (Arch == "R6000") && (OpSys == "AIX53")
# @ queue
