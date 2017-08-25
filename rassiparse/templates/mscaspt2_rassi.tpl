&gateway
 coord
  {{ coord_fn }}
 basis
  {{ basis }}
 ricd

&seward

>> copy {{ jobiph_fn }} $Project.JobIph

&caspt2
 multistate
  {{ multistate }}
 properties
 shift
  {{ shift }}

>> copy $Project.JobMix JOB001
&rassi
 cipr
 ejob

>> copy $Project.JobMix {{ curr_dir }}/$Project.{{ step_str }}.JobMix
