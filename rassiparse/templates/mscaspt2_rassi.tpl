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
  {{ multistate_str }}
 properties
 shift
  {{ shift }}

>> copy $Project.JobMix JOB001
&rassi
 cipr
 ejob

>> copy $Project.JobMix {{ curr_dir }}/$Project.{{ step_str }}.JobMix
{% for ms in range(1, multistate+1) %}
>> copy $Project.caspt2.molden.{{ ms }} {{ curr_dir }}/$Project.caspt2.{{ step_str }}.molden.{{ ms }}
{% endfor %}
