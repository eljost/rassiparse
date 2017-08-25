#!/bin/bash
# args:
#   $1 queue

{% for step_str in step_strs %}
    cd {{ step_str }}
    sbatch -p $1 submolcas.sh
    cd ..
{% endfor %}
