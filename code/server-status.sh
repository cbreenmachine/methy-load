#!/bin/bash
RUN_SERVERS=$(echo nebula-{1..10} | sed -e 's/\s\+/,/g')
parallel --nonall -S ${RUN_SERVERS} --workdir . 'bash summary.sh'