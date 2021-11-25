    hn=$(hostname)
    # https://stackoverflow.com/questions/9229333/how-to-get-overall-cpu-usage-e-g-57-on-linux
    cpu=$(top -bn1 | grep "Cpu(s)" | \
           sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | \
           awk '{print 100 - $1"%"}')
    mem=$(free -hm | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f 7)

    echo "${hn}   ${cpu}   ${mem}"