# Commands useful for debugging

```bash
%%bash
pg_ctl --version
ls 
cat /workspace/qcfractal/resources.yml
#pg_ctl -D /workspace/qcfractal/postgres stop -m immediate

# NOTE kill server when finished by removing the # and executing:
# !ps aux | grep qcfractal | awk '{ print $2 }' | xargs kill -9


%%bash
for pid in /proc/[0-9]*; do
    echo "${pid##*/}" `cat ${pid}/comm` >> procs.txt
done
grep -e 'postgres' ./procs.txt
!grep -e 'postgres' ./procs.txt | awk '{ print $1 }' | xargs kill -9

!ls qcfractal/
!cat ./qcfractal/qcfractal_database.log
!pidof qcfractal-server
!pidof postgresql
!top -b -n 1

```
