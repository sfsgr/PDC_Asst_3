import csv
import re
import statistics
import subprocess
import sys


def run(rank,n,distribution,oversample,runs):
  times=[]
  line=""
  for i in range(runs):
    command=["mpirun","--allow-run-as-root","-n",str(rank),"./parallel_sample_sort",
             str(n),str(distribution),str(oversample)]
    if rank>4:
      command=["mpirun","--allow-run-as-root","--oversubscribe","-n",str(rank),
               "./parallel_sample_sort",str(n),str(distribution),str(oversample)]
    result=subprocess.run(command,check=True,capture_output=True,text=True)
    line=result.stdout.strip()
    match=re.search(r"Total Time:\s*([0-9.]+)",line)
    times.append(float(match.group(1)))
  return statistics.median(times),line


if __name__=="__main__":
  n=int(sys.argv[1]) if len(sys.argv)>1 else 100000000
  runs=int(sys.argv[2]) if len(sys.argv)>2 else 3
  output=sys.argv[3] if len(sys.argv)>3 else "sample_sort_results.csv"
  rows=[]
  ranks=[1,2,4,8]
  for distribution in [0,1]:
    for oversample in [1,4]:
      for rank in ranks:
        median_time,line=run(rank,n,distribution,oversample,runs)
        fields={
          "N":n,
          "distribution":"uniform" if distribution==0 else "skewed",
          "ranks":rank,
          "oversample":oversample,
          "time":f"{median_time:.9f}",
        }
        for name,label in [
          ("local","Local Sort"),
          ("sampling","Sampling"),
          ("partition","Partition"),
          ("alltoall","Alltoall Counts"),
          ("alltoallv","Alltoallv Data"),
          ("final","Final Sort"),
          ("min","Min Elements"),
          ("max","Max Elements"),
          ("imbalance","Imbalance"),
        ]:
          match=re.search(r"%s:\s*([0-9.]+)"%label,line)
          fields[name]=match.group(1) if match else ""
        rows.append(fields)
        print(line)
  fieldnames=["N","distribution","ranks","oversample","time","local","sampling","partition",
              "alltoall","alltoallv","final","min","max","imbalance"]
  with open(output,"w",newline="") as file:
    writer=csv.DictWriter(file,fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
  print(output)
