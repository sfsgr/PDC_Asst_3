import csv
import re
import subprocess
import sys
import math


def run(command):
  result=subprocess.run(command,check=True,capture_output=True,text=True)
  line=result.stdout.strip()
  time_match=re.search(r"time=([0-9.]+)",line)
  exchange_match=re.search(r"exchange_time=([0-9.]+)",line)
  update_match=re.search(r"update_time=([0-9.]+)",line)
  if not time_match:
    raise RuntimeError(line)
  exchange=float(exchange_match.group(1)) if exchange_match else 0.0
  update=float(update_match.group(1)) if update_match else 0.0
  return float(time_match.group(1)),exchange,update,line


def median(values):
  values=sorted(values)
  middle=len(values)//2
  if len(values)%2:
    return values[middle]
  return (values[middle-1]+values[middle])/2


def mpi_command(ranks,arguments):
  return ["mpirun","--allow-run-as-root","-n",str(ranks),"./parallel"]+arguments


if __name__=="__main__":
  n=int(sys.argv[1]) if len(sys.argv)>1 else 1024
  steps=int(sys.argv[2]) if len(sys.argv)>2 else 500
  runs=int(sys.argv[3]) if len(sys.argv)>3 else 3
  ranks=[int(value) for value in sys.argv[4].split(",")] if len(sys.argv)>4 else [1,2,4,8]
  output=sys.argv[5] if len(sys.argv)>5 else "results.csv"

  serial_times=[]
  serial_line=""
  for run_index in range(runs):
    serial_times_value,serial_exchange,serial_update,serial_line=run(["./serial",str(n),str(steps),"0.1","50"])
    serial_times.append(serial_times_value)
  serial_median=median(serial_times)
  print(serial_line)

  rows=[]
  for rank in ranks:
    for decomp in ["1d","2d"]:
      for exchange in ["blocking","nonblocking"]:
        times=[]
        exchange_times=[]
        update_times=[]
        line=""
        oversubscribe="--oversubscribe" if rank>4 else ""
        for run_index in range(runs):
          command=["mpirun","--allow-run-as-root"]
          if oversubscribe:
            command.append(oversubscribe)
          command+=["-n",str(rank),"./parallel",str(n),str(steps),"0.1","50",decomp,exchange]
          time_value,exchange_time,update_time,line=run(command)
          times.append(time_value)
          exchange_times.append(exchange_time)
          update_times.append(update_time)
        median_time=median(times)
        median_exchange=median(exchange_times)
        median_update=median(update_times)
        speedup=serial_median/median_time
        efficiency=speedup/rank
        rows.append({
          "n":n,
          "t":steps,
          "ranks":rank,
          "decomp":decomp,
          "exchange":exchange,
          "serial_median":f"{serial_median:.9f}",
          "time":f"{median_time:.9f}",
          "exchange_time":f"{median_exchange:.9f}",
          "update_time":f"{median_update:.9f}",
          "speedup":f"{speedup:.6f}",
          "efficiency":f"{efficiency:.6f}",
        })

  for k in [1,50]:
    times=[]
    for run_index in range(runs):
      time_value,exchange_time,update_time,line=run(mpi_command(4,[str(n),str(steps),"0.1",str(k),"2d","blocking"]))
      times.append(time_value)
    print(f"K={k} time={median(times):.9f}")

  with open(output,"w",newline="") as file:
    fieldnames=["n","t","ranks","decomp","exchange","serial_median","time","exchange_time","update_time","speedup","efficiency"]
    writer=csv.DictWriter(file,fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

  weak_output="weak_"+output
  weak_rows=[]
  for rank in ranks:
    weak_n=max(n,int(round(n*math.sqrt(rank))))
    times=[]
    line=""
    command=["mpirun","--allow-run-as-root"]
    if rank>4:
      command.append("--oversubscribe")
    command+=["-n",str(rank),"./parallel",str(weak_n),str(steps),"0.1","50","2d","blocking"]
    for run_index in range(runs):
      time_value,exchange_time,update_time,line=run(command)
      times.append(time_value)
    median_time=median(times)
    weak_rows.append({"ranks":rank,"n":weak_n,"time":f"{median_time:.9f}"})
  with open(weak_output,"w",newline="") as file:
    writer=csv.DictWriter(file,fieldnames=["ranks","n","time"])
    writer.writeheader()
    writer.writerows(weak_rows)

  print(output)
  print(weak_output)
