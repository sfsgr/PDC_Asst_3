import csv
import re
import statistics
import subprocess
import sys


def run(length,k,mode,threads,skewed=0,runs=3):
  values=[]
  output=""
  for i in range(runs):
    completed=subprocess.run(
      ["./kmer",str(length),str(k),mode,str(threads),"12345",str(skewed)],
      check=True,capture_output=True,text=True)
    output=completed.stdout.strip()
    match=re.search(r"time=([0-9.]+)",output)
    values.append(float(match.group(1)))
  return statistics.median(values),output


if __name__=="__main__":
  length=int(sys.argv[1]) if len(sys.argv)>1 else 100000000
  runs=int(sys.argv[2]) if len(sys.argv)>2 else 3
  output=sys.argv[3] if len(sys.argv)>3 else "kmer_results.csv"
  rows=[]
  for k in [1,8,12]:
    for threads in [1,2,4,8]:
      for mode in ["serial","race","atomic","private","padded"]:
        if mode=="serial" and threads!=1:
          continue
        if mode=="padded" and k not in [1,8]:
          continue
        median_time,line=run(length,k,mode,threads,0,runs)
        match=re.search(r"throughput=([0-9.]+)M",line)
        throughput=float(match.group(1)) if match else 0.0
        rows.append({
          "N":length,
          "k":k,
          "mode":mode,
          "threads":threads,
          "time":f"{median_time:.9f}",
          "throughput":f"{throughput:.3f}",
          "correct": "1" if " correct=1 " in " "+line+" " else "0",
        })
        print(line)
  with open(output,"w",newline="") as file:
    writer=csv.DictWriter(file,fieldnames=["N","k","mode","threads","time","throughput","correct"])
    writer.writeheader()
    writer.writerows(rows)
  print(output)
