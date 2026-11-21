import csv
import sys

import matplotlib.pyplot as plt


def read_rows(path):
  with open(path, newline="") as file:
    return list(csv.DictReader(file))


def strong_plot(rows, output):
  ranks=[int(row["ranks"]) for row in rows]
  serial_time=float(rows[0]["t1"])
  times=[float(row["time"]) for row in rows]
  speedup=[serial_time/time for time in times]
  efficiency=[speedup[i]/ranks[i] for i in range(len(ranks))]

  fig,ax=plt.subplots()
  ax.plot(ranks,speedup,marker="o",label="Measured speedup")
  ax.plot(ranks,ranks,linestyle="--",label="Ideal speedup")
  ax.set_xlabel("MPI ranks")
  ax.set_ylabel("Speedup")
  ax.set_xticks(ranks)
  ax.grid(True)
  ax.legend()
  fig.tight_layout()
  fig.savefig(output)

  print("ranks,speedup,efficiency")
  for i,rank in enumerate(ranks):
    print(f"{rank},{speedup[i]:.6f},{efficiency[i]:.6f}")


def weak_plot(rows, output):
  ranks=[int(row["ranks"]) for row in rows]
  base_time=float(rows[0]["time"])
  times=[float(row["time"]) for row in rows]
  efficiency=[base_time/time for time in times]

  fig,ax=plt.subplots()
  ax.plot(ranks,efficiency,marker="o")
  ax.axhline(1.0,linestyle="--",label="Ideal efficiency")
  ax.set_xlabel("MPI ranks")
  ax.set_ylabel("Weak-scaling efficiency")
  ax.set_xticks(ranks)
  ax.grid(True)
  ax.legend()
  fig.tight_layout()
  fig.savefig(output)

  print("ranks,efficiency")
  for i,rank in enumerate(ranks):
    print(f"{rank},{efficiency[i]:.6f}")


if __name__=="__main__":
  if len(sys.argv)!=4:
    print("usage:python3 plot_scaling.py strong|weak CSV_FILE OUTPUT_IMAGE")
    raise SystemExit(1)

  mode=sys.argv[1]
  rows=read_rows(sys.argv[2])
  if mode=="strong":
    strong_plot(rows,sys.argv[3])
  elif mode=="weak":
    weak_plot(rows,sys.argv[3])
  else:
    print("mode must be strong or weak")
    raise SystemExit(1)
