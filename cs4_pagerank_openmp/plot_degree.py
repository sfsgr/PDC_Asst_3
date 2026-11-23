import csv
import matplotlib.pyplot as plt

degrees = []
counts = []

with open("degree_distribution.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        degrees.append(int(row["degree"]))
        counts.append(int(row["count"]))

plt.figure(figsize=(9, 6))
plt.loglog(degrees, counts, marker='.', linestyle='none')
plt.xlabel("In-degree (log scale)")
plt.ylabel("Number of vertices (log scale)")
plt.title("PageRank Graph In-Degree Distribution")
plt.grid(True, which="both", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("degree_distribution.png", dpi=300)
plt.show()
