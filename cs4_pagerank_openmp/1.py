
import csv
import matplotlib.pyplot as plt

degrees = []
counts = []

# Read degree distribution from CSV
with open("degree_distribution.csv", "r") as file:
    reader = csv.DictReader(file)

    for row in reader:
        degrees.append(int(row["degree"]))
        counts.append(int(row["count"]))

# Generate log-log plot
plt.figure(figsize=(10, 6))

plt.loglog(
    degrees,
    counts,
    marker=".",
    linestyle="none",
    markersize=5
)

plt.xlabel("In-degree (log scale)")
plt.ylabel("Number of vertices (log scale)")
plt.title("PageRank Graph In-Degree Distribution")

plt.grid(True, which="both", linestyle="--", alpha=0.4)
plt.tight_layout()

# Save plot
plt.savefig("degree_distribution.png", dpi=300)

plt.show()

print("Plot saved as degree_distribution.png")
