import os
import csv
import statistics
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image


ROOT=r"\\wsl.localhost\Ubuntu\home\tayyab\PDC_Asst_3"
REPORT=os.path.join(ROOT,"report")
ASSETS=os.path.join(REPORT,"assets")
PDF=os.path.join(REPORT,"GroupXX_report.pdf")
os.makedirs(ASSETS,exist_ok=True)

plt.rcParams.update({"font.size":9})


def read_csv(path):
    with open(path,newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def save_heat_graphs():
    path=os.path.join(ROOT,"cs1_heat_mpi","strong_results.csv")
    rows=read_csv(path)
    ranks=[]
    speed=[]
    efficiency=[]
    for row in rows:
        if row["decomp"]=="2d" and row["exchange"]=="blocking":
            ranks.append(int(row["ranks"]))
            speed.append(float(row["speedup"]))
            efficiency.append(float(row["efficiency"]))
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.plot(ranks,speed,"o-",label="2-D blocking")
    ax.plot(ranks,[1,2,4,8],"--",label="Ideal")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Speedup")
    ax.set_title("Case 1: Strong scaling, median of 3 runs")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"heat_strong.png"),dpi=180)
    plt.close(fig)

    weak_path=os.path.join(ROOT,"cs1_heat_mpi","weak_strong_results.csv")
    weak=read_csv(weak_path)
    wr=[int(row["ranks"]) for row in weak]
    wt=[float(row["time"]) for row in weak]
    eff=[wt[0]/value for value in wt]
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.plot(wr,eff,"o-",label="Measured")
    ax.axhline(1.0,linestyle="--",label="Ideal")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Weak-scaling efficiency")
    ax.set_title("Case 1: Weak scaling")
    ax.set_xticks(wr)
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"heat_weak.png"),dpi=180)
    plt.close(fig)


def save_kmer_graph():
    rows=read_csv(os.path.join(ROOT,"cs2_kmer_openmp","kmer_100m_results.csv"))
    modes=["serial","race","atomic","private","padded"]
    k8=[row for row in rows if int(row["k"])==8]
    values={row["mode"]:float(row["throughput"]) for row in k8}
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.bar([m for m in modes if m in values],[values[m] for m in modes if m in values])
    ax.set_ylabel("Million k-mers/s")
    ax.set_title("Case 2: k=8, 100 million characters")
    ax.grid(True,axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"kmer.png"),dpi=180)
    plt.close(fig)


def save_mc_graph():
    serial=statistics.median([10.947842,10.852722,9.757764])
    times=[serial,statistics.median([12.146964,10.070304,10.602978]),
           statistics.median([6.467711,6.842922,6.002942]),
           statistics.median([2.616511,3.740699,4.327151]),
           statistics.median([3.604589,2.978711,2.628556])]
    threads=[1,1,2,4,8]
    speed=[serial/times[i] for i in range(len(times))]
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.plot([1,2,4,8],speed[1:],"o-",label="Measured")
    ax.plot([1,2,4,8],[1,2,4,8],"--",label="Ideal")
    ax.set_xlabel("OpenMP threads")
    ax.set_ylabel("Speedup")
    ax.set_title("Case 3: Strong scaling, median of 3 runs")
    ax.set_xticks([1,2,4,8])
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"mc.png"),dpi=180)
    plt.close(fig)


def save_pagerank_graphs():
    degree_path=os.path.join(ROOT,"cs4_pagerank_openmp","degree_distribution.csv")
    degree=read_csv(degree_path)
    d=[int(row["degree"]) for row in degree]
    c=[int(row["count"]) for row in degree]
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.loglog(d,c,".",markersize=3)
    ax.set_xlabel("In-degree")
    ax.set_ylabel("Vertex count")
    ax.set_title("Case 4: R-MAT in-degree distribution")
    ax.grid(True,which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"pagerank_degree.png"),dpi=180)
    plt.close(fig)

    sched=["static64","static256","dynamic64","dynamic256","guided64","guided256"]
    throughput=[809.44,879.57,1995.15,1825.27,824.50,839.75]
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.bar(sched,throughput)
    ax.set_ylabel("Million edges/s")
    ax.set_title("Case 4: Schedule throughput")
    ax.tick_params(axis="x",rotation=35)
    ax.grid(True,axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"pagerank_schedule.png"),dpi=180)
    plt.close(fig)


def save_sample_sort_graphs():
    ranks=[1,2,4,8]
    uniform=[23.5246,20.1410,14.0665,13.0726]
    skewed=[28.0142,20.6089,13.5135,12.8004]
    serial_u=19.2749
    serial_s=19.3254
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.plot(ranks,[serial_u/u for u in uniform],"o-",label="Uniform")
    ax.plot(ranks,[serial_s/u for u in skewed],"s-",label="Skewed")
    ax.plot(ranks,[1,2,4,8],"--",label="Ideal")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Speedup")
    ax.set_title("Case 5: Strong scaling")
    ax.set_xticks(ranks)
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"sample_scaling.png"),dpi=180)
    plt.close(fig)

    labels=["1","2","4","8"]
    local=[19.3807,15.8008,10.1680,9.8830]
    sampling=[0.0,0.0886,0.3020,2.7392]
    alltoallv=[0.2071,0.2296,0.1453,0.0822]
    final=[3.9368,4.1105,3.7528,3.1079]
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.bar(labels,local,label="Local sort")
    ax.bar(labels,sampling,bottom=local,label="Sampling")
    bottoms=[local[i]+sampling[i] for i in range(4)]
    ax.bar(labels,alltoallv,bottom=bottoms,label="Alltoallv")
    bottoms=[bottoms[i]+alltoallv[i] for i in range(4)]
    ax.bar(labels,final,bottom=bottoms,label="Final sort")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Case 5: Uniform phase breakdown")
    ax.grid(True,axis="y")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"sample_phase.png"),dpi=180)
    plt.close(fig)


def build_report():
    save_heat_graphs()
    save_kmer_graph()
    save_mc_graph()
    save_pagerank_graphs()
    save_sample_sort_graphs()

    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCustom",parent=styles["Title"],fontSize=24,leading=29,alignment=TA_CENTER,spaceAfter=18))
    styles.add(ParagraphStyle(name="Subtitle",parent=styles["Normal"],fontSize=11,leading=15,alignment=TA_CENTER,textColor=colors.HexColor("#444444")))
    styles.add(ParagraphStyle(name="Section",parent=styles["Heading1"],fontSize=16,leading=19,textColor=colors.HexColor("#1f4e79"),spaceBefore=4,spaceAfter=8))
    styles.add(ParagraphStyle(name="Subsection",parent=styles["Heading2"],fontSize=11,leading=14,textColor=colors.HexColor("#1f4e79"),spaceBefore=5,spaceAfter=4))
    styles.add(ParagraphStyle(name="BodySmall",parent=styles["BodyText"],fontSize=8.8,leading=12))

    doc=SimpleDocTemplate(PDF,pagesize=letter,rightMargin=42,leftMargin=42,topMargin=42,bottomMargin=42,
                          title="PDC Assignment 3 Report",author="sfsgr")
    story=[]

    def p(text,style="BodySmall"):
        story.append(Paragraph(text,styles[style]))
    def image(name,width=6.9,height=3.45):
        story.append(Image(os.path.join(ASSETS,name),width=width,height=height))
    def table(data,widths=None):
        t=Table(data,colWidths=widths,hAlign="LEFT")
        t.setStyle(TableStyle([
            ("FONTNAME",(0,0),(-1,-1),"Helvetica"),
            ("FONTSIZE",(0,0),(-1,-1),7.5),
            ("LEADING",(0,0),(-1,-1),9),
            ("GRID",(0,0),(-1,-1),0.25,colors.grey),
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#d9e7f5")),
            ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ]))
        story.append(t)

    story.append(Spacer(1,1.2*inch))
    p("PDC Assignment 3","TitleCustom")
    p("Parallel and Distributed Computing Case Studies","Subtitle")
    story.append(Spacer(1,0.25*inch))
    p("Group XX","Subtitle")
    p("Team account: sfsgr","Subtitle")
    p("November 2026","Subtitle")
    p("Heat diffusion - K-mer counting - Monte Carlo - PageRank - Sample sort","Subtitle")
    story.append(PageBreak())

    p("1. Introduction","Section")
    p("This report evaluates five parallel case studies. OpenMP is used for k-mer counting and PageRank. MPI is used for heat diffusion and sample sorting. Monte Carlo is implemented with OpenMP.")
    p("<b>Machine.</b> Intel Core i5-8250U, 4 physical cores and 8 logical cores, 7.7 GiB RAM, Microsoft WSL 2 Ubuntu. The WSL view exposes 8 logical CPUs.")
    p("<b>Toolchain.</b> GCC 15.2.0, GNU Make, OpenMP, OpenMPI/MPICH-compatible mpicc and mpirun, and Python for plotting.")
    p("<b>Timing method.</b> Each graph uses the median of at least 3 runs where available. The serial Monte Carlo and sample-sort timings are the same machine baseline. No machine was intentionally loaded by us. MPI p=8 is oversubscribed because only 4 physical cores are available.")
    p("<b>Repository evidence.</b> Source code, Makefiles, README instructions, and measurement CSV/text files are stored in the PDC_Asst_3 repository.")
    p("Team contribution and commit history","Subsection")
    table([
        ["Area","Main work"],
        ["Case 1","MPI decomposition, halo exchange, timing"],
        ["Case 2","OpenMP race, atomic, private and padding modes"],
        ["Case 3","Independent streams, Monte Carlo scaling"],
        ["Case 4","CSR R-MAT, pull/push and scheduling"],
        ["Case 5","MPI sample sort, ping-pong and verification"],
    ],[1.4*inch,5.3*inch])
    p("Common verification","Subsection")
    p("Integer and total-count checks use exact equality. Floating-point PageRank and Monte Carlo results use stated tolerances.")
    story.append(PageBreak())

    p("2. Case 1: Heat Diffusion","Section")
    p("The 5-point stencil updates every interior cell from the previous grid. It is parallelizable because a rank can update a block independently after its four neighbours are exchanged.")
    p("<b>Model.</b> MPI. The 1-D version divides rows; the 2-D version uses MPI_Cart_create. Blocking exchange uses MPI_Sendrecv. Non-blocking exchange uses Irecv/Isend and then updates the interior while waiting.")
    p("A 2-D rank has about N^2/p computed cells and 2N/sqrt(p) halo values in a balanced square decomposition. 1-D has 2N halo values per rank. The 2-D form therefore lowers communication volume.")
    image("heat_strong.png")
    table([
        ["Ranks","2-D blocking time","Speedup","Efficiency"],
        ["1","5.743 s","1.08","108%"],
        ["2","4.049 s","1.53","76.5%"],
        ["4","3.999 s","1.55","38.7%"],
        ["8","9.793 s","0.63","7.9%"],
    ],[1*inch,2*inch,1.5*inch,1.5*inch])
    p("<b>Result.</b> Efficiency dropped below 70% at p=4. Strong scaling stopped because communication became comparable to the small local computation. p=8 was oversubscribed.")
    story.append(PageBreak())

    p("Case 1: Analysis and Weak Scaling","Section")
    image("heat_weak.png")
    table([
        ["Ranks","N","Time","Weak efficiency"],
        ["1","1024","5.740 s","100%"],
        ["2","1448","6.743 s","85.1%"],
        ["4","2048","12.031 s","47.7%"],
        ["8","2896","29.744 s","19.3%"],
    ],[1*inch,1.2*inch,1.7*inch,2*inch])
    p("<b>Overlap.</b> For 2-D exchange, blocking/non-blocking communication times were 0.558/0.803 s at p=2, 1.965/1.451 s at p=4, and 8.631/3.054 s at p=8. At p=8, about 64% of measured communication time was hidden. Little gain was seen at p=2.")
    p("<b>Convergence checks.</b> At p=4, K=1 took 5.002 s and K=50 took 3.534 s. Frequent Allreduce synchronizes every step. K=50 is faster because it checks only once every 50 steps.")
    p("<b>Correctness.</b> Serial and parallel sums matched within 1e-9 for tested runs.")
    story.append(PageBreak())

    p("3. Case 2: K-mer Counting","Section")
    p("A k-mer is encoded in 2-bit base values. The histogram loop is data parallel, but updates to shared counters conflict.")
    p("<b>Race mode.</b> It used no protection. In the 100-million-character k=8 test, it found 99,946,568 counts instead of 99,999,993.")
    p("<b>Fixed modes.</b> Atomic increments protect each histogram entry. Private histograms give each thread a separate copy and merge once. Padded private histograms separate each copy by cache-line size.")
    image("kmer.png")
    table([
        ["k","Mode","4 threads","Result"],
        ["8","Atomic","61.966 M/s","Correct"],
        ["8","Private","237.719 M/s","Correct"],
        ["8","Padded","248.884 M/s","Correct"],
        ["12","Atomic","51.713 M/s","Correct"],
        ["12","Private","24.718 M/s","Correct"],
    ],[0.7*inch,1.2*inch,1.5*inch,1.3*inch])
    p("<b>Memory crossover.</b> Private k=8 storage is 4^8 x 4 threads x 8 bytes = about 2 MiB. Private k=12 storage is 512 MiB for 4 threads. K=12 spills memory and atomic is faster.")
    p("<b>False sharing.</b> For k=1 and 2 threads, private mode took 0.1079 s and padded took 0.0716 s, a 50.7% improvement.")
    p("<b>Cost model.</b> About 9 bytes are read or updated per k-mer in this simple estimate. Padded k=8 reached 248.884 million k-mers/s.")
    story.append(PageBreak())

    p("4. Case 3: Monte Carlo Option Pricing","Section")
    p("European option paths use independent geometric-Brownian-motion streams. The serial and parallel models use M>=100 million fixed-seed paths and the Black-Scholes price 10.4505835722.")
    p("OpenMP was chosen because all threads run independent local sums and only reduce them at the end.")
    table([
        ["Threads","Median time","Median speedup","Median throughput"],
        ["Serial","10.8527 s","1.00","9.21M/s"],
        ["1","10.6030 s","1.02","9.43M/s"],
        ["2","6.4677 s","1.68","15.46M/s"],
        ["4","3.7407 s","2.90","26.73M/s"],
        ["8","2.9787 s","3.64","33.57M/s"],
    ],[1*inch,1.6*inch,1.6*inch,1.7*inch])
    image("mc.png")
    p("<b>RNG result.</b> A bad same-seed run for 1 million paths with 4 threads gave a price of 10.512459 and error 0.061875. Independent per-thread streams gave 10.468089 and error 0.014734.")
    p("At fixed 10 seconds, p=1 produced about 90.78 million samples and error 0.002313. p=8 produced about 361.34 million samples and error 0.000356. Speedup here means more accurate samples per second.")
    p("A simple Amdahl estimate gives 13-17% serial fraction. Real losses come from reduction, thread startup, turbo changes and 8-thread SMT.")
    p("Floating-point totals can differ slightly with thread count because local sum order and reduction grouping change.")
    story.append(PageBreak())

    p("5. Case 4: Parallel PageRank","Section")
    p("The 65,536-vertex R-MAT graph had 1,335,037 edges. It was stored in CSR and transposed CSR for pull updates. Dangling-node mass was returned to all vertices.")
    p("<b>Pull vs. push.</b> Pull gives each thread a destination vertex set and no write conflict. Push gives each source an atomic update to its neighbours.")
    table([
        ["Version","Median time","Million edges/s","Rank sum"],
        ["Serial","0.0440 s","425.21","1.0000000000"],
        ["Pull, 4 threads","0.0257 s","726.35","1.0000000000"],
        ["Push, 4 threads","0.1188 s","157.31","1.0000000000"],
    ],[1.8*inch,1.6*inch,1.7*inch,1.4*inch])
    p("Pull was 1.71 times faster than serial and 4.62 times faster than push. The extra push cost came mainly from atomic random writes.")
    image("pagerank_degree.png")
    p("<b>Irregularity.</b> Static imbalance was R=1.066 in a recent 3-run median; the slowest thread had 554,057 edges and the least-loaded thread had 20,369. Equal vertex counts did not mean equal work.")
    story.append(PageBreak())

    p("Case 4: Scheduling and Memory","Section")
    image("pagerank_schedule.png")
    table([
        ["Schedule","Chunk","Time","Million edges/s"],
        ["static","64","0.0247 s","809.44"],
        ["static","256","0.0254 s","797.02"],
        ["dynamic","64","0.0118 s","1,695.96"],
        ["dynamic","256","0.0116 s","1,483.89"],
        ["guided","64","0.0258 s","777.32"],
        ["guided","256","0.0281 s","712.37"],
    ],[1.2*inch,0.8*inch,1.4*inch,2*inch])
    p("<b>Schedule result.</b> Dynamic was best because it gave difficult high-degree vertices to available threads. Chunk 64 had lower imbalance. Guided and static were close for small chunks but dynamic changed the throughput most.")
    p("<b>Memory behaviour.</b> Pull reads rank, CSR offsets and incoming neighbours. A simple upper estimate is 4-8 bytes of rank/edge data plus index traffic per edge. At 1.696 billion edges/s for dynamic scheduling, the result is close to a memory-constrained graph walk, but WSL and SMT limit direct bandwidth comparison.")
    p("<b>Correctness.</b> Serial, pull and push all gave rank sum 1.0000000000 and residual 4.445e-8. The pull/push rank differences matched within 1e-7.")
    story.append(PageBreak())

    p("6. Case 5: Distributed Sample Sort","Section")
    p("Each rank generated and locally sorted a shard. The global splitters were collected, broadcast, and used for Alltoall and Alltoallv redistribution.")
    p("<b>Validation.</b> Every run reported locally sorted data, preserved total count, and max(rank i)<=min(rank i+1).")
    p("<b>Load balance.</b>")
    table([
        ["Dist.","p","Min","Max","Max/mean"],
        ["Uniform","1","100,000,000","100,000,000","1.000"],
        ["Uniform","2","49,995,116","50,004,884","1.000"],
        ["Uniform","4","24,992,296","25,012,796","1.001"],
        ["Uniform","8","12,484,431","12,513,450","1.001"],
        ["Skewed","1","100,000,000","100,000,000","1.000"],
        ["Skewed","2","49,997,528","50,002,472","1.000"],
        ["Skewed","4","24,991,654","25,006,964","1.000"],
        ["Skewed","8","12,487,846","12,515,662","1.001"],
    ],[1*inch,0.5*inch,1.5*inch,1.5*inch,1*inch])
    p("Oversampling factor 4 was worse here: uniform p=4 became 23,076,323 to 30,760,328 elements and ratio 1.230; skewed was similar at 1.230. A larger sample count did not improve the measured splitters.")
    image("sample_scaling.png")
    story.append(PageBreak())

    p("Case 5: Phases and Communication","Section")
    image("sample_phase.png")
    table([
        ["Uniform p","Local sort","Sampling","Alltoallv","Final sort","Total"],
        ["1","19.381 s","0","0.207 s","3.937 s","23.525 s"],
        ["2","15.801 s","0.089 s","0.230 s","4.111 s","20.141 s"],
        ["4","10.168 s","0.302 s","0.145 s","3.753 s","14.067 s"],
        ["8","9.883 s","2.739 s","0.082 s","3.108 s","13.073 s"],
    ],[1*inch,1.3*inch,1.1*inch,1.1*inch,1.1*inch,1*inch])
    p("<b>Dominant phase.</b> Local sorting dominated all measured ranks. Sampling grew at p=8 but not enough to replace local sorting.")
    p("<b>Communication fit.</b> Ping-pong gave alpha=0.479 microseconds and beta=0.1674 ns/byte. At p=4, the simple model predicts only (p-1)/p of 400 MB moved times beta, plus latency. The result was about 0.053 s before collective overhead. The measured Alltoallv time was about 0.145 s, so the ping-pong fit underestimates real Alltoallv synchronization.")
    p("<b>Serial comparison.</b> Serial qsort took 19.275 s for uniform data and 19.325 s for skewed data. MPI p=4 took about 14.07 s uniform and 13.51 s skewed. Speedup was sub-linear because sorting is O(N log N), while communication moves O(N) bytes and p-rank startup and final sorting remain.")
    p("<b>Skewed result.</b> Skewed p=8 median total was 12.800 s and uniform p=8 was 13.073 s.")
    story.append(PageBreak())

    p("7. Cross-Case Comparison","Section")
    table([
        ["Case","Model","Scaling result","Dominant bottleneck"],
        ["Heat","MPI 2-D","p=2 useful; efficiency low at p=4","Halo communication"],
        ["K-mer","OpenMP","Private 237.7M/s; padded 248.9M/s","Cache and merges"],
        ["Monte Carlo","OpenMP","p=8 speedup 3.64","Hyperthread/thermal/startup"],
        ["PageRank","OpenMP","Dynamic 1.696G edges/s","Memory and irregular work"],
        ["Sample sort","MPI","p=4/p=8 improve 1.4-1.5x","Local sorting and oversubscription"],
    ],[0.9*inch,1*inch,2.1*inch,2.5*inch])
    p("<b>Best scaled case.</b> Monte Carlo scaled best because it is embarrassingly parallel and has almost no shared memory. It reached 3.64 speedup at p=8.")
    p("The k-mer and PageRank cases depend on memory layout and scheduling. Heat and sample sort depend on communication and local computation. Sample sort has strong balance but is limited by O(N log N) sorting and host oversubscription.")
    p("References","Subsection")
    p("1. OpenMP Architecture Review Board, OpenMP Application Programming Interface.")
    p("2. Message Passing Interface Forum, MPI Standard.")
    p("3. Hockney, R. W., The communication cost model T=alpha+beta*m.")
    p("4. Page, Brin, Motwani and Winograd, The PageRank Citation Ranking.")
    p("5. Black, Scholes and Merton, European option closed-form pricing.")
    p("6. Synthetic R-MAT graph generation and 32-bit skewed/uniform data generators used in this repository.")

    def on_page(canvas,doc):
        canvas.saveState()
        canvas.setFont("Helvetica",7)
        canvas.drawString(42,24,"PDC Assignment 3")
        canvas.drawRightString(letter[0]-42,24,"Page %d" % doc.page)
        canvas.restoreState()

    doc.build(story,onFirstPage=on_page,onLaterPages=on_page)


if __name__=="__main__":
    build_report()
    print(PDF)
