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


def line_graph(name,title,xlabel,ylabel,series,x_ticks=None,ideal=None):
    fig,ax=plt.subplots(figsize=(6.9,3.25))
    for label,x,y,marker in series:
        ax.plot(x,y,marker=marker,label=label)
    if ideal is not None:
        ax.plot(ideal[0],ideal[1],"--",label="Ideal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    ax.grid(True,which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,name),dpi=180)
    plt.close(fig)


def bar_graph(name,title,xlabel,ylabel,labels,values):
    fig,ax=plt.subplots(figsize=(6.9,3.25))
    ax.bar(labels,values)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x",rotation=25)
    ax.grid(True,axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,name),dpi=180)
    plt.close(fig)


def build_assets():
    line_graph("heat_strong.png","Case 1: Strong scaling, 1024 and 4096 grids",
               "MPI ranks","Speedup",
               [("1024 fixed grid",[1,2,4,8],[1.0,1.529,1.548,0.632],"o"),
                ("4096 fixed grid",[1,2,4,8],[1.0,36.380129613/23.021887533,
                                                  36.380129613/19.538150376,
                                                  36.380129613/16.416894491],"s")],
               [1,2,4,8],([1,2,4,8],[1,2,4,8]))
    line_graph("heat_weak.png","Case 1: Weak-scaling efficiency",
               "MPI ranks","Efficiency",
               [("Measured",[1,2,4,8],[1.0,5.740460723/6.743372314,
                                       5.740460723/12.030873209,
                                       5.740460723/29.744415847],"o")],
               [1,2,4,8],([1,2,4,8],[1,1,1,1]))

    km_path=os.path.join(ROOT,"cs2_kmer_openmp","kmer_100m_results.csv")
    km=read_csv(km_path)
    k8={r["mode"]:float(r["throughput"]) for r in km if r["k"]=="8"}
    bar_graph("kmer.png","Case 2: k=8 throughput at 100 million characters",
              "Implementation","Million k-mers/s",
              ["Serial","Race","Atomic","Private","Padded"],
              [k8["serial"],k8["race"],k8["atomic"],k8["private"],k8["padded"]])
    private_times=[0.186190201,0.075566101,0.045952500,0.076444700]
    serial_kmer=0.172408302
    line_graph("kmer_scaling.png","Case 2: Private histogram strong scaling",
               "OpenMP threads","Speedup",
               [("Private", [1,2,4,8],[serial_kmer/private_times[i] for i in range(4)],"o")],
               [1,2,4,8],([1,2,4,8],[1,2,4,8]))

    serial_m=statistics.median([10.947842,10.852722,9.757764])
    medians=[serial_m,
             statistics.median([12.146964,10.070304,10.602978]),
             statistics.median([6.467711,6.842922,6.002942]),
             statistics.median([3.740699,2.616511,4.327151]),
             statistics.median([3.604589,2.978711,2.628556])]
    threads=[1,2,4,8]
    speed=[serial_m/medians[i] for i in range(1,5)]
    line_graph("mc.png","Case 3: Monte Carlo fixed-size strong scaling",
               "OpenMP threads","Speedup",
               [("Measured",threads,speed,"o")],
               threads,([1,2,4,8],[1,2,4,8]))

    degree=read_csv(os.path.join(ROOT,"cs4_pagerank_openmp","degree_distribution.csv"))
    d=[int(r["degree"]) for r in degree]
    c=[int(r["count"]) for r in degree]
    fig,ax=plt.subplots(figsize=(6.9,3.25))
    ax.loglog(d,c,".",markersize=3)
    ax.set_xlabel("In-degree")
    ax.set_ylabel("Number of vertices")
    ax.set_title("Case 4: R-MAT in-degree distribution")
    ax.grid(True,which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"pagerank_degree.png"),dpi=180)
    plt.close(fig)

    bar_graph("pagerank_schedule.png","Case 4: 3-run median schedule throughput",
              "Schedule and chunk","Million edges/s",
              ["Static64","Static256","Dynamic64","Dynamic256","Guided64","Guided256"],
              [809.44,788.98,1695.96,1726.18,777.32,712.37])

    bar_graph("pagerank_pull_push.png","Case 4: Serial, pull and push",
              "Implementation","Million edges/s",
              ["Serial","Pull, 4 threads","Push, 4 threads"],
              [425.21,726.35,157.31])

    ranks=[1,2,4,8]
    su=19.2749
    ss=19.3254
    uniform=[23.5246,20.1410,14.0665,13.0726]
    skewed=[28.0142,20.6089,13.5135,12.8004]
    line_graph("sample_scaling.png","Case 5: 100-million-key strong scaling",
               "MPI ranks","Speedup",
               [("Uniform",[1,2,4,8],[su/x for x in uniform],"o"),
                ("Skewed",[1,2,4,8],[ss/x for x in skewed],"s")],
               ranks,([1,2,4,8],[1,2,4,8]))

    labels=["1","2","4","8"]
    local=[19.3807,15.8008,10.1680,9.8830]
    sampling=[0.0,0.0886,0.3020,2.7392]
    alltoallv=[0.2071,0.2296,0.1453,0.0822]
    final=[3.9368,4.1105,3.7528,3.1079]
    fig,ax=plt.subplots(figsize=(6.9,3.25))
    ax.bar(labels,local,label="Local sort")
    ax.bar(labels,sampling,bottom=local,label="Sampling")
    bottom=[local[i]+sampling[i] for i in range(4)]
    ax.bar(labels,alltoallv,bottom=bottom,label="Alltoallv")
    bottom=[bottom[i]+alltoallv[i] for i in range(4)]
    ax.bar(labels,final,bottom=bottom,label="Final sort")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("Seconds")
    ax.set_title("Case 5: Uniform-data phase breakdown")
    ax.grid(True,axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"sample_phase.png"),dpi=180)
    plt.close(fig)

    bytes=[8,64,512,4096,32768,131072,524288,1048576,4194304,16777216]
    us=[0.479,0.264,0.515,3.775,7.924,17.616,53.662,151.441,725.160,2799.819]
    fig,ax=plt.subplots(figsize=(6.9,3.25))
    ax.loglog(bytes,us,"o-")
    ax.set_xlabel("Message size (bytes)")
    ax.set_ylabel("One-way time (microseconds)")
    ax.set_title("Case 5: MPI ping-pong result")
    ax.grid(True,which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS,"pingpong.png"),dpi=180)
    plt.close(fig)


def table(story,data,widths):
    t=Table(data,colWidths=widths,hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONTNAME",(0,0),(-1,-1),"Helvetica"),
        ("FONTSIZE",(0,0),(-1,-1),7.2),
        ("LEADING",(0,0),(-1,-1),8.5),
        ("GRID",(0,0),(-1,-1),0.25,colors.grey),
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#d9e7f5")),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ]))
    story.append(t)


def build_report():
    build_assets()
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name="CoverTitle",parent=styles["Title"],fontSize=24,leading=29,alignment=TA_CENTER,spaceAfter=18))
    styles.add(ParagraphStyle(name="CoverSub",parent=styles["Normal"],fontSize=11,leading=15,alignment=TA_CENTER,textColor=colors.HexColor("#444444")))
    styles.add(ParagraphStyle(name="SectionNew",parent=styles["Heading1"],fontSize=16,leading=19,textColor=colors.HexColor("#1f4e79"),spaceBefore=3,spaceAfter=7))
    styles.add(ParagraphStyle(name="SubNew",parent=styles["Heading2"],fontSize=11,leading=13,textColor=colors.HexColor("#1f4e79"),spaceBefore=4,spaceAfter=3))
    styles.add(ParagraphStyle(name="BodyNew",parent=styles["BodyText"],fontSize=8.7,leading=12))

    doc=SimpleDocTemplate(PDF,pagesize=letter,rightMargin=40,leftMargin=40,topMargin=40,bottomMargin=40,
                          title="PDC Assignment 3 Report",author="sfsgr")
    story=[]

    def p(text,style="BodyNew"):
        story.append(Paragraph(text,styles[style]))
    def img(name,width=6.9,height=3.25):
        story.append(Spacer(1,0.05*inch))
        story.append(Image(os.path.join(ASSETS,name),width=width,height=height))

    # Page 1
    story.append(Spacer(1,1.15*inch))
    p("PDC Assignment 3","CoverTitle")
    p("Parallel and Distributed Computing Case Studies","CoverSub")
    story.append(Spacer(1,0.25*inch))
    p("Group XX","CoverSub")
    p("Team member: Tayyab - GitHub: sfsgr","CoverSub")
    p("November 2026","CoverSub")
    p("Case studies: heat diffusion, k-mer counting, Monte Carlo, PageRank and sample sort","CoverSub")
    story.append(PageBreak())

    # Page 2
    p("1. Introduction","SectionNew")
    p("This report describes five parallel implementations and answers each critical-analysis question. OpenMP is used for k-mer counting, Monte Carlo and PageRank. MPI is used for heat diffusion and sample sort.")
    p("<b>Machine.</b> Intel Core i5-8250U with 4 physical cores and 8 logical cores, 7.7 GiB RAM, Microsoft WSL 2 Ubuntu.")
    p("<b>Toolchain.</b> GCC 15.2.0, GNU Make, OpenMP, OpenMPI 5 through mpicc/mpirun, and Python for plotting.")
    p("<b>Timing method.</b> Each graph uses the median of at least 3 runs where possible. The idle machine is not guaranteed to be thermally stable. MPI p=8 is oversubscribed because WSL exposes 8 logical CPUs but only 4 physical cores.")
    p("<b>Repository and contributions.</b> Source code, Makefiles, measurement files and this report are in PDC_Asst_3.")
    table(story,[
        ["Case","Main implementation","Key evidence"],
        ["1","1-D and 2-D MPI, blocking/non-blocking halo","Timing CSV and correctness sums"],
        ["2","Serial, race, atomic, private, padded OpenMP","Histogram totals and throughput"],
        ["3","Independent-thread Monte Carlo streams","Strong-scaling and fixed-time results"],
        ["4","CSR pull/push PageRank and schedules","Degree and timing output"],
        ["5","MPI regular sample sort and ping-pong","Balance, phases and Hockney fit"],
    ],[0.8*inch,2.4*inch,3.4*inch])
    p("<b>Common checks.</b> Integer results and element counts are exact. Floating-point PageRank and Monte Carlo values use a tolerance.")
    p("<b>Contribution evidence.</b> Git commits separate the case-study source, data, report and plot changes.")
    story.append(PageBreak())

    # Page 3
    p("2. Case 1: Heat Diffusion - Implementation","SectionNew")
    p("The stencil is a time-dependent 5-point update. It is parallel because each block uses values from the previous time step.")
    p("<b>Chosen model.</b> MPI. The first layout is 1-D row blocks. The improved layout is 2-D Cartesian blocks created with MPI_Cart_create.")
    p("<b>Communication.</b> A blocking rank uses MPI_Sendrecv. A non-blocking rank posts Irecv/Isend requests, updates interior cells, waits, then updates boundary-adjacent cells.")
    p("<b>1-D cost.</b> A rank computes N^2/p cells and exchanges 2N halo values per step.")
    p("<b>2-D cost.</b> With px and py ranks, a rank computes N^2/p and exchanges 2N/px + 2N/py halo values. For a square layout, this is about 2N/sqrt(p). This is less than 2N for p>1.")
    p("<b>Convergence.</b> Every K steps, local maximum change is sent with MPI_Allreduce.")
    p("<b>Strong scaling at N=1024.</b>")
    table(story,[
        ["p","2-D blocking time","Speedup","Efficiency"],
        ["1","5.743 s","1.00","100%"],
        ["2","4.049 s","1.529","76.5%"],
        ["4","3.999 s","1.548","38.7%"],
        ["8","9.793 s","0.632","7.9%"],
    ],[0.7*inch,2.1*inch,1.5*inch,1.5*inch])
    img("heat_strong.png")
    p("The p=1 speedup is 1.00 because it is the required serial baseline reference. Efficiency first drops below 70% at p=4.")
    story.append(PageBreak())

    # Page 4
    p("Case 1: Strong/Weak Scaling and Communication","SectionNew")
    p("<b>Strong scaling.</b> At p=4, the total time is only 1.548 times faster than serial. Communication and synchronization dominate. At p=8, oversubscription adds to the loss.")
    p("<b>Large-grid scaling.</b> A 4096 grid gives a median serial time of 36.380 s, compared with 23.022 s at p=2, 19.538 s at p=4 and 16.417 s at p=8. Median speedups are 1.58, 1.86 and 2.22.")
    table(story,[
        ["N","p","Median time","Speedup","Efficiency"],
        ["1024","1","5.743 s","1.00","100%"],
        ["1024","2","4.049 s","1.53","76.5%"],
        ["1024","4","3.999 s","1.55","38.7%"],
        ["1024","8","9.793 s","0.63","7.9%"],
        ["4096","1","36.380 s","1.00","100%"],
        ["4096","2","23.022 s","1.58","79.0%"],
        ["4096","4","19.538 s","1.86","46.6%"],
        ["4096","8","16.417 s","2.22","27.7%"],
    ],[0.8*inch,0.5*inch,1.5*inch,1*inch,1.2*inch])
    p("<b>Weak scaling.</b>")
    img("heat_weak.png")
    table(story,[
        ["p","N","Time","Weak efficiency"],
        ["1","1024","5.740 s","100%"],
        ["2","1448","6.743 s","85.1%"],
        ["4","2048","12.031 s","47.7%"],
        ["8","2896","29.744 s","19.3%"],
    ],[0.8*inch,1*inch,1.5*inch,2*inch])
    story.append(PageBreak())

    # Page 5
    p("Case 1: Critical Answers","SectionNew")
    p("<b>Surface to volume.</b> 1-D has N^2/p computed cells and 2N halo values per rank. 2-D has N^2/p computed cells and 2N/px+2N/py halo values. 2-D is better for large p because halo volume is divided in two directions.")
    p("<b>Overlap.</b> For 2-D blocking/non-blocking exchange, p=2 used 0.558/0.803 s, p=4 used 1.965/1.451 s, and p=8 used 8.631/3.054 s. The p=8 result hides about (8.631-3.054)/8.631=64.6% of blocking exchange time.")
    p("Little overlap was visible at p=2. Shared-memory transfer and small iteration work limit progress. MPI can also spend much of the wait in post or synchronization work.")
    p("<b>Why efficiency drops.</b> Tp is estimated by Tcomp/p+Tcomm. At p=4 and p=8 the halo, startup and Allreduce costs are no longer small compared with local update work.")
    p("<b>K=1 versus K=50.</b> At p=4, K=1 took 5.002 s and K=50 took 3.534 s. K=1 performs one Allreduce every step. K=50 performs 10 Allreduces for 500 steps.")
    p("<b>Correctness.</b> Serial and MPI sums matched within 1e-9. The 4096 run also preserved the same sum at p=2,4,8.")
    story.append(PageBreak())

    # Page 6
    p("3. Case 2: K-mer Counting","SectionNew")
    p("DNA bases are converted to 2-bit values and each k-mer is encoded as an integer. The serial baseline slides one window at a time.")
    p("<b>Race.</b> The parallel loop updated a shared histogram without protection. For k=8 and 100 million characters, the race result had 99,946,568 counts instead of 99,999,993.")
    p("<b>Atomic.</b> Each update is protected with omp atomic.")
    p("<b>Private.</b> Every thread has a histogram copy. Copies are merged after the parallel region.")
    p("<b>Padded private.</b> Thread-local histogram rows are padded to 64-byte cache lines.")
    p("<b>Chunk boundaries.</b> Whole k-mer windows are assigned to threads. Each chunk begins at encode_kmer(begin,k). Later windows use a rolling 2-bit index.")
    table(story,[
        ["k","Implementation","4 threads","Correct"],
        ["8","Atomic","61.966 M/s","Yes"],
        ["8","Private","237.719 M/s","Yes"],
        ["8","Padded","248.884 M/s","Yes"],
        ["12","Atomic","51.713 M/s","Yes"],
        ["12","Private","24.718 M/s","Yes"],
    ],[0.5*inch,1.5*inch,1.7*inch,1.2*inch])
    img("kmer.png")
    story.append(PageBreak())

    # Page 7
    p("Case 2: Critical Answers","SectionNew")
    p("<b>Race precisely.</b> Two threads can read the same counter value, add 1, and write the same result. The second write discards the other increment.")
    p("<b>Atomic versus private.</b>")
    table(story,[
        ["k","Private memory, 1 thread","Extra memory, 4 threads","Fastest measured"],
        ["8","512 KiB","2 MiB","Padded private"],
        ["12","128 MiB","512 MiB","Atomic"],
    ],[0.5*inch,2.1*inch,2.1*inch,1.5*inch])
    p("For k=8, private histograms fit cache better and atomic operations are avoided. For k=12, 16,777,216 counters per thread consume memory and private merge becomes expensive.")
    p("<b>False sharing.</b> For k=1 and 10 million characters with 2 threads, private mode took 0.1079 s and padded took 0.0716 s. Throughput rose from 92.698M/s to 139.690M/s, about 50.7% faster.")
    p("<b>Strong scaling.</b> The private k=8 10-million-character median data gives p=1 0.926 speedup, p=2 2.28, p=4 3.75 and p=8 2.25 because p=8 oversubscribes 4 physical cores.")
    img("kmer_scaling.png")
    p("<b>Memory or compute.</b> Each k-mer reads sequence bytes and updates a counter. A simple estimate is 9 bytes read plus histogram write traffic. Padded k=8 reached 248.884 million k-mers/s. Scaling flattens when threads compete for cache, memory bandwidth and merge cost.")
    story.append(PageBreak())

    # Page 8
    p("4. Case 3: Monte Carlo Option Pricing","SectionNew")
    p("A European call option is estimated from geometric Brownian motion. The parameters are S0=100, K=100, r=5%, sigma=20%, T=1 year. The Black-Scholes price is 10.4505835722.")
    p("<b>OpenMP model.</b> Each thread gets an independent stream and a private sum/square sum. The threads combine totals after the loop.")
    p("<b>Fixed-size strong scaling.</b>")
    table(story,[
        ["Threads","Median time","Speedup","Throughput","Serial fraction estimate"],
        ["Serial","10.8527 s","1.00","9.21M/s","-"],
        ["1","10.6030 s","1.02","9.43M/s","-"],
        ["2","6.4677 s","1.68","15.46M/s","19%"],
        ["4","3.7407 s","2.90","26.73M/s","13%"],
        ["8","2.9787 s","3.64","33.57M/s","17%"],
    ],[1*inch,1.4*inch,1*inch,1.3*inch,1.6*inch])
    img("mc.png")
    p("The p=8 value is not ideal because 8 threads use SMT, but it is close enough to show an embarrassingly parallel pattern.")
    story.append(PageBreak())

    # Page 9
    p("Case 3: Critical Answers","SectionNew")
    p("<b>RNG matters.</b> Independent per-thread streams produced 10.468089 at 4 threads for 1 million paths, with standard error 0.014734. Reusing the same seed gave 10.512458 and standard error 0.014780, but a larger price error 0.061875 against Black-Scholes.")
    p("The same seed creates correlated streams. It does not destroy output every time, but it makes the statistical error misleading.")
    p("<b>Fixed time.</b> The p=1 median result used 90.783744 million samples in 10.000135 s. The p=8 median result used 361.336832 million samples in 10.007727 s.")
    p("The p=1 median standard error was 0.001545, while p=8 was 0.000774. The sample count ratio is 3.98, and sqrt(3.98)=1.996, matching the standard-error reduction.")
    p("A speedup of 3.64 is fast time scaling. In fixed-time statistical terms, p=8 reduces SE by about half.")
    p("<b>Why it scales.</b> Almost all work is independent. Simple Amdahl estimates give 13-17% serial fraction. Remaining losses are reduction, startup, turbo/thermal changes and SMT.")
    p("<b>Reproducibility.</b> The same seed and thread count gave the same answer. Changing p changed the grouping of floating-point sums, so the last digits can differ. The result is reproducible for fixed p and should be compared with a tolerance.")
    story.append(PageBreak())

    # Page 10
    p("5. Case 4: Parallel PageRank","SectionNew")
    p("The R-MAT graph has 65,536 vertices and 1,335,037 edges. CSR is used for pull updates. Transposed CSR is used for incoming neighbour lists.")
    p("<b>Dangling nodes.</b> The sum of ranks for out-degree zero vertices is redistributed equally as part of the base rank.")
    p("<b>Pull.</b> One thread owns destination vertices and reads r_cur/u/outdeg(u). It has no write conflicts.")
    p("<b>Push.</b> One thread owns sources and adds to destination ranks with omp atomic.")
    table(story,[
        ["Version","Median time","Million edges/s","Rank sum"],
        ["Serial","0.0440 s","425.21","1.0000000000"],
        ["Pull, 4 threads","0.0257 s","726.35","1.0000000000"],
        ["Push, 4 threads","0.1188 s","157.31","1.0000000000"],
    ],[2*inch,1.4*inch,1.7*inch,1.4*inch])
    p("Pull is 1.71 times faster than serial and 4.62 times faster than push.")
    img("pagerank_degree.png")
    p("The degree plot is power-law-like. Most vertices have low in-degree, while a few have very high in-degree.")
    story.append(PageBreak())

    # Page 11
    p("Case 4: Critical Answers","SectionNew")
    p("<b>Irregularity.</b> In one 3-run median, static imbalance R was 1.066. The fastest thread had 554,057 edges in its range and the least-loaded thread had 20,369. Equal vertex counts produced 27.20 times more work for the busiest range.")
    p("<b>Scheduling.</b>")
    img("pagerank_schedule.png")
    table(story,[
        ["Schedule","Chunk","3-run median time","Million edges/s"],
        ["static","64","0.0247 s","809.44"],
        ["static","256","0.0254 s","788.98"],
        ["dynamic","64","0.0118 s","1695.96"],
        ["dynamic","256","0.0116 s","1726.18"],
        ["guided","64","0.0258 s","777.32"],
        ["guided","256","0.0281 s","712.37"],
    ],[1.2*inch,0.8*inch,2*inch,1.7*inch])
    p("Dynamic scheduling is best because a thread can take another chunk after finishing. Small chunks reduce imbalance but add scheduling work. Dynamic 64 is the measured best.")
    p("<b>Pull versus push.</b> Push spends time on atomic destination writes. Pull has random reads but no atomics.")
    img("pagerank_pull_push.png")
    p("<b>Memory behavior.</b> Pull moves about 20-32 bytes of rank, index and degree data per edge. Dynamic pull reached about 1.696 billion edges/s. This is close to memory-bandwidth pressure, but WSL cache and SMT prevent a precise RAM-bandwidth comparison.")
    p("<b>Correctness.</b> Serial, pull and push had rank sum 1.0000000000 and residual 4.445e-8.")
    story.append(PageBreak())

    # Page 12
    p("6. Case 5: Distributed Sample Sort","SectionNew")
    p("Each rank generates a local shard, sorts locally, takes regular samples, gathers splitters and partitions. MPI_Alltoall exchanges counts, then MPI_Alltoallv exchanges keys.")
    p("<b>Verification.</b> Each local array is sorted, max_i<=min_i+1 across neighbours, and MPI_Allreduce preserves the total key count.")
    p("<b>Uniform input, 100 million keys.</b>")
    table(story,[
        ["p","Serial/MPI total","Speedup","Min","Max","Max/mean"],
        ["1","23.525 s","0.82","100,000,000","100,000,000","1.000"],
        ["2","20.141 s","0.96","49,995,116","50,004,884","1.000"],
        ["4","14.067 s","1.37","24,992,296","25,012,796","1.001"],
        ["8","13.073 s","1.48","12,484,431","12,513,450","1.001"],
    ],[0.5*inch,1.7*inch,1*inch,1.5*inch,1.5*inch,1*inch])
    p("<b>Skewed input.</b>")
    table(story,[
        ["p","Time","Speedup","Min","Max","Max/mean"],
        ["1","28.014 s","0.69","100,000,000","100,000,000","1.000"],
        ["2","20.609 s","0.94","49,997,528","50,002,472","1.000"],
        ["4","13.514 s","1.43","24,991,654","25,006,964","1.000"],
        ["8","12.800 s","1.51","12,487,846","12,515,662","1.001"],
    ],[0.5*inch,1.4*inch,1*inch,1.5*inch,1.5*inch,1*inch])
    p("The p=1 time includes MPI start and the local two-sort structure. It is slower than direct serial qsort.")
    img("sample_scaling.png")
    story.append(PageBreak())

    # Page 13
    p("Case 5: Critical Answers","SectionNew")
    p("<b>Oversampling.</b> With factor 1, p=4 imbalance was 1.001 uniform and 1.000 skewed. Factor 4 produced 1.230 for both distributions. Oversampling gathered more samples, but it also chose a worse splitter pattern in these runs.")
    p("Oversampling normally helps by giving more candidates and rejecting edge splitters. It does not guarantee balance, especially with duplicates or a large skewed key cluster.")
    p("<b>Phase table, uniform.</b>")
    img("sample_phase.png")
    table(story,[
        ["p","Local sort","Sampling","Alltoallv","Final sort","Total"],
        ["1","19.381 s","0","0.207 s","3.937 s","23.525 s"],
        ["2","15.801 s","0.089 s","0.230 s","4.111 s","20.141 s"],
        ["4","10.168 s","0.302 s","0.145 s","3.753 s","14.067 s"],
        ["8","9.883 s","2.739 s","0.082 s","3.108 s","13.073 s"],
    ],[0.5*inch,1.3*inch,1*inch,1.1*inch,1.1*inch,1*inch])
    p("Local sort dominates at p=1 to p=8. Sampling becomes larger at p=8 but remains smaller than local sorting.")
    p("<b>Communication cost.</b>")
    img("pingpong.png")
    p("Ping-pong fitted alpha=0.479 microseconds and beta=0.1674 ns/byte, with 5.56 GB/s asymptotic bandwidth.")
    p("At p=4, (p-1)/p times 400 MB equals 300 MB. The simple Tcomm is 300 MB x 0.1674 ns/byte plus 3alpha, about 0.053 s. The measured Alltoallv time was 0.145 s. Collective synchronization and rank imbalance add the difference.")
    p("<b>Comparison to serial sort.</b> Serial qsort took 19.275 s uniform and 19.325 s skewed. MPI p=4 reached 14.067 s and 13.514 s. Speedup is sub-linear because qsort has O(N log N) work, while Alltoallv moves O(N) bytes and final sort repeats work.")
    story.append(PageBreak())

    # Page 14
    p("7. Cross-Case Comparison","SectionNew")
    table(story,[
        ["Case","Model","Best measured result","Dominant bottleneck"],
        ["Heat","MPI","4096 p=8, 2.22x speedup","Halo transfer and oversubscription"],
        ["K-mer","OpenMP","248.9M k-mers/s","Cache and private merge"],
        ["Monte Carlo","OpenMP","p=8, 3.64x speedup","SMT and startup/reduction"],
        ["PageRank","OpenMP","Dynamic 1.696G edges/s","Irregular memory work"],
        ["Sample sort","MPI","p=8, about 1.5x speedup","Local sort and sampling"],
    ],[1*inch,0.9*inch,2.2*inch,2.5*inch])
    p("<b>Best case.</b> Monte Carlo scaled best because its work is independent and only private sums are reduced.")
    p("<b>Weak points.</b> Heat is communication-limited. K-mer is histogram-cache limited. PageRank is irregular and memory-limited. Sample sort is sorting-limited and oversubscribed at p=8.")
    p("<b>References.</b>")
    p("1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.")
    p("2. Message Passing Interface Forum. MPI Standard.")
    p("3. R. W. Hockney. Communication model T=alpha+beta*m.")
    p("4. L. Page, S. Brin, R. Motwani and T. Winograd. The PageRank Citation Ranking.")
    p("5. F. Black, M. Scholes and R. Merton. European call option pricing.")
    p("6. GCC, OpenMPI and GNU Make documentation.")
    p("7. Repository measurement CSV/TXT files for all reported medians.")

    def on_page(canvas,doc):
        canvas.saveState()
        canvas.setFont("Helvetica",7)
        canvas.drawString(40,22,"PDC Assignment 3")
        canvas.drawRightString(letter[0]-40,22,"Page %d" % doc.page)
        canvas.restoreState()

    doc.build(story,onFirstPage=on_page,onLaterPages=on_page)


if __name__=="__main__":
    build_report()
    print(PDF)
