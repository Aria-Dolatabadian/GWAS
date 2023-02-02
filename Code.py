import pandas as pd
df = pd.read_csv('gwass.csv')
print(df)
import matplotlib.pyplot as plt
import geneview as gv
ax = gv.manhattanplot(data=df)
plt.show()
ax = gv.manhattanplot(data=df, xticklabel_kws={"rotation": "vertical"})
plt.show()
ax = gv.manhattanplot(data=df,
                   suggestiveline=None,  # Turn off suggestiveline
                   genomewideline=None,  # Turn off genomewideline
                   xticklabel_kws={"rotation": "vertical"})
plt.show()
# plot only results of chromosome 8.
gv.manhattanplot(data=df, CHR="chr8", xlabel="Chromosome 8")
plt.show()
ax = gv.manhattanplot(data=df,
                   sign_marker_p=1e-6,  # highline the significant SNP with ``sign_marker_color`` color.
                   is_annotate_topsnp=True,  # annotate the top SNP
                   xticklabel_kws={"rotation": "vertical"})
plt.show()
import matplotlib.pyplot as plt
import geneview as gv
# common parameters for plotting
plt_params = {
    "font.sans-serif": "Arial",
    "legend.fontsize": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14
}
plt.rcParams.update(plt_params)
# Create a manhattan plot
f, ax = plt.subplots(figsize=(12, 4), facecolor="w", edgecolor="k")
xtick = set(["chr" + i for i in list(map(str, range(1, 10))) + ["11", "13", "15", "18", "21", "X"]])
_ = gv.manhattanplot(data=df,
                     marker=".",
                     sign_marker_p=1e-6,  # Genome wide significant p-value
                     sign_marker_color="r",
                     snp="ID",  # The column name of annotation information for top SNPs.

                     title="Test",
                     xtick_label_set=xtick,

                     xlabel="Chromosome",
                     ylabel=r"$-log_{10}{(P)}$",

                     sign_line_cols=["#D62728", "#2CA02C"],
                     hline_kws={"linestyle": "--", "lw": 1.3},

                     is_annotate_topsnp=True,
                     ld_block_size=50000,  # 50000 bp
                     text_kws={"fontsize": 12,
                               "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                     ax=ax)
plt.show()
ax = gv.qqplot(data=df["P"])
plt.show()
#Show a better QQ plot
f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
_ = gv.qqplot(data=df["P"],
              marker="o",
              title="Test",
              xlabel=r"Expected $-log_{10}{(P)}$",
              ylabel=r"Observed $-log_{10}{(P)}$",
              ax=ax)
plt.show()
#Admixture plot
import matplotlib.pyplot as plt
from geneview.utils import load_dataset
from geneview import admixtureplot
f, ax = plt.subplots(1, 1, figsize=(14, 2), facecolor="w", constrained_layout=True, dpi=300)
admixtureplot(data=load_dataset("admixture_output.Q"),
              population_info=load_dataset("admixture_population.info"),
              ylabel_kws={"rotation": 45, "ha": "right"},
              ax=ax)
plt.show()
