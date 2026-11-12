# Trajectory timing summary per cluster (Fig 1 / S3): early, gap and late segment shares and peak age.
import os
import numpy as np
import pandas as pd

ROOT = os.environ["APA_ROOT"]
D2 = f"{ROOT}/data/02_traj"

sup = pd.read_csv(f"{D2}/01_support.tsv", sep="\t")
ga = sup.grid_age.values
am = sup.age_model.values
ok = sup.is_supported.astype(bool).values
segs = np.split(np.arange(len(ok)), np.where(np.diff(ok.astype(int)) != 0)[0] + 1)
segs = [s for s in segs if ok[s[0]]]

tl = pd.read_csv(f"{D2}/traj_logit.csv", index_col=0)
sc = pd.read_csv(f"{D2}/02_sites_clustered.tsv", sep="\t", usecols=["site_id", "cluster_id", "delta_logit"],
                 low_memory=False).dropna(subset=["cluster_id"]).set_index("site_id")
T = tl.loc[sc.index].values
e, l = segs[0], segs[1]
early = T[:, e[-1]] - T[:, e[0]]
gap = T[:, l[0]] - T[:, e[-1]]
late = T[:, l[-1]] - T[:, l[0]]
net = T[:, l[-1]] - T[:, e[0]]
tot = np.abs(early) + np.abs(gap) + np.abs(late)
df = pd.DataFrame({
    "cluster": sc.cluster_id.map({1.0: "C1", 2.0: "C2"}).values,
    "early": np.abs(early) / tot,
    "gap": np.abs(gap) / tot,
    "late": np.abs(late) / tot,
    "net": net,
})
d = np.diff(T, axis=1) / np.diff(am)
mid = (ga[:-1] + ga[1:]) / 2
okm = ok[:-1] & ok[1:]
d[:, ~okm] = 0
dirn = np.sign(net)[:, None]
peak = mid[np.argmax(d * dirn, axis=1)]
df["peak_age"] = peak
df["peak_prenatal"] = peak < 0
summ = df.groupby("cluster").agg(
    n=("net", "size"),
    net_sign_pos=("net", lambda x: (x > 0).mean()),
    early_share=("early", "median"),
    gap_share=("gap", "median"),
    late_share=("late", "median"),
    peak_age_median=("peak_age", "median"),
    frac_peak_prenatal=("peak_prenatal", "mean"),
    frac_peak_late=("peak_age", lambda x: (x >= 11).mean()),
)
summ.to_csv(f"{D2}/s3_timing_sum.tsv", sep="\t")
