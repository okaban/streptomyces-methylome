"""領域構成を 3 つの分母で比べる: ゲノム上の全出現 / 測定可能 (dep>=30) / 修飾あり (lvl>=0.5)。
正準スクリプト nf_r3_recompute_261101.py の関数をそのまま再利用する。"""
import importlib.util, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import chi2_contingency

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("nfr3", os.path.join(HERE, "nf_r3_recompute_261101.py"))
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)

SERIES = list(R.SERIES) + [
    ("GCCGGC", "4mC", [2], "GCCGGC C3"),
    ("GATC",   "6mA", [1], "GATC A2"),
]
CLASSES = R.CLASSES

def base_positions(genome, base):
    """鎖指定 0-based 位置。plus 鎖は base、minus 鎖は相補塩基の位置。"""
    g = np.frombuffer(genome.encode(), dtype="S1")
    comp = {"A": b"T", "C": b"G"}[base]
    pos_p = np.flatnonzero(g == base.encode())
    pos_m = np.flatnonzero(g == comp)
    return (np.concatenate([pos_p, pos_m]),
            np.concatenate([np.ones(len(pos_p), bool), np.zeros(len(pos_m), bool)]))

def run(timepoints):
    genome = R.M.read_reference(R.M.REF)
    tss = pd.read_csv(R.TSSCSV); assert len(tss) == 3570
    cds_mask = np.zeros(len(genome), bool)
    for line in open(R.GFF):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[2] != "CDS": continue
        cds_mask[int(f[3]) - 1:int(f[4])] = True

    rows = []
    for tp in timepoints:
        store = {m: R.pooled(R.TIMEPOINTS[tp], c) for m, c in R.CODE.items()}
        series = SERIES + [("<all A>", "6mA", None, "background all A"),
                           ("<all C>", "4mC", None, "background all C")]
        for motif, mark, offs, label in series:
            t0 = time.time()
            if offs is None:
                pos0, plus = base_positions(genome, "A" if mark == "6mA" else "C")
            else:
                pos0, plus = R.motif_positions(genome, motif, offs)
            s = store[mark]
            idx = pd.MultiIndex.from_arrays([pos0, plus], names=["start", "plus"])
            d = s.reindex(idx)
            dep = d["dep"].to_numpy(); lvl = d["lvl"].to_numpy()
            meas = np.nan_to_num(dep, nan=0) >= 30
            called = meas & (np.nan_to_num(lvl, nan=0) >= 0.5)
            dist = R.signed_tss_distance(pos0 + 1, tss["tss_position"].to_numpy(),
                                         tss["strand"].to_numpy())
            cls = R.classify(pos0, dist, cds_mask)
            for c in CLASSES:
                rows.append(dict(timepoint=tp, mark=mark, series=label, region=c,
                                 n_occurrence=int(cls[c].sum()),
                                 n_measurable=int((cls[c] & meas).sum()),
                                 n_modified=int((cls[c] & called).sum()),
                                 occ_total=int(len(pos0)),
                                 meas_total=int(meas.sum()),
                                 mod_total=int(called.sum())))
            print(f"  {tp} {label:22} occ={len(pos0):>8} meas={meas.sum():>8} mod={called.sum():>7}  {time.time()-t0:.1f}s", flush=True)
    df = pd.DataFrame(rows)
    for col, tot in (("n_occurrence", "occ_total"), ("n_measurable", "meas_total"), ("n_modified", "mod_total")):
        df["pct_" + col.split("_")[1]] = (100 * df[col] / df[tot].replace(0, np.nan)).round(2)
    return df

if __name__ == "__main__":
    tps = sys.argv[1:] or ["T1"]
    df = run(tps)
    out = os.path.join(HERE, "tables", "nf_r3_occurrence_denominator_261102.tsv")
    df.to_csv(out, sep="\t", index=False)
    print("\nwrote", out)
    for tp in tps:
        sub = df[df.timepoint == tp]
        print(f"\n=== {tp}: promoter の割合（%) ===")
        p = sub[sub.region == "promoter"][["series", "pct_occurrence", "pct_measurable", "pct_modified",
                                           "n_occurrence", "n_measurable", "n_modified"]]
        print(p.to_string(index=False))
        print(f"\n=== {tp}: 修飾あり 対 全出現 の分布差 (chi2) ===")
        for lab in sub.series.unique():
            s = sub[sub.series == lab].set_index("region").loc[CLASSES]
            tab = np.vstack([s.n_modified.to_numpy(), s.n_occurrence.to_numpy() - s.n_modified.to_numpy()])
            tab = tab[:, tab.sum(0) > 0]
            chi2, pv, dof, _ = chi2_contingency(tab)
            print(f"  {lab:22} chi2={chi2:10.2f} dof={dof} p={pv:.3g}")
