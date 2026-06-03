import os
from collections import Counter, defaultdict
from datetime import datetime
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import json
from moepy import lowess

class Fastq:
    def read_fastq(self, path, amp=None, filter_seq=False):
        with open(path, "r") as f:
            counts = Counter(
                line.rstrip().upper()
                for _, line, _, _ in zip(f, f, f, f)
            )
        counts = dict(counts)

class General:
    def rt_seq(self, seq):
        rt = ""
        dic = {
            "A": "T",
            "G": "C",
            "T": "A",
            "C": "G",
            "N": "N",
            "a": "t",
            "g": "c",
            "t": "a",
            "c": "g",
        }

        for i in seq[::-1]:
            rt += dic[i]

        return rt

class DataParsing:
    def dfdic(self, df):
        dflst = self.divide_dataframe(df)
        dic = {}

        for eachdf in dflst:
            sb = eachdf.iloc[0, 0]
            dic[sb] = eachdf

        return dic

    def divide_dataframe(self, df):
        bclst = self._barcode_indexing(df)

        t1 = datetime.now()
        n = 0
        lst = []

        while n < len(bclst) - 1:
            tup0 = bclst[n]
            tup1 = bclst[n + 1]

            idx = tup0[1]
            idx1 = tup1[1]

            eachdf = df[idx:idx1]
            lst.append(eachdf)

            if n == len(bclst) - 2:
                lastdf = df[idx1:]
                lst.append(lastdf)

            n += 1

        t2 = datetime.now()
        print("Make DF list:", t2 - t1)

        return lst

    def _barcode_indexing(self, df):
        first_column = list(df.columns)[0]
        df = df.sort_values(by=first_column)

        t1 = datetime.now()

        n = 1
        sb = list(df.iloc[:, 0])
        lst = [(sb[0], 0)]

        while n < df.shape[0]:
            if sb[n] != sb[n - 1]:
                lst.append((sb[n], n))
            n += 1

        t2 = datetime.now()
        print("Barcode_Indexing:", t2 - t1)

        return lst

    def dflst(self, df):
        bclst = self._barcode_indexing(df)

        t1 = datetime.now()
        n = 0
        lst = []

        while n < len(bclst) - 1:
            tup0 = bclst[n]
            tup1 = bclst[n + 1]

            idx = tup0[1]
            idx1 = tup1[1]

            eachdf = df[idx:idx1]
            lst.append(eachdf)

            if n == len(bclst) - 2:
                lastdf = df[idx1:]
                lst.append(lastdf)

            n += 1

        t2 = datetime.now()
        print("Make DF list:", t2 - t1)

        return lst

    def divide_list(self, lst, num):
        chunk_size = len(lst) // num
        final = []

        for i in range(num):
            if i == num - 1:
                eachlst = lst[chunk_size * i :]
            else:
                eachlst = lst[chunk_size * i : chunk_size * (i + 1)]

            final.append(eachlst)

        return final

class Codon:
    def __init__(self):
        self.codon = {
            "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
            "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
            "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
            "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
            "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
            "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
            "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
            "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
            "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
            "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
            "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
            "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
            "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
            "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
            "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
            "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
        }

        self.codon_short = {
            "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
            "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
            "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
            "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
            "STOP": "*",
        }

        self.codon_abb = {
            "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
            "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
            "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
        }

        self.nucleotide = {
            "A": ["A"],
            "C": ["C"],
            "G": ["G"],
            "T": ["T"],
            "R": ["A", "G"],
            "Y": ["C", "T"],
            "S": ["G", "C"],
            "W": ["A", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "B": ["C", "G", "T"],
            "D": ["A", "G", "T"],
            "H": ["A", "C", "T"],
            "V": ["A", "C", "G"],
            "N": ["A", "T", "C", "G"],
        }

def run_fastp_workflow(root_folder, output_folder=None):
    root_path = Path(root_folder)
    output_path = Path(output_folder) if output_folder else None

    sample_dirs = [
        d
        for d in root_path.rglob("*")
        if d.is_dir()
        and "sample_" in d.name.lower()
    ]

    if not sample_dirs:
        return

    summary_stats = []
    files_to_cleanup = []

    for folder in sample_dirs:
        sample_id = (
            folder.name
            .replace("sample_", "")
            .replace("Sample_", "")
        )

        r1_files = list(folder.glob("*_1.fq.gz"))

        for r1_file in r1_files:
            r2_file = folder / r1_file.name.replace(
                "_1.fq.gz",
                "_2.fq.gz",
            )

            if not r2_file.exists():
                continue

            parent_dir = (
                output_path
                if output_path
                else folder.parent
            )

            merged_out = parent_dir / f"{sample_id}.fastq"

            unmerged_r1 = (
                parent_dir
                / f"{sample_id}_unmerged_1.fq.gz"
            )

            unmerged_r2 = (
                parent_dir
                / f"{sample_id}_unmerged_2.fq.gz"
            )

            html_report = (
                parent_dir
                / f"{sample_id}_report.html"
            )

            json_report = (
                parent_dir
                / f"{sample_id}_report.json"
            )

            command = [
                fastp_path,
                "-i", str(r1_file),
                "-I", str(r2_file),
                "--merge",
                "--merged_out", str(merged_out),
                "--out1", str(unmerged_r1),
                "--out2", str(unmerged_r2),
                "--disable_adapter_trimming",
                "--overlap_len_require", "12",
                "--overlap_diff_limit", "0",
                "--overlap_diff_percent_limit", "0",
                "--qualified_quality_phred", "20",
                "--html", str(html_report),
                "--json", str(json_report),
                "--thread", "16",
            ]

            try:
                subprocess.run(
                    command,
                    check=True,
                    capture_output=True,
                    text=True,
                )

                if json_report.exists():
                    with open(json_report, "r") as f:
                        data = json.load(f)

                    filt = data["filtering_result"]

                    total = (
                        filt["passed_filter_reads"]
                        + filt["low_quality_reads"]
                        + filt["too_many_N_reads"]
                    )

                    summary_stats.append(
                        {
                            "id": sample_id,
                            "total": total,
                            "passed": filt["passed_filter_reads"],
                            "low_qual": filt["low_quality_reads"],
                            "many_n": filt["too_many_N_reads"],
                        }
                    )

                files_to_cleanup.extend(
                    [
                        unmerged_r1,
                        unmerged_r2,
                        html_report,
                        json_report,
                    ]
                )

            except subprocess.CalledProcessError:
                pass

    print_summary_table(summary_stats)

    for temp_file in files_to_cleanup:
        if temp_file.exists():
            temp_file.unlink()

def load_reference_sequences(json_file, exon):
    with open(json_file) as f:
        reference_dict = json.load(f)

    cds = reference_dict.get(f"{exon}_cds")
    amp = reference_dict.get(f"{exon}_amp")

    if cds is not None:
        cds = cds.upper()

    if amp is not None:
        amp = amp.upper()

    return cds, amp

def calculate_log10_or(
    df,
    d0_count_col,
    wt_count_col,
    d0_wt_raw,
    wt_wt_raw,
):
    return df.apply(
        lambda r: np.log10(
            ((r[d0_count_col] + 1) / (d0_wt_raw + 1))
            /
            ((r[wt_count_col] + 1) / (wt_wt_raw + 1))
        ),
        axis=1,
    )

def calculate_fisher_p(
    df,
    d0_count_col,
    wt_count_col,
    d0_wt_raw,
    wt_wt_raw,
):
    return df.apply(
        lambda r: fisher_exact(
            [
                [int(r[d0_count_col]), int(r[wt_count_col])],
                [int(d0_wt_raw), int(wt_wt_raw)],
            ]
        )[1],
        axis=1,
    )

def calculate_fdr(p_values):
    _, fdr, _, _ = multipletests(
        p_values,
        method="fdr_bh",
    )
    return fdr

def add_rpm_columns(df, sample_sums):
    for sample_name, total_sum in sample_sums.items():
        df[f"{sample_name}_RPM"] = (
            df[f"{sample_name}_Count"] / total_sum * 1e6
        )

    return df

def filter_variants(
    df,
    d0_rpm_col,
    rpm_cutoff=0.5,
    or_cutoff=3,
    fdr_cutoff=0.05,
):
    rpm_filter = df[d0_rpm_col] > rpm_cutoff
    or_filter = df["log10_OR"] > np.log10(or_cutoff)
    fdr_filter = df["FDR"] < fdr_cutoff

    combined_filter = (
        rpm_filter
        & or_filter
        & fdr_filter
    )

    filtered_df = df.loc[combined_filter].copy()

    filtered_df = filtered_df[
        ~filtered_df["ID"].isin(["WT", "FALSE"])
    ]

    return filtered_df

def calculate_lowess(fit_position, fit_value, predict_position, frac):
    fit_df = (
        pd.DataFrame({
            "position": fit_position,
            "value": fit_value
        })
        .sort_values("position")
    )

    model = lowess.Lowess()

    model.fit(
        fit_df["position"].to_numpy(),
        fit_df["value"].to_numpy(),
        frac=frac,
    )

    return model.predict(np.asarray(predict_position))

def calculate_weighted_lfc(group, lfc_col, rpm_col):
    lfc = pd.to_numeric(group[lfc_col], errors="coerce")
    rpm = pd.to_numeric(group[rpm_col], errors="coerce")

    valid = lfc.notna() & rpm.notna()

    lfc = lfc[valid]
    rpm = rpm[valid]

    if rpm.sum() == 0:
        return np.nan

    return (lfc * rpm).sum() / rpm.sum()

def classify_variants(
    df,
    r1_col,
    r2_col,
    type_col="type",
):
    syn = df[df[type_col] == "synonymous"].copy()

    if len(syn) < 5:
        return ["Intermediate"] * len(df)

    thr_res1 = np.percentile(
        syn[r1_col].dropna(),
        99.87,
    )
    thr_res2 = np.percentile(
        syn[r2_col].dropna(),
        99.87,
    )

    thr_sens1 = np.percentile(
        syn[r1_col].dropna(),
        97.72,
    )
    thr_sens2 = np.percentile(
        syn[r2_col].dropna(),
        97.72,
    )

    classifications = []

    for _, row in df.iterrows():
        v1 = row[r1_col]
        v2 = row[r2_col]

        if pd.isna(v1) or pd.isna(v2):
            classifications.append("Intermediate")

        elif v1 >= thr_res1 and v2 >= thr_res2:
            classifications.append("Resistant")

        elif v1 <= thr_sens1 and v2 <= thr_sens2:
            classifications.append("Sensitive")

        else:
            classifications.append("Intermediate")

    return classifications
