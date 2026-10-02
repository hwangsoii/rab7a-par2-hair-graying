# Stage 2, checkpoint A helper (2026-10-01): summarise the official NODE
# OEP002321 FASTQ manifest per FASTQ-ID prefix (read-only; no download).
# Output: ../outputs/A0_node_manifest_by_prefix.tsv
import os, re, sys, collections, openpyxl, warnings
warnings.filterwarnings("ignore")
root = os.getcwd()
man = os.path.join(root, "OneDrive_2025-09-11/Data/cell discovery/OEP00002321_Data_1760192653404.xlsx")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "A0_node_manifest_by_prefix.tsv")
mapping = {"ryg035": "F18", "ryg029": "F59", "ryg047": "F62B", "ryg048": "F62W",
           "black-1": "F31B", "white-2": "F31W"}  # official id2FastqID.xlsx (OED00760915)
unit = {"KB": 1e3, "MB": 1e6, "GB": 1e9}
rows = list(openpyxl.load_workbook(man, read_only=True).worksheets[0].iter_rows(values_only=True))
hdr = rows[0]
agg = collections.defaultdict(lambda: {"n": 0, "gb": 0.0, "runs": set(), "flowcell_tags": set(), "style": set()})
for r in rows[1:]:
    d = dict(zip(hdr, r))
    name = d["data_name"]
    if not name or not re.search(r"\.(fq|fastq)\.gz$", name):
        continue
    pre = next((p for p in mapping if name.startswith(p + "-") or name.startswith(p + "_")), None)
    if pre is None:
        sys.exit("unmapped FASTQ prefix: " + name)
    num, u = d["data_size"].split()
    a = agg[pre]
    a["n"] += 1
    a["gb"] += float(num) * unit[u] / 1e9
    a["runs"].add(d["run_id"])
    a["flowcell_tags"].update(re.findall(r"BKDL(\d{4})", name))
    a["style"].add("Illumina_S_L_R" if re.search(r"_S\d+_L\d{3}_[RI]\d_001", name) else "lane_1_2")
with open(out, "w") as fh:
    fh.write("fastq_prefix\tofficial_sample\tn_fastq\ttotal_size_gb\tn_run_ids\tBKDL_tag_YYMM\tname_style\n")
    for p in mapping:
        a = agg[p]
        fh.write(f"{p}\t{mapping[p]}\t{a['n']}\t{a['gb']:.2f}\t{len(a['runs'])}\t"
                 f"{','.join(sorted(a['flowcell_tags']))}\t{','.join(sorted(a['style']))}\n")
print(open(out).read())
