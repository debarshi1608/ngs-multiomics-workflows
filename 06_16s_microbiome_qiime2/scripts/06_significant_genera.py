# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(E) Step 8 (run after exporting ANCOM-BC results)
import csv
lfc_file = "ancombc-genus-export/lfc_slice.csv"
q_file = "ancombc-genus-export/q_val_slice.csv"
with open(lfc_file, newline="") as f:
    lfc = {row["id"]: float(row["subjectsubject-2"]) for row in csv.DictReader(f)}
with open(q_file, newline="") as f:
    q = {row["id"]: float(row["subjectsubject-2"]) for row in csv.DictReader(f)}
results = []
for feature in lfc:
    if q[feature] < 0.05:
        genus = feature.split(";")[-1].replace("g__", "")
        results.append((genus, lfc[feature], q[feature]))
results.sort(key=lambda x: x[2])
print("Significant genera (q < 0.05)")
print("-----------------------------")
print("Genus\tLFC\tq-value")
for genus, lfc_value, q_value in results:
    direction = "subject-2 enriched" if lfc_value > 0 else "subject-1 enriched"
    print(f"{genus}\t{lfc_value:.4f}\t{q_value:.4g}\t{direction}")
print("\nTotal significant genera:", len(results))
