import gzip, csv, os

# 10 significant mediation genes from Phase 5
target_genes = ["GDF11", "CLU", "HNF1A", "MYL6B", "HDGFRP2", "ZNF655",
                "FAM195A", "NDUFAF2", "RP5-1115A15.1"]

# Read Phase 5 mediation results to get eQTLGen direction
print("=== Reading Phase 5 mediation results ===")
with open("results/tables/06_mediation_results.csv") as f:
    med = list(csv.DictReader(f))
sig_med = [r for r in med if r["sig_mediation"] == "TRUE"]
print(f"  {len(sig_med)} significant mediation paths")

# Get gene -> step1 beta (smoking -> expression in eQTLGen blood)
eqtlgen_dir = {}
for r in sig_med:
    g = r["gene"]
    if g not in eqtlgen_dir:
        eqtlgen_dir[g] = []
    eqtlgen_dir[g].append({
        "outcome": r["outcome"],
        "beta_step1": float(r["beta_step1"]),
        "p_step1": float(r["p_step1"]),
        "beta_step2": float(r["beta_step2"]),
        "p_step2": float(r["p_step2"]),
        "mediation_effect": float(r["mediation_effect"]),
        "prop_mediated": float(r["prop_mediated"])
    })

print(f"  Genes: {list(eqtlgen_dir.keys())}")

# Tissues to check
tissues = ["Lung", "Whole_Blood"]
gtex_dir = "data/gtex_lung/GTEx_Analysis_v8_eQTL"

results = []

for tissue in tissues:
    print(f"\n=== Processing {tissue} ===")

    # Step 1: Build gene_name -> gene_id mapping from egenes file
    egenes_file = f"{gtex_dir}/{tissue}.v8.egenes.txt.gz"
    gene_map = {}  # gene_name -> (gene_id, qval, slope, pval_nominal)
    with gzip.open(egenes_file, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gname = row["gene_name"]
            if gname in target_genes:
                qval = float(row["qval"]) if row["qval"] != "NA" else 1.0
                gene_map[gname] = {
                    "gene_id": row["gene_id"],
                    "qval": qval,
                    "slope": float(row["slope"]),
                    "pval_nominal": float(row["pval_nominal"]),
                    "is_egenie": qval < 0.05,
                    "lead_variant": row.get("variant_id", ""),
                    "rsid": row.get("rs_id_dbSNP151_GRCh38p7", "")
                }

    print(f"  Found {len(gene_map)}/{len(target_genes)} genes in {tissue} egenes")
    for g, info in gene_map.items():
        print(f"    {g} ({info['gene_id']}): qval={info['qval']:.2e}, "
              f"egenie={'YES' if info['is_egenie'] else 'no'}, "
              f"beta={info['slope']:.4f}, p={info['pval_nominal']:.2e}")

    # Step 2: Get all significant eQTLs for target genes
    sig_file = f"{gtex_dir}/{tissue}.v8.signif_variant_gene_pairs.txt.gz"
    target_ensembl_ids = set(info["gene_id"] for info in gene_map.values())

    gene_eqtls = {}  # gene_name -> list of eQTLs
    with gzip.open(sig_file, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["gene_id"] in target_ensembl_ids:
                # Find which gene this is
                for gname, info in gene_map.items():
                    if info["gene_id"] == row["gene_id"]:
                        if gname not in gene_eqtls:
                            gene_eqtls[gname] = []
                        gene_eqtls[gname].append({
                            "variant_id": row["variant_id"],
                            "maf": float(row["maf"]),
                            "pval": float(row["pval_nominal"]),
                            "slope": float(row["slope"]),
                            "slope_se": float(row["slope_se"])
                        })
                        break

    # Step 3: For each gene, get lead eQTL and compare direction
    for gname in target_genes:
        if gname not in gene_map:
            results.append({
                "gene": gname, "tissue": tissue,
                "found": False, "is_egenie": False,
                "qval": "", "n_eqtls": 0, "lead_variant": "", "lead_p": "",
                "lead_beta": "", "lead_maf": "",
                "eqtlgen_beta_step1": "", "gtex_direction": "",
                "mediation_outcomes": ""
            })
            continue

        info = gene_map[gname]
        eqtls = gene_eqtls.get(gname, [])
        # Sort by p-value to get lead eQTL
        eqtls.sort(key=lambda x: x["pval"])
        lead = eqtls[0] if eqtls else None

        # Get eQTLGen direction for this gene
        eqtlgen_info = eqtlgen_dir.get(gname, [{}])[0]
        eqtlgen_beta = eqtlgen_info.get("beta_step1", None)

        # Direction: GTEx slope sign vs eQTLGen beta sign
        # Note: eQTLGen beta is for smoking->expression, GTEx slope is for variant->expression
        # They're not directly comparable in sign, but we can note the GTEx direction
        dir_consistent = ""
        if lead and eqtlgen_beta is not None:
            # Both are effect sizes on expression, but from different exposures
            # Just record the GTEx direction
            dir_consistent = "positive" if lead["slope"] > 0 else "negative"

        results.append({
            "gene": gname,
            "tissue": tissue,
            "found": True,
            "is_egenie": info["is_egenie"],
            "qval": f"{info['qval']:.2e}",
            "n_eqtls": len(eqtls),
            "lead_variant": lead["variant_id"] if lead else "",
            "lead_p": f"{lead['pval']:.2e}" if lead else "",
            "lead_beta": f"{lead['slope']:.4f}" if lead else "",
            "lead_maf": f"{lead['maf']:.4f}" if lead else "",
            "eqtlgen_beta_step1": f"{eqtlgen_beta:.4f}" if eqtlgen_beta is not None else "",
            "gtex_direction": dir_consistent,
            "mediation_outcomes": "; ".join(r["outcome"] for r in eqtlgen_dir.get(gname, []))
        })

# Save results
out_file = "results/tables/10_gtex_multitissue_validation.csv"
with open(out_file, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "gene", "tissue", "found", "is_egenie", "qval",
        "n_eqtls", "lead_variant", "lead_p", "lead_beta", "lead_maf",
        "eqtlgen_beta_step1", "gtex_direction", "mediation_outcomes"
    ])
    w.writeheader()
    for r in results:
        w.writerow(r)

print(f"\n\n=== Summary ===")
print(f"Results saved to {out_file}")
print(f"\n--- Lung tissue ---")
lung_results = [r for r in results if r["tissue"] == "Lung"]
for r in lung_results:
    if r["found"]:
        status = "eGene" if r["is_egenie"] else "not eGene"
        print(f"  {r['gene']:15s}: {status}, qval={r.get('qval','NA')}, "
              f"n_eqtls={r['n_eqtls']}, lead_beta={r['lead_beta']}, lead_p={r['lead_p']}")
    else:
        print(f"  {r['gene']:15s}: NOT FOUND in GTEx Lung")

print(f"\n--- Whole Blood tissue ---")
blood_results = [r for r in results if r["tissue"] == "Whole_Blood"]
for r in blood_results:
    if r["found"]:
        status = "eGene" if r["is_egenie"] else "not eGene"
        print(f"  {r['gene']:15s}: {status}, qval={r.get('qval','NA')}, "
              f"n_eqtls={r['n_eqtls']}, lead_beta={r['lead_beta']}, lead_p={r['lead_p']}")
    else:
        print(f"  {r['gene']:15s}: NOT FOUND in GTEx Whole Blood")

# Cross-tissue comparison
print(f"\n--- Cross-tissue eGene status ---")
for g in target_genes:
    lung = next((r for r in lung_results if r["gene"] == g), {})
    blood = next((r for r in blood_results if r["gene"] == g), {})
    lung_s = "eGene" if lung.get("is_egenie") else ("not-eGene" if lung.get("found") else "NA")
    blood_s = "eGene" if blood.get("is_egenie") else ("not-eGene" if blood.get("found") else "NA")
    print(f"  {g:15s}: Lung={lung_s:12s}  Blood={blood_s:12s}")
