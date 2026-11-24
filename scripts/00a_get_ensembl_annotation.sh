#!/usr/bin/env bash
# ==============================================================================
# 00a_get_ensembl_annotation.sh — Ensembl release 115 gene table
# Downloads the GTF that includes ALT haplotypes/patches and extracts one row
# per gene: ensembl_gene_id, chromosome_name, external_gene_name, gene_biotype.
# Output: data/annotation/ensembl115_genes.tsv  (used by 01a and 01b)
# ==============================================================================
set -euo pipefail

mkdir -p data/annotation
GTF=data/annotation/ensembl115.gtf.gz
OUT=data/annotation/ensembl115_genes.tsv

if [ ! -s "$GTF" ]; then
  curl -fL -o "$GTF" \
    https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz
fi

zcat "$GTF" | awk -F'\t' 'BEGIN{OFS="\t"; print "ensembl_gene_id","chromosome_name","external_gene_name","gene_biotype"}
  $3=="gene"{a=$9; id=a; sub(/.*gene_id "/,"",id); sub(/".*/,"",id);
  nm="NA"; if(a ~ /gene_name "/){nm=a; sub(/.*gene_name "/,"",nm); sub(/".*/,"",nm)};
  bt=a; sub(/.*gene_biotype "/,"",bt); sub(/".*/,"",bt); print id,$1,nm,bt}' > "$OUT"

echo "Genes written: $(($(wc -l < "$OUT") - 1))"
grep -E "ENSG00000204580|ENSG00000230456" "$OUT"
