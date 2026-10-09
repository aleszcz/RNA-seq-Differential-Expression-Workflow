# RNA-seq Differential Expression Workflow

Reference, QC'd, and brought up to current (2026) best practice. Corrected your HISAT2 syntax, added the missing index-build and SAM→BAM steps, swapped in a couple of tools that have mostly replaced the older ones in practice, and noted where a pipeline manager (Nextflow/nf-core) saves you from hand-running all of this.

## Pipeline at a glance

```
SRA/GEO accession
      │
      ▼
[1] Download raw reads  ──────────────  FASTQ (gzipped, paired-end)
      │
      ▼
[2] Read QC + trimming  ──────────────  fastp  (QC report + cleaned FASTQ)
      │
      ▼
[3] Reference download  ──────────────  genome FASTA + GTF  (Ensembl)
      │
      ▼
[4] Build HISAT2 index  ──────────────  genome.*.ht2  (+ splice-site info from GTF)
      │
      ▼
[5] Align reads         ──────────────  SAM → sorted, indexed BAM
      │
      ▼
[6] Alignment QC        ──────────────  samtools flagstat / RSeQC
      │
      ▼
[7] Assign reads to genes ────────────  featureCounts (or HTSeq-count)
      │
      ▼
[8] Differential expression ──────────  DESeq2 (R)
      │
      ▼
[9] Downstream           ─────────────  shrunk LFCs, PCA, volcano, GO/pathway enrichment
```

Your original outline had the right backbone (FASTQ → align → count → DESeq2). The fixes below are: a broken HISAT2 command, two missing steps (index build, SAM→BAM conversion), and a few tool swaps that reflect what people actually run now instead of 2015-era defaults.

---

## What changed vs. your draft, and why

| Step | What you had | What's improved | Why |
|---|---|---|---|
| Download | "SRA explorer / GEO files" | ENA direct FASTQ links where possible, `prefetch` + `fasterq-dump` otherwise | ENA mirrors SRA and hosts FASTQ directly — no `.sra`→FASTQ conversion step needed, which is the single biggest time-saver in step 1 |
| QC/trim | *(missing)* | `fastp` | Does adapter trimming + QC in one pass, much faster than the old FastQC+Trimmomatic combo, and writes an HTML/JSON report |
| HISAT2 command | `hisat2_[options]* -x<bt2-idx> {-1<m1>-2<m2>1.4cr} [-s<sam>` | corrected flag syntax below | the flags were run together and missing the index-build step entirely |
| Splice awareness | *(not used)* | `--known-splicesite-infile` built from the GTF | meaningfully improves junction alignment accuracy over vanilla HISAT2, and takes one extra script you already have (bundled with HISAT2) |
| SAM handling | *(missing)* | `samtools sort` + `samtools index` into BAM | HTSeq/featureCounts and all downstream QC tools need sorted, indexed BAM, not raw SAM |
| Counting | HTSeq | **featureCounts** (Subread), HTSeq as fallback | faster, handles paired-end and multi-mapping reads more robustly, and is more actively maintained; HTSeq is still fine if your lab's existing scripts expect its output format |
| LFC estimates | plain DESeq2 output | `lfcShrink()` with `apeglm` | apeglm shrinkage is now the DESeq2-recommended default for ranking/visualizing fold changes — raw LFCs from low-count genes are noisy |
| Whole pipeline | manual bash, step by step | optional: **nf-core/rnaseq** (Nextflow) | if you'll run this on more than a couple of sample sets, nf-core/rnaseq wraps steps 1–7 in one containerized, resumable pipeline with QC built in — worth knowing about even if you run the manual version below for learning purposes |

---

## Step 1 — Get the raw reads

Your notes list three ways to do this — here's each one corrected/filled in, in order of how little setup they need.

**Option A — SRA Explorer** ([sra-explorer.info](https://sra-explorer.info/)). Paste a study/SRA identifier (GSE, SRP, PRJNA, or SRR), it resolves the runs and gives you, per sample: a direct download link to raw FASTQ, a ready-made bash script that downloads all of them with their SRR identifiers as filenames, and plain copyable links. The links point at the ENA FTP mirror, which hosts FASTQ directly — no `.sra`→FASTQ conversion step needed, which is what makes this option the fastest for a handful of samples.

On a cluster: copy the generated bash script (or individual `wget`/`curl` links), paste it into your home or scratch directory on the login/head node (`user@header:~/GEO_data/` or similar), and run it there — no need for an interactive compute session just to download, since it's only network I/O.
```bash
# what SRA Explorer's bash script does per sample, for reference:
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/<ACC_PREFIX>/<ACCESSION>/<ACCESSION>_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/<ACC_PREFIX>/<ACCESSION>/<ACCESSION>_2.fastq.gz
```
I can fetch these directly if you give me the accession — no need to open SRA Explorer by hand first.

**Option B — SRA Toolkit** (the "official" NCBI route). Install from the [sra-tools GitHub](https://github.com/ncbi/sra-tools) (releases page has prebuilt binaries — building from source isn't necessary). On a shared cluster this one *does* need an interactive session (`srun --pty bash` / `qsub -I` or equivalent) since `prefetch`/`fasterq-dump` do real CPU/disk work, not just a login-node download.
```bash
# prefetch grabs the compressed .sra object first (resumable, checksum-verified)
prefetch SRR_ACCESSION

# fasterq-dump converts it to FASTQ — this is the slow/CPU step
# --split-files: separate R1/R2 for paired-end
# --progress: show progress bar
# -e 8: threads, parallelizes the dump itself
fasterq-dump SRR_ACCESSION --split-files --progress -e 8

# fasterq-dump doesn't gzip its output, so compress after:
pigz -p 8 SRR_ACCESSION_1.fastq SRR_ACCESSION_2.fastq   # pigz = parallel gzip, much faster than plain gzip
```
Name output files by their SRR identifier (the default) so they stay traceable back to the GEO/SRA accession table.

**Option C — SRA Run Selector / "SRA downloader" (step-through GUI).** NCBI's [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) lets you browse a study's runs in a table, tick the ones you want, and either get an accession list to feed into `prefetch`, or use the browser-based cloud download. Good for cherry-picking a subset of samples from a large study when you don't want everything SRA Explorer would hand you at once.

**Which to use:** SRA Explorer for speed when you want everything in a GSE; SRA Toolkit when SRA Explorer's ENA mirror is missing a run (happens occasionally for newer/embargoed submissions) or your institution requires the "official" NCBI tool; Run Selector when you need to pick specific runs out of a big study rather than all of them.

## Step 2 — QC and trim

```bash
fastp \
  -i sample_1.fastq.gz -I sample_2.fastq.gz \
  -o sample_1.trim.fastq.gz -O sample_2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --thread 8 \
  --json sample.fastp.json --html sample.fastp.html
```

Run `multiqc .` afterward once you've got several samples' `fastp` reports to get one combined QC summary.

## Step 3 — Reference genome and annotation (Ensembl)

```bash
# Genome FASTA (use primary_assembly or toplevel, not individual chromosomes)
wget https://ftp.ensembl.org/pub/release-XXX/fasta/<species>/dna/<Species>.<assembly>.dna.primary_assembly.fa.gz

# GTF annotation — this is what tells you which gene each read belongs to
wget https://ftp.ensembl.org/pub/release-XXX/gtf/<species>/<Species>.<assembly>.XXX.gtf.gz

gunzip *.gz
```
Replace `release-XXX`, `<species>`, `<Species>`, `<assembly>` with the current Ensembl release and your organism (e.g. `release-114`, `mus_musculus`, `Mus_musculus`, `GRCm39`). I can pull the exact current URLs for your organism if you tell me which one.

## Step 4 — Build the HISAT2 index (this step was missing)

```bash
# Extract splice sites and exons from the GTF — improves alignment accuracy
hisat2_extract_splice_sites.py genome.gtf > splicesites.txt
hisat2_extract_exons.py genome.gtf > exons.txt

# Build the index, incorporating known splice sites
hisat2-build -p 8 \
  --ss splicesites.txt --exon exons.txt \
  genome.fa genome_index
```

## Step 5 — Align (your corrected command)

Your draft: `hisat2_[options]* -x<bt2-idx> {-1<m1>-2<m2>1.4cr} [-s<sam>` — flags need spaces and the right dashes, and you need to pipe straight to sorted BAM rather than leaving a loose SAM file around.

```bash
hisat2 -p 8 \
  -x genome_index \
  --known-splicesite-infile splicesites.txt \
  -1 sample_1.trim.fastq.gz -2 sample_2.trim.fastq.gz \
  --summary-file sample.align_summary.txt \
  | samtools sort -@ 8 -o sample.sorted.bam -

samtools index sample.sorted.bam
```

(`-x` = index prefix, `-1`/`-2` = paired mates, piping into `samtools sort` avoids ever writing an uncompressed SAM to disk.)

## Step 6 — Alignment QC

```bash
samtools flagstat sample.sorted.bam > sample.flagstat.txt
# optional, more detail:
# read_distribution.py -i sample.sorted.bam -r genome.bed   (RSeQC)
```

## Step 7 — Assign reads to genes

```bash
featureCounts -p --countReadPairs -T 8 \
  -a genome.gtf \
  -o counts.txt \
  sample1.sorted.bam sample2.sorted.bam sample3.sorted.bam ...
```
`-p --countReadPairs` tells it the BAMs are paired-end. If you specifically need HTSeq-compatible output instead:
```bash
htseq-count -f bam -r pos -s no -i gene_id \
  sample.sorted.bam genome.gtf > sample.htseq.txt
```

## Step 8 — Differential expression (DESeq2, R)

```r
library(DESeq2)

counts <- read.table("counts.txt", header = TRUE, row.names = 1, skip = 1)
counts <- counts[, 6:ncol(counts)]  # drop featureCounts' annotation columns

colData <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "control", "treated", "treated", "treated")
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# filter near-zero-count genes before running
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treated", "control"))

# current recommended practice: shrink fold changes for ranking/plotting
res_shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")

res_ordered <- res_shrunk[order(res_shrunk$padj), ]
write.csv(as.data.frame(res_ordered), "DE_results.csv")
```

## Step 9 — Downstream (optional but worth doing)

Same three sub-steps as before, written out in full rather than one-liners — these all run right after Step 8, using the `dds` and `res_shrunk` objects already in your R session.

### 9a. PCA / sample QC

Checks that your samples cluster by condition rather than by batch, before you trust any DE result.

```r
# variance-stabilizing transform — like log2, but handles low counts better
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")

# if you have a batch/replicate variable, color by it too to spot batch effects:
# plotPCA(vsd, intgroup = c("condition", "batch"))
```
If samples don't separate cleanly by `condition` on PC1/PC2, that's worth investigating (mislabeled sample, batch effect, outlier) before interpreting the DE results — don't skip this and go straight to the gene list.

### 9b. Volcano plot

Plots every tested gene by effect size (shrunk log2FoldChange) vs. significance (`padj`), same convention as the DESeq2 docs.

```r
library(ggplot2)

res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)
res_df$sig <- with(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "log2 fold change (shrunk)", y = "-log10 adjusted p-value",
       title = "Treated vs control") +
  theme_minimal()
```
Using the **shrunk** `log2FoldChange` here (not raw) matters — it's the same `apeglm`-shrunk values from Step 8, so low-count genes with wild raw fold changes don't dominate the plot.

### 9c. Pathway / GO enrichment

Takes your significant gene list and asks what biological processes/pathways are over-represented in it.

```r
library(clusterProfiler)
library(org.Hs.eg.db)   # swap for your organism's annotation package, e.g. org.Mm.eg.db for mouse

sig_genes <- rownames(res_shrunk)[which(res_shrunk$padj < 0.05)]

# over-representation test on your significant gene set
ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",     # match this to whatever IDs are in your counts (Ensembl gene IDs from featureCounts)
  ont           = "BP",          # Biological Process; also "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot(ego, showCategory = 15)
write.csv(as.data.frame(ego), "GO_enrichment.csv")

# ranked version (GSEA) if you want to use the whole gene list, not just the significant cutoff:
ranked <- res_shrunk$log2FoldChange
names(ranked) <- rownames(res_shrunk)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

gsea <- gseGO(geneList = ranked, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
```
`enrichGO` tests only your significant genes against a background; `gseGO` instead ranks *every* gene by fold change and tests whether a pathway's genes skew toward one end — more sensitive when you have a lot of small, consistent changes rather than a few huge ones. Worth running both if you're not sure which failure mode you're in.

---

## If you'll run this more than once: nf-core/rnaseq

Everything from step 1 (raw reads) through step 7 (gene counts), plus QC at every stage, is wrapped into one pipeline here: https://nf-co.re/rnaseq. It runs in Docker/Singularity/Conda, resumes from where it failed, and produces a MultiQC report automatically. Worth switching to once you're past understanding each step by hand — the manual version above is what it's doing internally.

## On downloads

I can fetch the Ensembl genome/GTF files or pull FASTQ for a given SRA/GEO accession directly in this conversation if you give me the accession numbers or species — no need to go through SRA Explorer by hand first.
