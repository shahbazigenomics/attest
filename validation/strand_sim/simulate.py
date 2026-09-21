"""Simulated libraries for calibrating the strandedness check.

A random 2 Mb genome with 400 single-exon genes on both strands (15% silent),
60 of them in antisense-overlapping pairs (the case that separates stranded from
unstranded counting). Reads are 75 bp single-end, drawn from genes by a
log-normal expression profile, in one of three protocols:

  reverse    - dUTP / TruSeq Stranded: the read is the reverse complement of the RNA
  forward    - e.g. Ligation / SMART: the read has the RNA's sequence
  unstranded - either, 50:50

An optional share of reads comes from outside every gene (intronic/intergenic
signal, rRNA-like background), to show what that does to the numbers.
Real STAR and featureCounts are then run on the output (run.sh).
"""
import random, gzip, os, sys, zlib

random.seed(20260921)
G = 2_000_000
comp = str.maketrans("ACGT", "TGCA")
genome = "".join(random.choice("ACGT") for _ in range(G))

# genes: (id, start0, end0, strand); some in antisense-overlapping pairs
genes, pos = [], 5_000
while len(genes) < 400 and pos < G - 20_000:
    L = random.randint(800, 4000)
    s = random.choice("+-")
    genes.append((f"gene{len(genes)+1:04d}", pos, pos + L, s))
    if len(genes) <= 60 and len(genes) % 2 == 1:            # antisense partner, overlapping
        o = random.randint(200, L // 2)
        genes.append((f"gene{len(genes)+1:04d}", pos + o, pos + o + L, "-" if s == "+" else "+"))
        pos += o + L + random.randint(2_000, 6_000)
    else:
        pos += L + random.randint(2_000, 6_000)

os.makedirs("ref", exist_ok=True)
with open("ref/genome.fa", "w") as f:
    f.write(">chrS\n")
    for i in range(0, G, 60):
        f.write(genome[i:i+60] + "\n")
with open("ref/genes.gtf", "w") as f:
    for gid, a, b, s in genes:
        attr = f'gene_id "{gid}"; transcript_id "{gid}.t1";'
        f.write(f"chrS\tsim\texon\t{a+1}\t{b}\t.\t{s}\t.\t{attr}\n")

# 15% of genes silent, as in any real tissue - a matrix where every gene has
# reads looks pre-filtered to attest_completeness(), correctly
expr = {g[0]: (0.0 if random.random() < 0.15 else random.lognormvariate(0, 1.5)) for g in genes}
tot = sum(expr.values())
covered = sorted((a, b) for _, a, b, _ in genes)

def outside():
    while True:
        p = random.randint(0, G - 80)
        if not any(a - 80 < p < b for a, b in covered):
            return p

def reads(protocol, n, background, sample_seed):
    rng = random.Random(sample_seed)
    ids = [g for g in genes]
    w = [expr[g[0]] * rng.lognormvariate(0, 0.3) for g in ids]  # sample-to-sample variation
    if sum(w) == 0: w = [1.0] * len(ids)
    RL = 75
    for i in range(n):
        if rng.random() < background:
            p = outside(); seq = genome[p:p+RL]
            if rng.random() < 0.5: seq = seq.translate(comp)[::-1]
        else:
            gid, a, b, s = rng.choices(ids, weights=w)[0]
            p = rng.randint(a, b - RL)
            rna = genome[p:p+RL] if s == "+" else genome[p:p+RL].translate(comp)[::-1]
            if protocol == "forward":   seq = rna
            elif protocol == "reverse": seq = rna.translate(comp)[::-1]
            else:                       seq = rna if rng.random() < 0.5 else rna.translate(comp)[::-1]
        yield f"@r{i}\n{seq}\n+\n{'I'*RL}\n"

if __name__ == "__main__":
    os.makedirs("fastq", exist_ok=True)
    plan = [("reverse", 0.0), ("forward", 0.0), ("unstranded", 0.0), ("reverse", 0.6)]
    for protocol, bg in plan:
        for k in range(3):                                          # 3 samples each
            name = f"{protocol}{'_bg' if bg else ''}_s{k+1}"
            with gzip.open(f"fastq/{name}.fq.gz", "wt") as f:
                n_reads = random.randint(110_000, 190_000)          # libraries differ in depth
                for r in reads(protocol, n_reads, bg, zlib.crc32(name.encode())):   # hash() is salted per process
                    f.write(r)
            print(name, flush=True)
