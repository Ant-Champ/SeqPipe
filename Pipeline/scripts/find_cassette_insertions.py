#!/usr/bin/env python3
"""find_cassette_insertions -- antibiotic-cassette insertion sites from a sorted,
indexed BAM mapped against genome contigs + cassette contigs.

Cassettes usually REPLACE a gene; reads spanning a genome/cassette boundary
(supplementary / split alignments) give the junctions. Cassettes can share
homologous ends (e.g. Keio FRT/priming sites), so the cassette part of a
junction read is often MAPQ 0: the junction POSITION is still used, but the
cassette IDENTITY is only claimed from uniquely mapped (MAPQ >= min-mapq)
cassette evidence on the same read or, for short reads, on its mate. Reads
without such evidence are 'ambiguous' (candidates = the junction segment's
contig + its secondary / XA alternatives); they are never assigned silently.

Steps: junctions per read -> clusters by (contig, genome side) -> left/right
junctions paired into insertion sites (span <= max-span) -> per-cassette counts.
A read gives a junction when
  * split read: query-adjacent genome + cassette segments, each >= min-segment
    aligned bp (supplementary alignments; ONT and short reads), or
  * soft clip (short reads, needs --fasta): a uniquely mapped genome read whose
    soft clip (>= --min-clip bp) continues into a cassette sequence. Most short
    reads over a junction are only clipped (e.g. 48M27S, no supplementary
    alignment), so this is the main evidence at moderate depth. Identity: the
    cassettes containing the clipped bases, else the mate if it maps uniquely.
Read pairs with one mate unique on the genome and the other on a cassette are
counted per genome kb and printed to the log ("mate links"). When the junction
reads cannot name the cassette (shared ends), the pairs ending within
--mate-window bp of the junction (genome side) whose other mate is unique on a
cassette name it (counted in support, not in n_unique). A cassette with reads but
no junction and no mate links is not next to the reference genome (e.g. on a
plasmid, or the wrong reference).

Hotspots (--out-hotspots): every read on a cassette is followed back to the
genome -- its mate (unique on the genome), its soft clip (the clipped bases
located once, exactly, in the genome; any read on a cassette, primary or
supplementary), and the junction reads above. Evidence closer than
--hotspot-gap bp is merged into one hotspot; hotspots with >= --hotspot-min-reads
reads are reported (snapshot_regions.py draws each), even when no precise
junction pair could be called.

Coordinates are 0-based half-open. A junction 'pos' is the boundary on the
genome: genome_side 'left' = genome lies at < pos (cassette continues right),
'right' = genome lies at >= pos. For a paired site start = left pos, end = right
pos, so length = end - start = genome bases replaced. Outputs are always
written (header-only when empty).

--out-junctions columns ('#' header line):
  #contig  start  end  genome_side  cassette  cassette_end  orientation  n_reads
  n_unique  n_ambiguous  support  junction_id
  cassette = name | mixed(A:n,B:m) | ambiguous(A|B); support = unique reads per
  cassette ("A:12,B:0"; mates near the junction when the junction reads had none); orientation '+' = genome and cassette on same strand.
--out-insertions columns:
  #contig  start  end  length  status(paired|left_only|right_only)  cassette
  orientation  n_reads_left  n_reads_right  left_junction  right_junction
  (first three columns stay contig/start/end so a features column can be appended)
--out-hotspots columns:
  #contig  start  end  length  status(hotspot)  cassette  n_reads  n_mates
  n_genome_clip  n_cassette_clip  n_split  junctions (top positions, pos:reads)
--out-counts columns:
  #cassette  length  reads_unique  reads_ambiguous  breadth_unique_pct
  mean_depth_unique   (primary alignments with MAPQ >= / < min-mapq; breadth and
  depth from MAPQ >= min-mapq reads only)

Python >= 3.8; needs pysam.
"""
import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REF_OPS = (0, 2, 3, 7, 8)       # M D N = X consume reference
QRY_OPS = (0, 1, 4, 5, 7, 8)    # M I S H = X span the original read


def _clips(r):
    """(ref_start, ref_end, qstart, qend) of an alignment, query interval in the
    ORIGINAL read orientation (cigar is in reference orientation: flip if reverse)."""
    ct = r.cigartuples
    lc = ct[0][1] if ct[0][0] in (4, 5) else 0
    rc = ct[-1][1] if ct[-1][0] in (4, 5) else 0
    qlen = sum(n for op, n in ct if op in QRY_OPS)
    qs, qe = lc, qlen - rc
    if r.is_reverse:
        qs, qe = qlen - qe, qlen - qs
    return r.reference_start, r.reference_end, qs, qe


def _xa_contigs(r):
    """Contigs listed in a bwa XA tag (chr,+pos,CIGAR,NM;...)."""
    if not r.has_tag("XA"):
        return set()
    return {x.split(",")[0] for x in r.get_tag("XA").split(";") if x}


def _read_key(r):
    n = 1 if r.flag & 0x40 else 2 if r.flag & 0x80 else 0
    return (r.query_name, n)


def _cluster(events, tol):
    """Group dicts (key 'pos') whose sorted positions chain within tol."""
    if not events:
        return []
    events = sorted(events, key=lambda e: e["pos"])
    clusters, cur = [], [events[0]]
    for e in events[1:]:
        if e["pos"] - cur[-1]["pos"] <= tol:
            cur.append(e)
        else:
            clusters.append(cur)
            cur = [e]
    clusters.append(cur)
    return clusters


def _call_set(call):
    """Cassette names behind a call string + its kind (named|mixed|amb)."""
    m = re.match(r"(mixed|ambiguous)\((.*)\)$", call)
    if not m:
        return {call}, "named"
    if m.group(1) == "mixed":
        return {x.split(":")[0] for x in m.group(2).split(",")}, "mixed"
    return set(m.group(2).split("|")) - {"."}, "amb"


def _compatible(c1, c2):
    s1, k1 = _call_set(c1)
    s2, k2 = _call_set(c2)
    if k1 == k2 == "named":
        return s1 == s2
    return bool(s1 & s2)


def _merge_call(c1, c2):
    s1, k1 = _call_set(c1)
    s2, k2 = _call_set(c2)
    for c, k in ((c1, k1), (c2, k2)):               # a named call wins over ambiguous
        if k == "named":
            return c
    for c, k in ((c1, k1), (c2, k2)):
        if k == "mixed":
            return c
    return "ambiguous(%s)" % "|".join(sorted(s1 & s2)) if s1 & s2 else c1


def _junctions_of(segs, key, args, cass_set, secs, read_segs):
    """Junction events of one read. segs: list of
    [contig, rs, re, strand, mapq, qs, qe, xa_set], sorted by qs."""
    out = []
    # identity evidence: unique cassette segments on this read (+ its mate for short reads)
    srcs = [segs]
    if args.platform == "short" and key[1] in (1, 2):
        srcs.append(read_segs.get((key[0], 3 - key[1]), []))
    uniq = {s[0] for src in srcs for s in src
            if s[0] in cass_set and s[4] >= args.min_mapq and s[2] - s[1] >= args.min_segment}
    for a, b in zip(segs, segs[1:]):                # query-adjacent pairs
        if (a[0] in cass_set) == (b[0] in cass_set):
            continue
        if a[2] - a[1] < args.min_segment or b[2] - b[1] < args.min_segment:
            continue
        g_first = a[0] not in cass_set
        g, c = (a, b) if g_first else (b, a)
        if g[4] < args.min_mapq:                    # genome position must be unique
            continue
        # boundary of the genome segment facing the cassette: read-end side if the
        # genome part comes first in the read, read-start side otherwise;
        # on + strand read-end = ref_end, on - strand read-end = ref_start
        at_ref_end = (g[3] == "+") == g_first
        pos = g[2] if at_ref_end else g[1]
        side = "left" if at_ref_end else "right"
        # same for the cassette segment: its boundary faces the genome part
        # (read-start side if the cassette comes second)
        c_at_end = (c[3] == "+") == (not g_first)
        cpos = c[2] if c_at_end else c[1]
        clen = args.cass_len[c[0]]
        cend = "5'" if cpos < clen - cpos else "3'"
        orient = "+" if g[3] == c[3] else "-"
        if len(uniq) == 1:
            name, cands = next(iter(uniq)), set(uniq)
        elif uniq:
            name, cands = None, set(uniq)           # conflicting
        else:
            name = None
            cands = {c[0]} | (c[7] & cass_set)
            for sc, sqs, sqe in secs.get(key, []):  # secondary alignments overlapping this segment
                if sc in cass_set and sqs < c[6] and sqe > c[5]:
                    cands.add(sc)
        out.append(dict(contig=g[0], pos=pos, side=side, cend=cend, orient=orient,
                        key=key, name=name, cands=cands, kind="split"))
    return out


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _read_fasta(path, names):
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                cur = cur if cur in names else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def _in_seqs(kmer, seqs):
    """kmer (either strand) occurs in one of seqs."""
    rc = kmer.translate(_COMP)[::-1]
    return any(kmer in s or rc in s for s in seqs.values())


def _clip_events(r, args, cass_seqs, cass_unique, genome_seqs=None):
    """Junction events from the soft clips of one uniquely mapped genome read."""
    ct, q = r.cigartuples, r.query_sequence
    if not ct or not q:
        return []
    out = []
    for at_end in (True, False):
        op, n = ct[-1] if at_end else ct[0]
        if op != 4 or n < args.min_clip:
            continue
        k = min(n, args.clip_kmer)
        # clipped bases next to the junction, in reference orientation
        kmer = q[len(q) - n:len(q) - n + k] if at_end else q[n - k:n]
        hits, orient, cpos = set(), None, None
        for c, seq in cass_seqs.items():
            i = seq.find(kmer)
            if i >= 0:
                hits.add(c)
                orient, cpos = "+", i + k / 2
                continue
            i = seq.find(kmer.translate(_COMP)[::-1])
            if i >= 0:
                hits.add(c)
                orient, cpos = "-", i + k / 2
        if not hits:
            continue
        if genome_seqs and _in_seqs(kmer, genome_seqs):
            continue                                   # the clip is genome sequence too (e.g. reads over the origin)
        clen = len(cass_seqs[next(iter(hits))])
        name = next(iter(hits)) if len(hits) == 1 else None
        if name is None and r.is_paired and cass_unique.get(r.query_name) in hits:
            name = cass_unique[r.query_name]           # identity from the uniquely mapped mate
            hits = {name}
        out.append(dict(contig=r.reference_name,
                        pos=r.reference_end if at_end else r.reference_start,
                        side="left" if at_end else "right",
                        cend="5'" if cpos < clen / 2 else "3'", orient=orient,
                        key=_read_key(r), name=name, cands=hits, kind="genome_clip"))
    return out


def _find_unique(kmer, genome_seqs):
    """(contig, index, strand) of the only exact occurrence of kmer (either
    strand) in the genome, else None."""
    hit, rc = None, kmer.translate(_COMP)[::-1]
    for contig, seq in genome_seqs.items():
        for q, strand in ((kmer, "+"), (rc, "-")):
            i = seq.find(q)
            while i >= 0:
                if hit is not None:
                    return None
                hit = (contig, i, strand)
                i = seq.find(q, i + 1)
    return hit


def _cassette_clip_events(r, args, genome_seqs, cache, cass_set, cass_seqs=None):
    """Junction events from the soft clips of a read aligned on a cassette: the
    clipped bases next to the cassette alignment are located in the genome."""
    ct, q = r.cigartuples, r.query_sequence
    if not ct or not q:
        return []
    out = []
    for at_end in (True, False):
        op, n = ct[-1] if at_end else ct[0]
        if op != 4 or n < args.min_clip:
            continue
        k = min(n, args.clip_kmer)
        kmer = q[len(q) - n:len(q) - n + k] if at_end else q[n - k:n]
        if kmer not in cache:
            # clipped bases that are still cassette sequence (a mismatch, or a part
            # of the cassette also found in the genome) say nothing about the genome
            cache[kmer] = None if cass_seqs and _in_seqs(kmer, cass_seqs) else _find_unique(kmer, genome_seqs)
        hit = cache[kmer]
        if hit is None:
            continue
        contig, i, strand = hit
        # read = cassette | genome (clip after the cassette part) or genome | cassette
        if at_end:
            pos, side = (i, "right") if strand == "+" else (i + k, "left")
        else:
            pos, side = (i + k, "left") if strand == "+" else (i, "right")
        cpos = r.reference_end if at_end else r.reference_start
        clen = args.cass_len[r.reference_name]
        name = r.reference_name if r.mapping_quality >= args.min_mapq else None
        out.append(dict(contig=contig, pos=pos, side=side, cend="5'" if cpos < clen / 2 else "3'",
                        orient=strand, key=_read_key(r), name=name,
                        cands={r.reference_name} | (_xa_contigs(r) & cass_set), kind="cassette_clip"))
    return out


def _homology(cass_seqs, genome_seqs, k=25, step=5):
    """Parts of the cassettes identical to the genome: [(cassette, cs, ce, contig, gs, ge)]."""
    out = []
    for c, cseq in cass_seqs.items():
        run = None                                  # [cs, ce, contig, gs, ge]
        for i in range(0, max(len(cseq) - k + 1, 0), step):
            kmer = cseq[i:i + k]
            hit = None
            for contig, gseq in genome_seqs.items():
                j = gseq.find(kmer)
                if j < 0:
                    j = gseq.find(kmer.translate(_COMP)[::-1])
                if j >= 0:
                    hit = (contig, j)
                    break
            if hit and run and hit[0] == run[2] and abs(hit[1] - run[4]) < 2 * step + k:
                run[1], run[3], run[4] = i + k, min(run[3], hit[1]), max(run[4], hit[1] + k)
            elif hit:
                if run:
                    out.append((c, *run))
                run = [i, i + k, hit[0], hit[1], hit[1] + k]
            elif run:
                out.append((c, *run))
                run = None
        if run:
            out.append((c, *run))
    return out


def _hotspots(evidence, args, cass_list):
    """Merge evidence [contig, start, end, kind, read, cassette or None] into hotspots."""
    by_contig = defaultdict(list)
    for e in evidence:
        by_contig[e[0]].append(e)
    spots = []
    for contig, ev in by_contig.items():
        ev.sort(key=lambda e: e[1])
        cur, end = [ev[0]], ev[0][2]
        groups = []
        for e in ev[1:]:
            if e[1] - end <= args.hotspot_gap:
                cur.append(e)
                end = max(end, e[2])
            else:
                groups.append(cur)
                cur, end = [e], e[2]
        groups.append(cur)
        for g in groups:
            reads = {e[4] for e in g}
            if len(reads) < args.hotspot_min_reads:
                continue
            kinds = {k: len({e[4] for e in g if e[3] == k})
                     for k in ("mate", "genome_clip", "cassette_clip", "split")}
            named = Counter()
            for read in reads:                       # one vote per read
                names = {e[5] for e in g if e[4] == read and e[5]}
                if len(names) == 1:
                    named[names.pop()] += 1
            if named and named.most_common(1)[0][1] >= 0.8 * sum(named.values()):
                call = named.most_common(1)[0][0]
            elif len(named) >= 2:
                call = "mixed(%s)" % ",".join(f"{k}:{v}" for k, v in named.most_common())
            else:
                call = "ambiguous"
            jpos = Counter(e[1] for e in g if e[3] != "mate")
            top = ",".join(f"{p}:{n}" for p, n in jpos.most_common(3)) or "."
            s, e_ = min(e[1] for e in g), max(e[2] for e in g)
            spots.append((contig, s, e_, e_ - s, "hotspot", call, len(reads), kinds["mate"],
                          kinds["genome_clip"], kinds["cassette_clip"], kinds["split"], top))
    return spots


def run(args):
    import pysam
    for p in (args.out_junctions, args.out_insertions, args.out_counts):
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    bam = pysam.AlignmentFile(args.bam, "rb")
    refs = list(bam.references)
    lens = dict(zip(refs, bam.lengths))
    cass_list = list(args.cassettes)
    cass_set = set(cass_list)
    for c in cass_list:
        if c not in lens:
            print(f"find_cassette_insertions: warning: cassette {c} not in BAM", file=sys.stderr)
    args.cass_len = lens
    genome = [c for c in refs if c not in cass_set]

    genome_seqs = _read_fasta(args.fasta, set(genome)) if args.fasta else {}
    cass_seqs = _read_fasta(args.fasta, cass_set) if args.fasta else {}
    homology = _homology(cass_seqs, genome_seqs) if genome_seqs else []
    for c, cs, ce, g, gs, ge in homology:
        print(f"WARNING: {c}:{cs}-{ce} ({ce - cs} bp) is identical to the genome at {g}:{gs}-{ge}: "
              f"evidence there needs reads that map uniquely on the cassette")
    # genome windows where the cassette copy can attract reads of the genome copy
    homol_win = defaultdict(list)
    for c, cs, ce, g, gs, ge in homology:
        homol_win[g].append((gs - args.mate_window, ge + args.mate_window))
    cass_clip_events, kmer_cache = [], {}
    read_segs = defaultdict(list)       # (qname, mate) -> segments
    secs = defaultdict(list)            # (qname, mate) -> [(contig, qs, qe)] secondary alns
    interest = set()                    # qnames with any record on a cassette contig
    cass_unique = {}                    # qname -> cassette of a primary alignment with MAPQ >= min-mapq

    def add(r):
        if r.flag & 0x604 or not r.cigartuples:     # unmapped, secondary excluded below
            return
        rs, re_, qs, qe = _clips(r)
        read_segs[_read_key(r)].append(
            [r.reference_name, rs, re_, "-" if r.is_reverse else "+", r.mapping_quality,
             qs, qe, _xa_contigs(r)])

    def add_secondary(r):
        if r.is_unmapped or not r.cigartuples or r.flag & 0x600:
            return
        rs, re_, qs, qe = _clips(r)
        secs[_read_key(r)].append((r.reference_name, qs, qe))

    # ---- pass 1: cassette contigs (counts + reads of interest)
    counts, cass_depth = [], {}
    for c in cass_list:
        if c not in lens:
            counts.append((c, 0, 0, 0, 0.0, 0.0))
            continue
        n_u = n_a = 0
        for r in bam.fetch(c):
            if r.is_unmapped or r.flag & 0x600:
                continue
            interest.add(r.query_name)
            if r.flag & 0x100:
                add_secondary(r)
                continue
            add(r)
            if genome_seqs and "S" in (r.cigarstring or ""):
                cass_clip_events.extend(_cassette_clip_events(r, args, genome_seqs, kmer_cache, cass_set, cass_seqs))
            if not r.flag & 0x800:
                if r.mapping_quality >= args.min_mapq:
                    n_u += 1
                    cass_unique[r.query_name] = c
                else:
                    n_a += 1
        cov = bam.count_coverage(
            c, quality_threshold=0,
            read_callback=lambda r: not (r.flag & 0xF04) and r.mapping_quality >= args.min_mapq)
        depth = [sum(x) for x in zip(*cov)]
        L = lens[c]
        # depth of all primary reads (any MAPQ: shared cassette ends are MAPQ 0), for the report
        cov_all = bam.count_coverage(c, quality_threshold=0, read_callback=lambda r: not (r.flag & 0xF04))
        cass_depth[c] = [sum(x) for x in zip(*cov_all)]
        counts.append((c, L, n_u, n_a, 100.0 * sum(1 for d in depth if d > 0) / max(L, 1),
                       sum(depth) / max(L, 1)))

    # ---- pass 2: genome contigs: reads with a cassette alignment, soft-clipped
    # reads (short reads) and mate links
    clip_cass_seqs = cass_seqs if args.platform == "short" else {}
    clip_events, mate_links = [], Counter()
    unique_mates = defaultdict(list)    # genome contig -> [(start, end, cassette)] mate MAPQ >= min-mapq
    mate_reads = []                     # (contig, start, end, qname, unique cassette or None)
    for g in genome:
        for r in bam.fetch(g):
            if r.is_unmapped or r.flag & 0x200:
                continue
            primary_unique = not r.flag & 0x900 and r.mapping_quality >= args.min_mapq
            if clip_cass_seqs and primary_unique and "S" in (r.cigarstring or ""):
                clip_events.extend(_clip_events(r, args, clip_cass_seqs, cass_unique, genome_seqs))
            if primary_unique and r.is_paired and not r.mate_is_unmapped \
                    and r.next_reference_name in cass_set:
                mate_links[(g, r.reference_start // 1000, r.next_reference_name)] += 1
                mate_reads.append((g, r.reference_start, r.reference_end, r.query_name,
                                   cass_unique.get(r.query_name)))
                if cass_unique.get(r.query_name) == r.next_reference_name:
                    unique_mates[g].append((r.reference_start, r.reference_end, r.next_reference_name))
            if r.query_name not in interest:
                continue
            if r.flag & 0x100:
                add_secondary(r)
            else:
                add(r)
    bam.close()

    # ---- junctions per read
    events = list(clip_events) + cass_clip_events
    n_clip, n_cclip = len(clip_events), len(cass_clip_events)
    for key, segs in read_segs.items():
        if len(segs) < 2:
            continue
        segs.sort(key=lambda s: (s[5], s[6]))
        events.extend(_junctions_of(segs, key, args, cass_set, secs, read_segs))
    # near a genome region identical to a cassette, only reads naming a cassette count
    events = [e for e in events if e["name"]
              or not any(a <= e["pos"] < b for a, b in homol_win.get(e["contig"], ()))]

    by_site = defaultdict(list)
    for e in events:
        by_site[(e["contig"], e["side"])].append(e)
    order = {c: i for i, c in enumerate(refs)}

    junctions = []                      # dicts, sorted later
    for (contig, side), ev in by_site.items():
        for cl in _cluster(ev, args.cluster_tol):
            # one entry per read (a named call wins: split read + its own clip)
            per_read = {}
            for e in cl:
                if e["key"] not in per_read or (e["name"] and not per_read[e["key"]]["name"]):
                    per_read[e["key"]] = e
            reads = list(per_read.values())
            if len(reads) < args.min_reads:
                continue
            pos = sorted(Counter(e["pos"] for e in reads).items(), key=lambda x: (-x[1], x[0]))[0][0]
            sup = Counter(e["name"] for e in reads if e["name"])
            n_u = n_u_reads = sum(sup.values())
            cands = set()
            for e in reads:
                cands |= e["cands"]
            if n_u < args.min_reads and len(sup) <= 1:
                # junction reads cannot tell (shared cassette ends): use the read
                # pairs on the genome side of the junction whose other mate is
                # unique on a cassette (within --mate-window bp)
                for a, b, c in unique_mates.get(contig, []):
                    near = (pos - args.mate_window <= a and b <= pos + 5) if side == "left" \
                        else (pos - 5 <= a and b <= pos + args.mate_window)
                    if near:
                        sup[c] += 1
                n_u = sum(sup.values())
            if n_u >= args.min_reads and sup.most_common(1)[0][1] >= 0.8 * n_u:
                call = sup.most_common(1)[0][0]
            elif len(sup) >= 2:
                call = "mixed(%s)" % ",".join(f"{k}:{v}" for k, v in sup.most_common())
            else:
                call = "ambiguous(%s)" % ("|".join(sorted(cands)) or ".")
            junctions.append(dict(
                contig=contig, pos=pos, side=side, call=call, n_reads=len(reads),
                n_unique=n_u_reads, n_amb=len(reads) - n_u_reads,
                cend=Counter(e["cend"] for e in reads).most_common(1)[0][0],
                orient=Counter(e["orient"] for e in reads).most_common(1)[0][0],
                support=",".join(f"{c}:{sup.get(c, 0)}" for c in cass_list)))
    junctions.sort(key=lambda j: (order[j["contig"]], j["pos"], j["side"]))
    for i, j in enumerate(junctions, 1):
        j["id"] = f"J{i}"

    # ---- pair left/right junctions into sites (smallest span first)
    cand_pairs = []
    for l in junctions:
        if l["side"] != "left":
            continue
        for r in junctions:
            if r["side"] != "right" or r["contig"] != l["contig"]:
                continue
            d = r["pos"] - l["pos"]
            if -args.cluster_tol <= d <= args.max_span and l["orient"] == r["orient"] \
                    and _compatible(l["call"], r["call"]):
                cand_pairs.append((abs(d), l["pos"], l, r))
    cand_pairs.sort(key=lambda x: x[:2])
    used, sites = set(), []
    for _, _, l, r in cand_pairs:
        if l["id"] in used or r["id"] in used:
            continue
        used |= {l["id"], r["id"]}
        s, e = min(l["pos"], r["pos"]), max(l["pos"], r["pos"])
        sites.append((l["contig"], s, e, e - s, "paired", _merge_call(l["call"], r["call"]),
                      l["orient"], l["n_reads"], r["n_reads"], l["id"], r["id"]))
    for j in junctions:
        if j["id"] in used:
            continue
        if j["side"] == "left":
            sites.append((j["contig"], j["pos"], j["pos"] + 1, 0, "left_only", j["call"],
                          j["orient"], j["n_reads"], 0, j["id"], "."))
        else:
            sites.append((j["contig"], j["pos"], j["pos"] + 1, 0, "right_only", j["call"],
                          j["orient"], 0, j["n_reads"], ".", j["id"]))
    sites.sort(key=lambda s: (order[s[0]], s[1], s[2]))

    # ---- write
    with open(args.out_junctions, "w") as out:
        out.write("#contig\tstart\tend\tgenome_side\tcassette\tcassette_end\torientation\t"
                  "n_reads\tn_unique\tn_ambiguous\tsupport\tjunction_id\n")
        for j in junctions:
            out.write(f"{j['contig']}\t{j['pos']}\t{j['pos']+1}\t{j['side']}\t{j['call']}\t"
                      f"{j['cend']}\t{j['orient']}\t{j['n_reads']}\t{j['n_unique']}\t"
                      f"{j['n_amb']}\t{j['support']}\t{j['id']}\n")
    with open(args.out_insertions, "w") as out:
        out.write("#contig\tstart\tend\tlength\tstatus\tcassette\torientation\t"
                  "n_reads_left\tn_reads_right\tleft_junction\tright_junction\n")
        for s in sites:
            out.write("\t".join(str(x) for x in s) + "\n")
    evidence = [[e["contig"], e["pos"], e["pos"] + 1, e["kind"], e["key"][0], e["name"]] for e in events]
    evidence += [[g, a, b, "mate", qn, nm] for g, a, b, qn, nm in mate_reads]
    n_ev = len(evidence)
    evidence = [e for e in evidence if e[5] or not any(a <= e[1] < b for a, b in homol_win.get(e[0], ()))]
    if n_ev - len(evidence):
        print(f"{n_ev - len(evidence)} evidence read(s) dropped: near a genome region identical to a "
              f"cassette, without a read mapping uniquely on the cassette")
    spots = sorted(_hotspots(evidence, args, cass_list), key=lambda x: (order[x[0]], x[1]))
    if args.out_hotspots:
        with open(args.out_hotspots, "w") as out:
            out.write("#contig\tstart\tend\tlength\tstatus\tcassette\tn_reads\tn_mates\t"
                      "n_genome_clip\tn_cassette_clip\tn_split\tjunctions\n")
            for sp in spots:
                out.write("\t".join(str(x) for x in sp) + "\n")
    print(f"{len(spots)} hotspot(s) (>= {args.hotspot_min_reads} reads, gap <= {args.hotspot_gap} bp):")
    for sp in spots:
        print(f"  {sp[0]}:{sp[1]}-{sp[2]} {sp[5]} reads={sp[6]} mates={sp[7]} genome_clip={sp[8]} "
              f"cassette_clip={sp[9]} split={sp[10]} junctions={sp[11]}")
    if args.out_depth:
        with open(args.out_depth, "w") as out:
            out.write("#contig\tstart\tend\tmean_depth\n")
            for c, d in cass_depth.items():
                for b in range(0, len(d), args.depth_bin):
                    w = d[b:b + args.depth_bin]
                    out.write(f"{c}\t{b + 1}\t{b + len(w)}\t{sum(w) / len(w):.1f}\n")
    with open(args.out_counts, "w") as out:
        out.write("#cassette\tlength\treads_unique\treads_ambiguous\t"
                  "breadth_unique_pct\tmean_depth_unique\n")
        for c, L, nu, na, br, md in counts:
            out.write(f"{c}\t{L}\t{nu}\t{na}\t{br:.1f}\t{md:.2f}\n")

    # ---- summary
    print(f"evidence: {n_clip} genome read(s) clipped into a cassette, {n_cclip} cassette read(s) "
          f"whose clip maps once in the genome, {len(events) - n_clip - n_cclip} split-read junction(s), "
          f"{len(mate_reads)} pair(s) with a unique genome mate")
    if mate_links:
        print("mate links (unique genome mate <-> cassette mate), top genome kb:")
        for (g, kb, c), n in mate_links.most_common(10):
            print(f"  {g}:{kb * 1000}-{(kb + 1) * 1000}  {c}  {n} pair(s)")
    elif args.platform == "short":
        print("mate links: none (no pair with one mate on the genome and one on a cassette)")
    n_pair = sum(1 for s in sites if s[4] == "paired")
    print(f"{len(junctions)} junction(s); {len(sites)} site(s) ({n_pair} paired, "
          f"{len(sites) - n_pair} single-sided)")
    for s in sites:
        print(f"  {s[0]}:{s[1]}-{s[2]} {s[4]} {s[5]} {s[6]} len={s[3]}")
    for c, L, nu, na, br, md in counts:
        print(f"  {c}: unique={nu} ambiguous={na} breadth={br:.1f}% depth={md:.1f}x")
    print(f"done -> {args.out_junctions}, {args.out_insertions}, {args.out_counts}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="sorted, indexed BAM (genome + cassette contigs)")
    ap.add_argument("--cassettes", required=True, nargs="+",
                    help="names of the cassette contigs in the BAM (all others are genome)")
    ap.add_argument("--platform", required=True, choices=["ont", "short"],
                    help="ont: long reads; short: paired reads (mate evidence + XA tag)")
    ap.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality for unique evidence [20]")
    ap.add_argument("--min-segment", type=int, default=30,
                    help="minimum aligned ref length of each side of a junction (bp) [30]")
    ap.add_argument("--cluster-tol", type=int, default=10,
                    help="position tolerance when clustering junctions (bp) [10]")
    ap.add_argument("--min-reads", type=int, default=2, help="minimum reads supporting a junction [2]")
    ap.add_argument("--max-span", type=int, default=20000,
                    help="maximum genome span between paired left/right junctions (bp) [20000]")
    ap.add_argument("--fasta", default=None,
                    help="genome + cassettes FASTA: enables soft-clip evidence")
    ap.add_argument("--min-clip", type=int, default=12,
                    help="minimum soft clip (bp) used as junction evidence [12]")
    ap.add_argument("--clip-kmer", type=int, default=20,
                    help="clipped bases (max) matched against the cassettes [20]")
    ap.add_argument("--mate-window", type=int, default=1000,
                    help="short reads: mates this close to a junction name its cassette [1000]")
    ap.add_argument("--out-junctions", required=True, help="junctions TSV")
    ap.add_argument("--out-insertions", required=True, help="insertion sites TSV")
    ap.add_argument("--out-counts", required=True, help="per-cassette read counts TSV")
    ap.add_argument("--out-hotspots", default=None, help="hotspots TSV (all evidence)")
    ap.add_argument("--hotspot-gap", type=int, default=1000,
                    help="merge hotspot evidence closer than this (bp) [1000]")
    ap.add_argument("--hotspot-min-reads", type=int, default=3,
                    help="minimum reads in a hotspot [3]")
    ap.add_argument("--out-depth", default=None,
                    help="depth along each cassette, all primary reads (depth_bins format)")
    ap.add_argument("--depth-bin", type=int, default=10, help="bin size of --out-depth (bp) [10]")
    args = ap.parse_args()
    try:
        run(args)
    except ImportError as e:
        sys.exit(f"find_cassette_insertions: missing dependency: {e}")
    except (OSError, ValueError) as e:
        sys.exit(f"find_cassette_insertions: error: {e}")


if __name__ == "__main__":
    main()
