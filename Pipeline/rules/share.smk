# ======================================================================
# rules/share.smk :: one folder to share the run with the lab
# (run_share_bundle; runs last, after the report)
#
#   <share_dir>/bams/<sample>.bam(.bai)      every analysed sample
#   <share_dir>/references/<reference>.gb    one per reference: the staged
#                                            copy (contig names = the BAMs')
#   <share_dir>/<run>_bams.tar               BAMs are already compressed:
#                                            bundled, not re-gzipped
#   <share_dir>/<run>_references.tar.gz
#   <share_dir>/report.html, summary.tsv
# Files in bams/ and references/ are hard links when the filesystem allows
# it (no extra disk space), else copies.
# ======================================================================

SHARE_DIR = os.path.join(config.get("share_dir") or OUT_DIR + "share", "")
RUN_NAME  = re.sub(r"[^\w.-]", "_", os.path.basename(os.path.abspath(os.getcwd())))


def _ref_name(sample):
    return re.sub(r"[^\w.-]", "_", os.path.splitext(os.path.basename(ref_gb_path(sample)))[0])


def _share_inputs(wildcards):
    OK = ok_samples()
    first = {}                                  # reference name -> first sample using it
    for s in OK:
        first.setdefault(_ref_name(s), s)
    d = {
        "bams": expand(OUT_DIR + "{sample}/{sample}.sorted.bam", sample=OK),
        "bais": expand(OUT_DIR + "{sample}/{sample}.sorted.bam.bai", sample=OK),
        "gbs":  [OUT_DIR + f"{s}/{s}.ref.gb" for s in first.values()],
        "summary": OUT_DIR + "summary.tsv",
    }
    if RUN_REPORT:
        d["report"] = OUT_DIR + "report.html"    # also makes this the last step
    return d


rule share_bundle:
    input:
        unpack(_share_inputs)
    output:
        bams_tar = SHARE_DIR + RUN_NAME + "_bams.tar",
        refs_tar = SHARE_DIR + RUN_NAME + "_references.tar.gz",
    log:
        OUT_DIR + "logs/share_bundle.log",
    params:
        share   = SHARE_DIR,
        samples = lambda wc: ok_samples(),
        refs    = lambda wc: {_ref_name(s): s for s in reversed(ok_samples())},
    run:
        import shutil

        def place(src, dst):
            if os.path.lexists(dst):
                os.remove(dst)
            try:
                os.link(src, dst)                   # hard link: no extra space
            except OSError:
                shutil.copy2(src, dst)

        bam_dir, ref_dir = params.share + "bams", params.share + "references"
        for d in (bam_dir, ref_dir):
            shutil.rmtree(d, ignore_errors=True)    # no leftovers of an earlier run
            os.makedirs(d)
        with open(log[0], "w") as lg:
            for s in params.samples:
                place(OUT_DIR + f"{s}/{s}.sorted.bam", f"{bam_dir}/{s}.bam")
                place(OUT_DIR + f"{s}/{s}.sorted.bam.bai", f"{bam_dir}/{s}.bam.bai")
            for ref, s in sorted(params.refs.items()):
                place(OUT_DIR + f"{s}/{s}.ref.gb", f"{ref_dir}/{ref}.gb")
            for f in (input.summary, input.get("report")):
                if f:
                    place(f, params.share + os.path.basename(f))
            lg.write(f"{len(params.samples)} BAM(s), {len(params.refs)} reference(s) -> {params.share}\n")
        shell("tar -C {params.share} -cf {output.bams_tar} bams 2>> {log}")
        shell("tar -C {params.share} -czf {output.refs_tar} references 2>> {log}")
