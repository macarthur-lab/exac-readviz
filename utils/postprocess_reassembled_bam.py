import argparse
import logging
import pysam
import re


def interval_union(a, b):
    """Returns a 2-tuple representing the union of two intervals.
    Args:
        a:  2-tuple containing integer coordinates representing a half-open interval.
        b:  2-tuple containing integer coordinates representing the other half-open interval.
    Returns:
        2-tuple containing integer coordinates representing the union(a, b) as a half-open interval.
    """

    return min(a[0], b[0]), max(a[1], b[1])


def do_intervals_intersect(a, b):
    """Returns true if the given 2-tuples overlap.

    Args:
        a:  2-tuple containing integer coordinates representing a half-open interval.
        b:  2-tuple containing integer coordinates representing the other half-open interval.
    Returns:
        True if a and b overlap.
    """

    return a[0] < b[1] and a[1] > b[0]


def read_name_hash(read_name):
    """Takes a read name from an input bam and shortens & obfuscates it"""
    return str(abs(hash(read_name)) % 10**9)  # 9-digit read name


def postprocess_bam(input_bam_path, output_bam_path, chrom, pos, ref, alt):
    """Copies the input_bam to the output_bam while discarding extraneous
    or sensitive information. Leaves only the minimum required header,
    obfuscates and downsizes read names, discards all tags (including read groups).

    If the input bam doesn't have any reads, the output bam won't be written.

    Args:
        input_bam_path:  input bam path
        output_bam_path:  output bam path
        chrom: chromosome (eg. '1' or 'X')
        pos: minrep'ed variant position integer (eg. 12345)
        ref: minrep'ed ref allele (eg. 'A', 'ACT', etc.)
        alt: minrep'ed alt allele (eg. 'GCT', 'C', etc.)

    Return:
        2-tuple with (is_empty, artificial_haplotype_counter) where
           is_empty: is True if the input_bam was empty
           artificial_haplotype_counter: the number of artificial haplotypes found in the input bam
           artificial_haplotypes_deleted_counter: the number of artificial haplotypes discarded because they
                overlap other artificial haplotypes in a way that might cause double-counting of reads.

    """

    # compute variant start, end reference coords (half-open)
    variant_start = pos
    variant_end = pos + len(ref)

    # This counter is used as a sanity check that HC added at least one artificial haplotype (typically it adds
    # 2*n of these where n is the number of SNPs in the region).
    artificial_haplotype_counter = 0

    # artificial haplotype coords are half-open (eg. (start=83, end=93) has length 10)
    union_of_artificial_haplotypes_that_overlap_variant = (1e9, 0)  # union of genomic intervals spanned by artificial haplotypes that overlap the variant

    # iterate over the reads
    raw_reads = []  # maps each artificial haplotype id (eg. HC tag value) to the list of reads assigned to this haplotype (eg. that have this id in their HC tag)
    ibam = pysam.AlignmentFile(input_bam_path, "rb")
    for r in ibam:
        tags = dict(r.tags)
        if tags.get('RG') == "ArtificialHaplotype":
            # handle reads that are actually artificial haplotypes
            artificial_haplotype_counter += 1

            # check whether the artificial haplotype overlaps the variant
            if r.reference_start <= variant_end and r.reference_end >= variant_start:
                union_of_artificial_haplotypes_that_overlap_variant = interval_union(
                    (r.reference_start, r.reference_end),
                    union_of_artificial_haplotypes_that_overlap_variant)
        else:
            # this is a regular read, not an artificial haplotype
            raw_reads.append(r)

    if not raw_reads:
        is_bam_empty = True
        return (is_bam_empty, artificial_haplotype_counter, artificial_haplotypes_deleted_counter)

    # write out the bam
    reference_sequences = []
    for reference_id in range(len(ibam.header['SQ'])):
        d = {}  # dcitionary of ref. sequences
        reference_sequences.append(d)
        for key, value in ibam.header['SQ'][reference_id].items():
            if key in ["SN", "LN"]:
                d[key] = value

    header = {'HD': { 'VN': '1.4', 'SO': 'coordinate' },
              'SQ': reference_sequences,
              'RG': [],
              }

    is_bam_empty = True
    obam = pysam.AlignmentFile(output_bam_path, "wb", header=header)
    for r in raw_reads:
        if not (r.reference_start >= union_of_artificial_haplotypes_that_overlap_variant[0] and \
                r.reference_end <= union_of_artificial_haplotypes_that_overlap_variant[1]):
            continue  # the read must be fully contained within the haplotype interval to avoid the issue with coverage bumps

        # copy info from r to s
        s = pysam.AlignedSegment()
        s.query_name = read_name_hash(r.query_name)
        s.query_sequence = r.query_sequence
        s.flag = r.flag
        s.reference_id = r.reference_id  # since the bam should only have reads from one chromosome, there will always be just 1 reference entry in the output bam header
        s.reference_start = r.reference_start
        s.mapping_quality = r.mapping_quality
        s.cigar = r.cigar
        s.next_reference_id = r.next_reference_id
        s.next_reference_start = r.next_reference_start
        s.template_length = r.template_length
        s.query_qualities = r.query_qualities
        
        obam.write(s)
        is_bam_empty = False

    if obam is not None:
        obam.close()

    return (is_bam_empty, artificial_haplotype_counter)


if __name__ == "__main__":
    p = argparse.ArgumentParser("Takes an HC output bam and discards non-essential header fields and tags, obfuscates read names, etc.")
    p.add_argument("-i", "--input-bam", help=".bam output from HaplotypeCaller", required=True)
    p.add_argument("-o", "--output-bam", help="Postprocessed bam", required=True)
    args = p.parse_args()

    match = re.search('([0-9XY]{1,2})-([0-9]{1,9})-([ACGTN]+)-([ACGTN]+)', args.input_bam)
    chrom = match.group(1)
    pos = int(match.group(2))
    ref = match.group(3)
    alt = match.group(4)
    postprocess_bam(args.input_bam, args.output_bam, chrom, pos, ref, alt)
