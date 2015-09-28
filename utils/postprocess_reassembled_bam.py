import argparse
import pysam

def read_name_hash(read_name):
    """Takes a read name from an input bam and shortens & obfuscates it"""
    return str(abs(hash(read_name)) % 10**9)  # 9-digit read name

def postprocess_bam(input_bam_path, output_bam_path):
    """Copies the input_bam to the output_bam while discarding extraneous
    or sensitive information. Leaves only the minimum required header,
    obfuscates and downsizes read names, discards all tags (including read groups).

    If the input bam doesn't have any reads, the output bam won't be written.

    Return:
        2-tuple with (is_empty, artificial_haplotype_counter) where
           is_empty: is True if the input_bam was empty
           artificial_haplotype_counter: the number of artificial haplotypes found in the input bam

    """

    # This counter is used as a sanity check that HC added at least one read
    # representing the assembled haplotype (typically it adds 2*n such reads
    # where n is the number of SNPs in the region).
    artificial_haplotype_counter = 0

    ibam = pysam.AlignmentFile(input_bam_path, "rb")

    reference_sequences = []
    for reference_id in range(len(ibam.header['SQ'])):
        d = {}
        reference_sequences.append(d)
        for key, value in ibam.header['SQ'][reference_id].items():
            if key in ["SN", "LN"]:
                d[key] = value

    header = {'HD': { 'VN': '1.4', 'SO': 'coordinate' },
              'SQ': reference_sequences,
              'RG': [],
              }
    obam = pysam.AlignmentFile(output_bam_path, "wb", header=header)

    # iterate over the reads
    is_input_bam_empty = True
    for r in ibam:
        #if r.get_tag('RG') == "ArtificialHaplotype":  # this doesn't work with the old version of pysam installed on the cluster
        if dict(r.tags).get('RG') == "ArtificialHaplotype":
            artificial_haplotype_counter += 1
            continue

        # copy info from r to s
        s = pysam.AlignedSegment()
        s.query_name = read_name_hash(r.query_name)
        s.query_sequence = r.query_sequence
        s.flag = r.flag
        s.reference_id = r.reference_id  # since the bam should only have reads from one chromosome, there will always be just 1 chromosome entry in the header, and so this reference_id can always be 0.
        s.reference_start = r.reference_start
        s.mapping_quality = r.mapping_quality
        s.cigar = r.cigar
        s.next_reference_id = r.next_reference_id
        s.next_reference_start = r.next_reference_start
        s.template_length = r.template_length
        s.query_qualities = r.query_qualities

        obam.write(s)

        is_input_bam_empty = False

    if obam is not None:
        obam.close()

    if not is_input_bam_empty:
        assert artificial_haplotype_counter > 0, \
            "Expected HaplotypeCaller to add at least one record with " \
            "RG == 'ArtificialHaplotype'. " \
            "%(input_bam_path)s => %(output_bam_path)" % locals()

    return (is_input_bam_empty, artificial_haplotype_counter)

if __name__ == "__main__":
    p = argparse.ArgumentParser("Takes an HC output bam and discards non-essential header fields and tags, obfuscates read names, etc.")
    p.add_argument("-i", "--input-bam", help=".bam output from HaplotypeCaller", required=True)
    p.add_argument("-o", "--output-bam", help="Postprocessed bam", required=True)
    args = p.parse_args()

    postprocess_bam(args.input_bam, args.output_bam)
