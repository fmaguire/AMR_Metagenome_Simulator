#!/usr/bin/env python

import pysam
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Add proper reference names to ARTs BAM")
    parser.add_argument('-i', '--input_art_bam', required = True, help="BAM generated by ART sim")
    parser.add_argument('-o', '--output_bam', required = True, help="Path to output fixed BAM")
    args = parser.parse_args()

    art_bam = pysam.AlignmentFile(args.input_art_bam)

    fixed_header = art_bam.header.to_dict()
    fixed_reference_names = []
    for ref in fixed_header['SQ']:
        #modified the headers to remove : and - as well as they are reserved characters in samtools view
        fixed_reference_names.append({'SN': ref['SN'].split()[0].replace(":","|").replace("-","|"),
                                      'LN': ref['LN']})
    fixed_header['SQ'] = fixed_reference_names
    fixed_header = pysam.AlignmentHeader.from_dict(fixed_header)
    fixed_bam = pysam.AlignmentFile(args.output_bam, "wb", header=fixed_header)

    for record in art_bam:
        record_dict = record.to_dict()
        #modified the headers to remove : and - as well as they are reserved characters in samtools view
        #changed split to rsplit in case of multiple - in refname
        record_dict['ref_name'] = record_dict['name'].rsplit("-",1)[0].replace(":","|").replace("-","|")
        fixed_record = pysam.AlignedSegment.from_dict(record_dict, header=fixed_header)
        fixed_bam.write(fixed_record)

#
