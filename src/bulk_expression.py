import argparse
import os
import logging
import subprocess
from pathlib import Path

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('fq_1', type=str,
    help='tumor rna fq 1. if multiple fastqs then seperate them with "," character')

parser.add_argument('fq_2', type=str,
    help='tumor rna fq 2. if multiple fastqs then seperate them with "," character')

parser.add_argument('--out-dir', type=str, default='output',
    help='output directory')

parser.add_argument('--star-index', type=str,
    help='location of star index directory')

parser.add_argument('--gtf', type=str,
    help='Location of gtf')

parser.add_argument('--gene-info', type=str,
    help='location of gene info file')

parser.add_argument('--cpu', type=int, default=16,
    help='location of gene info file')

parser.add_argument('--compress-featurecounts-script', type=str,
    help='location of compress featurecounts script')

parser.add_argument('--generate-fpkm-script', type=str,
    help='location of generate fpkm script')


args = parser.parse_args()


def star_align(fq_1, fq_2, star_index, out_dir, cpu):
    out_prefix = out_dir + '/' if out_dir[-1] != '/' else out_dir
    pieces = [
        "STAR ",
        f"--readFilesIn {fq_1} {fq_2} ",
        # Most parameters follow GDC
        "--alignIntronMax 1000000 ",
        "--alignIntronMin 20 ",
        "--alignMatesGapMax 1000000 ",
        "--alignSJDBoverhangMin 1 ",
        "--alignSJoverhangMin 8 ",
        "--alignSoftClipAtReferenceEnds Yes ",

        # Follow arriba's recommendation regarding chimera parameters
        # Ref: https://arriba.readthedocs.io/en/latest/workflow/
        "--chimJunctionOverhangMin 10 ",
        "--chimMainSegmentMultNmax 1 ",
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip ",
        "--chimOutJunctionFormat 1 ",
        "--chimSegmentMin 10 ",
        "--chimScoreMin 1",
        "--chimScoreDropMax 30 ",
        "--chimScoreJunctionNonGTAG 0 ",
        "--chimScoreSeparation 1 ",
        "--alignSJstitchMismatchNmax 5 -1 5 5 ",
        "--chimSegmentReadGapMax 3 ",

        f"--genomeDir {star_index} ",
        "--genomeLoad NoSharedMemory ",
        "--limitBAMsortRAM 0 ",
        "--limitSjdbInsertNsj 1200000 ",
        f"--outFileNamePrefix {out_prefix} ",
        "--outFilterIntronMotifs None ",
        "--outFilterMatchNminOverLread 0.33 ",
        "--outFilterMismatchNmax 999 ",
        "--outFilterMismatchNoverLmax 0.1 ",
        "--outFilterMultimapNmax 20 ",
        "--outFilterScoreMinOverLread 0.33 ",
        "--outFilterType BySJout ",
        "--outSAMattributes NH HI AS nM NM ch ",
        "--outSAMstrandField intronMotif ",
        "--outSAMtype BAM SortedByCoordinate ",
        "--outSAMunmapped Within ",
        "--quantMode TranscriptomeSAM GeneCounts ",
        "--readFilesCommand zcat ",
        f"--runThreadN {cpu} ",
        "--twopassMode Basic "
    ]
    return ''.join(pieces)


def index_bam(bam):
    return f'samtools index {bam} {bam}.bai'


def featurecounts_unstranded_readcount(bam, gtf, output_fp, cpu):
    pieces = [
        'featureCounts ',
        '-g gene_id ',  # feature id (-i in htseq)
        '-t exon ',  # feature type (-t in htseq)
        f'-T {cpu} ',
        '-Q 10 ',  # htseq set this minimal mapping quality by default
        '-p ',  # pair-end reads are considered one fragment; default HTSeq behavior
        '-B ',  # both reads of a read pair need to be mapped
        f'-a {gtf} ',
        f'-o {output_fp} {bam}'
    ]
    return ''.join(pieces)


def compress_featurecounts(in_fp, out_fp, compress_featurecounts_script):
    return f'python {compress_featurecounts_script} {in_fp} | gzip -9 -c > {out_fp}'


def generate_fpkm(generate_fpkm_script, gene_info, in_fp, out_fp):
    return f'python {generate_fpkm_script} {gene_info} {in_fp} {out_fp}'


def run_bulk_expression(
        fq_1, fq_2, star_index, gtf, gene_info, out_dir, cpu,
        compress_featurecounts_script, generate_fpkm_script):
    cpu = str(cpu)

    logging.info('aligning fastqs')
    star_out_dir = './star'
    Path(star_out_dir).mkdir(parents=True, exist_ok=True)
    cmd = star_align(fq_1, fq_2, star_index, star_out_dir, cpu)
    logging.info(f'executing command: {cmd}')
    subprocess.check_output(cmd, shell=True)

    logging.info('indexing bam')
    bam_fp = os.path.join(star_out_dir, 'Aligned.sortedByCoord.out.bam')
    cmd = index_bam(bam_fp)
    logging.info(f'executing command: {cmd}')
    subprocess.check_output(cmd, shell=True)

    logging.info('generating featurecounts')
    featurecounts_fp = os.path.join(
        out_dir, 'featurecounts_unstranded_readcount.tsv')
    cmd = featurecounts_unstranded_readcount(
        bam_fp, gtf, featurecounts_fp, cpu)
    logging.info(f'executing command: {cmd}')
    subprocess.check_output(cmd, shell=True)

    logging.info('compressing featurecounts')
    featurecounts_compressed_fp = f'{featurecounts_fp}.gz'
    cmd = compress_featurecounts(
        featurecounts_fp, featurecounts_compressed_fp,
        compress_featurecounts_script)
    logging.info(f'executing command: {cmd}')
    subprocess.check_output(cmd, shell=True)

    logging.info('generate fpkm')
    fpkm_fp = os.path.join(out_dir, 'readcount_and_fpkm.tsv.gz')
    cmd = generate_fpkm(
        generate_fpkm_script, gene_info, featurecounts_compressed_fp, fpkm_fp)
    logging.info(f'executing command: {cmd}')
    subprocess.check_output(cmd, shell=True)


def main():
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    run_bulk_expression(
        args.fq_1, args.fq_2, args.star_index, args.gtf, args.gene_info,
        args.out_dir, args.cpu, args.compress_featurecounts_script,
        args.generate_fpkm_script)


if __name__ == '__main__':
    main()
