import click
import contextlib
import errno
import itertools
import logging
import os
import pandas as pd
import shutil
import subprocess
import sys
import tempfile


class Metadata(object):
    def __init__(self, table, sample_id, group_id, barcode, header=True):
        """
        Args:
            table (str): file path to metadata table; csv or tsv
            sample_id (str or int): sample id column header label or 0-based column
            group_id (str or int): group id column header label or 0-based column
            barcode (str or int): barcode column header label or 0-based column
            header (boolean): header presence

        Notes:
            Exits if metadata table does not validate, e.g. sample_id header label was not found.
        """
        self.table = os.path.abspath(table)
        self.sample_id = sample_id
        self.group_id = group_id
        self.barcode = barcode
        self.dataframe = self.validate_table(header)
        logging.info("Found {samples} samples across {groups} groups within {table}".format(
            samples=len(self.dataframe), groups=len(self.dataframe[self.group_id].unique()),
            table=table))

    def validate_table(self, header):
        """Validates input metadata table and verifies that sample_id, group_id, and barcode
        exist within the table.

        Args:
            header (boolean): header presence in metadata table

        Returns:
            pandas.DataFrame
        """
        df = pd.read_table(self.table, header=None if not header else 0, sep=None, engine="python")
        if self.sample_id not in df.columns:
            logging.critical("--sample-id '%s' was not found in the header of %s %s" %
                            (self.sample_id, self.table, list(df.columns)))
            sys.exit(1)
        if self.group_id not in df.columns:
            logging.critical("--group-id '%s' was not found in the header of %s %s" %
                            (self.group_id, self.table, list(df.columns)))
            sys.exit(1)
        if self.barcode not in df.columns:
            logging.critical("--barcode '%s' was not found in the header of %s %s" %
                            (self.barcode, self.table, list(df.columns)))
            sys.exit(1)
        return df

    def barcodes(self, out_file):
        """Generates barcode TSV of <sample> <sequence> from the metadata table.

        Args:
            out_file (str): file path to barcodes TSV file

        Returns:
            str
        """
        self.dataframe.to_csv(out_file, columns=[self.sample_id, self.barcode], header=False,
                              index=False, sep="\t")
        return os.path.abspath(out_file)


class FastqMultx(object):
    def __init__(self, r1, barcodes, r2=None, i1=None, out_dir=None, mismatches=0, distance=2,
                 quality=0):
        """
        Args:
            tsv (str): file path of barcodes TSV with name<tab>barcode
            i1 (str): file path of Index Fastq
            r1 (str): file path of Read 1
            r2 (str): file path of Read 2
            out_dir (str): directory path for demultiplexed reads
            mismatches (Optional[int]): allowable mismatches in barcode as long as they are unique
            distance (Optional[int]): minimum distance between best and next best match
            quality (Optional[int]): require a minimum phred quality to accept a barcode base
        """
        self.r1 = os.path.abspath(r1)
        self.barcodes = barcodes

        if not r2 and not i1:
            if "_R1" not in r1:
                logging.critical("'_R1' ('_R2' and '_I1') must exist within the read name(s).")
                raise OSError(errno.ENOENT, "Invalid naming convention", r1)
        if not r2:
            r2 = r1.replace("_R1", "_R2")
            if not os.path.exists(r2):
                logging.critical("R2 was not found. It should be located in the same directory as R1.")
                raise OSError(errno.ENOENT, os.strerror(errno.ENOENT), r2)
        if not i1:
            i1 = r1.replace("_R1", "_I1")
            if not os.path.exists(i1):
                logging.critical("I1 was not found. It should be located in the same directory as R1.")
                raise OSError(errno.ENOENT, os.strerror(errno.ENOENT), i1)

        self.r2 = os.path.abspath(r2)
        self.i1 = os.path.abspath(i1)

        if not out_dir:
            out_dir = os.path.dirname(self.r1)

        self.out_dir = os.path.abspath(out_dir)
        self.stats_file = os.path.join(self.out_dir, "demultiplex_stats.txt")
        self.mismatches = mismatches
        self.distance = distance
        self.quality = quality
        self.output_files = self.expected_output_files()

    def expected_output_files(self):
        """Given the barcodes TSV and output directory for fastq-multx, return the list of output
        files we're expecting.

        Expected demux format is out_dir/%_R1.fastq, out_dir/%_R2.fastq

        Returns:
            list: [[r1 file path, r2 file path],]
        """
        out_files = []
        with open(self.barcodes, 'rU') as fh:
            for line in fh:
                # sample name\tbarcode
                toks = line.strip().split("\t")
                if len(toks) < 2:
                    raise OSError(errno.ENOENT, "Invalid barcodes file; expects <sample>\\t<seq>", r1)
                out_files.append([os.path.join(self.out_dir, "%s_R1.fastq" % toks[0]),
                                  os.path.join(self.out_dir, "%s_R2.fastq" % toks[0])])
        return out_files

    def run(self):
        """Runs `fastq-multx` across barcodes and paired-end reads, removes unmatched fastqs, and
        validates that expected output files are present.
        """
        if not os.path.exists(self.out_dir):
            os.path.makedirs(self.out_dir)

        # demultiplexed using pattern: %_R{1,2}.fastq
        cmd = ("fastq-multx -B {barcodes} {i1} {r1} {r2} "
               "-o {dir}/%_I1.fastq -o {dir}/%_R1.fastq -o {dir}/%_R2.fastq "
               "-x -m {m} -d {d} -q {q} > {stats} 2> {null}").format(barcodes=self.barcodes,
                   i1=self.i1, r1=self.r1, r2=self.r2, dir=self.out_dir, m=self.mismatches,
                   d=self.distance, q=self.quality, stats=self.stats_file, null=os.devnull)
        logging.info(("Demultiplexing (mismatches={mismatches}, distance={distance}, "
                      "quality={quality})".format(mismatches=self.mismatches,
                                                  distance=self.distance, quality=self.quality)))
        logging.debug(cmd)
        subprocess.check_call(cmd, shell=True)
        # remove unmatched R1 and R2
        if os.path.exists(os.path.join(self.out_dir, "unmatched_R1.fastq")):
            os.remove(os.path.join(self.out_dir, "unmatched_R1.fastq"))
        if os.path.exists(os.path.join(self.out_dir, "unmatched_R2.fastq")):
            os.remove(os.path.join(self.out_dir, "unmatched_R2.fastq"))
        # validate that we have the right files here
        for (out_r1, out_r2) in self.output_files:
            if not os.path.exists(out_r1):
                logging.error("Expected file %s was not found!" % out_r1)
            if not os.path.exists(out_r2):
                logging.error("Expected file %s was not found!" % out_r2)


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("r1", type=click.Path(exists=True, dir_okay=False))
@click.argument("table", type=click.Path(exists=True, dir_okay=False))
@click.option("--i1", type=click.Path(exists=True, dir_okay=False),
              help=("index read FASTQ; if the file is located in the same directory as 'R1' and "
                    "contains '_I1' in the read name, specifying is optional"))
@click.option("--r2", type=click.Path(exists=True, dir_okay=False),
              help=("second read (R2) FASTQ; if the file is located in the same directory as 'R1' "
                    "and contains '_R2' in the read name, specifying is optional"))
@click.option("-a", "--output-action", type=click.Choice(["undetermined", "groupid"]),
              default="undetermined", show_default=True,
              help=("affects the final output file; 'undetermined' write Undetermined I1, R1, and "
                    "R2 within group named directories, whereas 'groupid' outputs files prefixed "
                    "with individual group IDs"))
@click.option("-o", "--out", default=".", show_default=True,
              help="output directory for demultiplexed files")
@click.option("--header/--no-header", default=True, help=("metadata table header: when there's no "
              "header, sample-id and group-id should be 0-based column numbers"))
@click.option("--sample-id", default="sampleid", show_default=True,
              help="metadata header label for sample id column or 0-based column number")
@click.option("--group-id", default="groupid", show_default=True,
              help="metadata header label for group id column or 0-based column number")
@click.option("--barcode", default="barcode", show_default=True,
              help="metadata header label for barcode column or 0-based column number")
@click.option("-m", "--barcode-mismatches", type=int, default=1, show_default=True,
              help="allowed mismatches as long as barcodes are unique")
@click.option("-d", "--distance", type=int, default=2, show_default=True,
              help=("helps in determining uniqueness when allowing mismatches; minimum distance "
                    "between best and next best match"))
@click.option("-q", "--quality", type=int, default=0, show_default=True,
              help="require a minimum phred quality to accept a barcode base")
@click.option("--stats-file", default=os.devnull,
              help="file to save fastq-multx stdout which contains per sample read counts")
def group_demux(r1, table, i1, r2, output_action, out, header, sample_id, group_id, barcode,
                barcode_mismatches, distance, quality, stats_file):
    """Demultiplex FASTQs using `fastq-multx` (`conda install -c bioconda fastq-multx`) into grouped
    FASTQs rather than individual samples. Effectively, this subsets 'Undetermined' into smaller
    groups of 'Undetermined' files.

    Reads `R1` and the `metadata` sheet -- sample and group ID mapping to barcode sequences -- and
    depends on the user to define `--sample-id`, `--group-id`, and `--barcode` to match the header
    or column number.
    """
    logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M",
                        format="[%(asctime)s %(levelname)s] %(message)s")
    if not header:
        try:
            sample_id = int(sample_id)
            group_id = int(group_id)
            barcode = int(barcode)
        except ValueError:
            logging.warning(("When not using a header, the sample, group, and barcode columns "
                             "should be integers"))
            sys.exit(1)
    if not os.path.exists(out):
        os.makedirs(out)
    metadata = Metadata(table, sample_id, group_id, barcode, header)
    group_files = {}
    with temp_dir() as td:
        barcodes = metadata.barcodes(os.path.join(td, "barcodes.tsv"))
        fmultx = FastqMultx(r1, barcodes, r2, i1, td, barcode_mismatches, distance, quality)
        fmultx.run()
        # save the stats if the user would like a copy
        if not stats_file == os.devnull:
            shutil.copy(fmultx.stats_file, stats_file)
        # concatenate the files
        cmds = []
        for gid, gdf in metadata.dataframe.groupby([group_id]):
            if output_action == "groupid":
                r1_result_file = os.path.join(out, "{group}_R1.fastq".format(group=gid))
                r2_result_file = os.path.join(out, "{group}_R2.fastq".format(group=gid))
                i1_result_file = os.path.join(out, "{group}_I1.fastq".format(group=gid))
            # "undetermined"
            else:
                r1_result_file = os.path.join(out, gid, "Undetermined_R1.fastq")
                r2_result_file = os.path.join(out, gid, "Undetermined_R2.fastq")
                i1_result_file = os.path.join(out, gid, "Undetermined_I1.fastq")
                if not os.path.exists(os.path.join(out, gid)):
                    os.makedirs(os.path.join(out, gid))
            group_files[gid] = r1_result_file
            # per sample fastq for r1, r2, and i1
            r1_fastqs = [os.path.join(td, "%s_R1.fastq" % sample) for sample in gdf[sample_id]]
            r2_fastqs = [os.path.join(td, "%s_R2.fastq" % sample) for sample in gdf[sample_id]]
            i1_fastqs = [os.path.join(td, "%s_I1.fastq" % sample) for sample in gdf[sample_id]]
            # per group cat command
            cmds.append("cat {fastqs} > {result}".format(fastqs=" ".join(r1_fastqs),
                                                         result=r1_result_file))
            cmds.append("cat {fastqs} > {result}".format(fastqs=" ".join(r2_fastqs),
                                                         result=r2_result_file))
            cmds.append("cat {fastqs} > {result}".format(fastqs=" ".join(i1_fastqs),
                                                         result=i1_result_file))
        logging.info("Joining reads across groups")
        # parallelize concat (at least) across R1, R2, and I1 per group
        running_groups = [(subprocess.Popen(cmd, shell=True) for cmd in cmds)] * 3
        for processes in itertools.zip_longest(*running_groups):
            for p in filter(None, processes):
                p.wait()
                if p.returncode != 0:
                    logging.warning("Joining the reads failed")
        # validate the group FASTQs
        validate_group_fastqs(group_files, fmultx.stats_file, metadata.dataframe, sample_id,
                              group_id)
    logging.info("Processing complete")


@contextlib.contextmanager
def temp_dir():
    try:
        td = tempfile.mkdtemp()
        yield td
    except:
        shutil.rmtree(td, ignore_errors=True)
        raise
    else:
        shutil.rmtree(td, ignore_errors=True)


def read_count(fastq):
    """Very fast read counter (assuming the FASTQ is formatted properly).

    Args:
        fastq (str): fastq file path

    Returns:
        int
    """
    return int(subprocess.check_output("awk '{n++}END{print n/4}' " + fastq, shell=True).decode())


def validate_group_fastqs(group_files, demultiplexing_stats, metadata_df, sample_id, group_id):
    """Validate the counts of the group total versus the individual sample files of the groups.
    Columns 'Id' and 'Count' come from `fastq-multx` STDOUT contained in demultiplexing_stats.

    Args:
        group_files (dict): group id: file path for each group FASTQ
        demultiplexing_stats (str): file path to demultiplexing stats
        metadata_df (pandas.DataFrame): dataframe of metadata containing group and sample id cols
        sample_id (str or int): sample id column in metadata_df
        group_id (str or int): group id column in metadata_df
    """
    logging.info("Validating group read counts with sample counts")
    stats = pd.read_table(demultiplexing_stats, index_col=False, usecols=[0, 1], header=0)
    merged = metadata_df.merge(stats, how="left", left_on=sample_id, right_on="Id")
    sums = merged.groupby([group_id])["Count"].sum()
    for gid, filepath in group_files.items():
        if not read_count(filepath) == sums[gid]:
            logging.warning(("{group} did not validate; read counts across samples do not add up "
                             "to the group total").format(group=gid))


if __name__ == "__main__":
    group_demux()
