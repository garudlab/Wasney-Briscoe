""" 
Constructs a species union file for a given host.

JCM 20220824 - Modified by MWW 20230424
"""


import sys
import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import logging
import time
import argparse
import warnings

import numpy
import pandas as pd


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""

    def error(self, message):
        """Print error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


class ConstructSpeciesUnion():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    def constructSpeciesUnionParser(self):
        """Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                'Given a list of individuals and their corresponding '
                'accession numbers, this creates the necessary species '
                ' files needed for MIDAS processing.'),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument(
            'host_srs_list', type=self.ExistingFile,
            help=('Comma separated values file containing a list of hosts'
                  ' and their associated SRA files'))
        parser.add_argument(
            'input_host', type=str,
            help=('String representing a subject id associated with '
                  'at least one SRA file.'))
        parser.add_argument(
            'midas_path', type=str,
            help=('The path to where the MIDAS directories live.'))
        parser.add_argument(
            '--min_cov', type=float,
            help=('Float describing the  minimum coverage to include in '
                  'the species unions.'), default=3.0)
        return parser

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.constructSpeciesUnionParser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        input_host = args['input_host'].rstrip()
        min_cov = args['min_cov']
        host_srs_list = args['host_srs_list']
        midas_path = args['midas_path']
        host_srs_df = pd.read_csv(host_srs_list)
        host_srs_dict = {}
        for index, row in host_srs_df.iterrows():
            subject_id = str(row['subject_id'])
            SRS = row['SRS']
            if subject_id not in host_srs_dict.keys():
                host_srs_dict[subject_id] = [SRS]
            else:
                host_srs_dict[subject_id].append(SRS)

        species_set = set()
        for file in host_srs_dict[input_host]:
            indir = midas_path + 'midas_output/' + str(file) + '/species/'
            profile_df = pd.read_csv(indir + 'species_profile.txt',
                                               sep='\t')
            # print(profile_df)
            coverage_df = profile_df.loc[profile_df['coverage'] >= min_cov]
            this_set = set(coverage_df['species_id'].to_list())
            species_set = species_set.union(this_set)
            del this_set

        species_union_list = list(species_set)
        print(species_union_list)

        for file in host_srs_dict[input_host]:
            # create output directory if needed
            outprefix = midas_path + 'species_union/' + str(file) + '/'
            outdir = os.path.dirname(outprefix)
            if outdir:
                if not os.path.isdir(outdir):
                    if os.path.isfile(outdir):
                        os.remove(outdir)
                    os.mkdir(outdir)

            # Output files: logfile
            # Remove output files if they already exist
            underscore = '' if outprefix[-1] == '/' else '_'
            output_plot = \
                '{0}{1}likelihood_surface.jpg'.format(
                    outprefix, underscore)
            species_union = \
                '{0}{1}species_union.txt'.format(
                    outprefix, underscore)
            logfile = '{0}{1}log.log'.format(outprefix, underscore)

            to_remove = [logfile, species_union]
            for f in to_remove:
                if os.path.isfile(f):
                    os.remove(f)

            # Set up to log everything to logfile.
            logging.shutdown()
            logging.captureWarnings(True)
            logging.basicConfig(
                format='%(asctime)s - %(levelname)s - %(message)s',
                level=logging.INFO)
            logger = logging.getLogger(prog)
            warning_logger = logging.getLogger("py.warnings")
            logfile_handler = logging.FileHandler(logfile)
            logger.addHandler(logfile_handler)
            warning_logger.addHandler(logfile_handler)
            formatter = logging.Formatter(
                '%(asctime)s - %(levelname)s - %(message)s')
            logfile_handler.setFormatter(formatter)
            logger.setLevel(logging.INFO)

            # print some basic information
            logger.info('Beginning execution of {0} in directory {1}\n'.format(
                prog, os.getcwd()))
            logger.info('Progress is being logged to {0}\n'.format(logfile))
            logger.info('Parsed the following arguments:\n{0}\n'.format(
                '\n'.join(['\t{0} = {1}'.format(*tup) for tup in args.items()])))
            logging.shutdown()
            output_df = pd.DataFrame()
            output_df['species_id'] = species_union_list
            output_df['count_reads'] = [''] * len(species_union_list)
            output_df['coverage'] = [''] * len(species_union_list)
            output_df['relative_abundance'] = [''] * len(species_union_list)
            output_df.to_csv(species_union, sep='\t', header=True, index=False)


if __name__ == '__main__':
    ConstructSpeciesUnion().main()
  
