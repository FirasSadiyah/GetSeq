#! /usr/bin/env python3

import argparse
import json
import logging
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
import requests


class EnsemblRestClient():

    def __init__(self, server="https://rest.ensembl.org", reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs, params, regions):
        """ construct the POST or GET request"""
        if params:
            endpoint += '?' + urllib.parse.urlencode(params)
        else:
            endpoint += '?'
        data = None
        # check if rate limit is needed
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        # submit the POST or GET request to Ensembl REST API server
        # and try to catch errors returned by the server

        if regions:
            request = requests.post(
                self.server + endpoint, headers=hdrs, data=json.dumps(regions))
        else:
            request = requests.get(self.server + endpoint, headers=hdrs)

        if not request.ok:
            request.raise_for_status()
            sys.exit()

        data = request.json()
        self.req_count += 1

        return data

    def get_genomes(self):
        """ list available genomes from Ensembl"""

        # construct the GET request
        hdrs = {"Content-Type": "application/json"}
        params = {}
        endpoint = '/info/species'

        # process the GET request
        releases = self.perform_rest_action(
            endpoint, hdrs, params, regions=None)

        # list available genomes sorted by name
        if releases:
            sorted_releases = sorted(releases['species'], key=lambda k: k['name'])
            print("%-40s %s" % ('Name', 'Display name'))
            print('-' * 30, ' ' * 9, '-' * 30)
            for i in sorted_releases:
                print("%-40s %s" % (i['name'], i['display_name']))

    def get_assemblies(self, species):
        """ list available assemblies of specific species from Ensembl """

        # construct the GET request
        hdrs = {"Content-Type": "application/json"}
        params = {}
        endpoint = '/info/assembly/{0}'.format(species)

        # process the GET request
        genomes = self.perform_rest_action(
            endpoint, hdrs, params, regions=None)
        if genomes:
            genomes_list = genomes['coord_system_versions']
            for item in genomes_list:
                print(item)

    def get_sequences(self, species, assembly, regions, upstream, downstream):
        """ return the DNA sequences for the submitted regions """

        # construct the POST request
        hdrs = {'Content-Type': 'application/json',
                'Accept': 'application/json'}
        params = {'coord_system_version': assembly,
                  'expand_3prime': downstream, 'expand_5prime': upstream}
        endpoint = '/sequence/region/{0}'.format(species)

        # process the POST request
        sub_sequences = self.perform_rest_action(
            endpoint, hdrs, params, regions=regions)
        return sub_sequences

    def submit_regions(self, species, assembly, regions, upstream, downstream):
        """ iterate through the list of regions and process them individually """
        sequences = []
        for key, value in regions.items():
            sub_sequences = self.get_sequences(species, assembly, {'regions': value}, upstream, downstream)
            sequences.extend(sub_sequences)
        return sequences

    def get_regions(self, bed):
        """ prepare the regions as documented in Ensembl REST APIs """

        # checking for headers in the input bed file
        try:

            if bed:
                reads = pd.read_csv(bed, sep='\t', usecols=[
                    0, 1, 2], header=None)

            else:
                reads = pd.read_csv(sys.stdin, header=None, usecols=[0, 1, 2],
                                    delim_whitespace=True)

            start_cell = reads.iloc[0][1]

            if str(start_cell).isdigit():
                logging.info('\n No header detected, proceeding ...')
            elif start_cell in {'start', 'Start'}:
                reads = reads.iloc[1:]
                logging.info('\n Header detected, proceeding ...')
            else:
                logging.info('\n Error: Please check the columns headers of the input bed file')
                sys.exit(1)

        # raise an error if bed file is badly formatted
        except ValueError as e:
            if 'Usecols do not match columns' in str(e):
                print(
                    '\n Error: Please check the indentation of the columns in input bed file')
                sys.exit(1)

        # construct regions as defined in the API, e.g. '1:64807745..65123006'
        # respect the API limitation of 50 regions per request
        regions = {}

        reads.columns = ['chrom', 'start', 'end']

        logging.info('\n Printing top 10 input regions ... \n')
        logging.info(reads.iloc[:10, :3])

        if reads['chrom'].astype(str).str.contains('chr', na=False).any():
            # change 'chrM' to 'MT' and remove the 'chr' prefix
            reads['chrom'] = reads['chrom'].str.replace('chrM', 'MT')
            reads['chrom'] = reads['chrom'].str.replace('chr', '')
        else:
            pass

        pd.options.mode.chained_assignment = None
        for k, subreads in reads.groupby(np.arange(len(reads)) // 50):
            subreads['region'] = subreads['chrom'].astype(
                str) + ':' + subreads['start'].astype(str) + '..' + subreads['end'].astype(str)

            regions.update(
                {subreads.iloc[0]['region']: subreads['region'].tolist()})

        logging.info('\n Printing top 10 parsed regions ... \n')
        regions_list = list(regions.values())
        regions_topten = regions_list[0][:10]
        for region in regions_topten:
            logging.info(region)
        return regions


def retrieve_seq(args):
    getseq_logger(args.log)
    client = EnsemblRestClient()
    regions = client.get_regions(args.bed)
    regions_count = sum(len(v) for v in regions.values())

    logging.info(
        '\n Retrieving %s DNA sequences from Ensembl (%s/%s) ...' % (regions_count, args.species, args.assembly))

    sequences = client.submit_regions(args.species, args.assembly, regions, args.upstream, args.downstream)

    # save the sequences to output file
    if sequences:
        if args.output:
            # write results to output
            with open(args.output, "w") as f:
                for sequence in sequences:
                    print(">{id}\n{seq}".format(**sequence), file=f)
        else:
            # write results to sys.stdout
            for sequence in sequences:
                sys.stdout.write(">{id}\n{seq}\n".format(**sequence))

    logging.info('\n Done')


def retrieve_genomes(args):
    client = EnsemblRestClient()

    print('\n Retrieving a list of available genomes ... \n')
    client.get_genomes()


def retrieve_assemblies(args):
    client = EnsemblRestClient()

    print('\n Retrieving a list of available assemblies for ', args.species, '\n')
    client.get_assemblies(args.species)


def getseq_logger(log_file):
    logging.basicConfig(level=logging.INFO, format='%(message)s', filename=log_file, filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler())


def getseq():
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='getseq',
                                     description='Parse genomic regions from \
                                    a bed file and retrieve their DNA sequences')
    subparsers = parser.add_subparsers(dest='command')

    # create the parser for the "sequences" command
    parser_sequences = subparsers.add_parser('sequences')
    parser_sequences.add_argument('species', type=str,
                                  help='the species, e.g. human', metavar='species')
    parser_sequences.add_argument('assembly', type=str,
                                  help='the genome assembly, e.g. GRCh37', metavar='assembly')
    parser_sequences.add_argument('-b', '--bed', type=str,
                                  help='source bed file', metavar='')
    parser_sequences.add_argument('-o', '--output', type=str,
                                  help='output file', metavar='')
    parser_sequences.add_argument('-l', '--log', type=str, default='getseq.log',
                                  help='log file', metavar='')
    parser_sequences.add_argument('-u', '--upstream', type=int, default=0,
                                  help='the number of bp to extend upstream', metavar='')
    parser_sequences.add_argument('-d', '--downstream', type=int, default=0,
                                  help='the number of bp to extend downstream', metavar='')
    parser_sequences.set_defaults(func=retrieve_seq)

    # create the parser for the "genomes" command
    parser_genomes = subparsers.add_parser('genomes')
    parser_genomes.set_defaults(func=retrieve_genomes)

    # create the parser for the "assemblies" command
    parser_assemblies = subparsers.add_parser('assemblies')
    parser_assemblies.add_argument('species', type=str,
                                   help='the species, e.g. human', metavar='species')
    parser_assemblies.set_defaults(func=retrieve_assemblies)

    # parse the args or show help if no argument provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        args = parser.parse_args()
        args.func(args)


if __name__ == '__main__':
    getseq()
