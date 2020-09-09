#!/usr/bin/env /scratch/Shares/layer/workspace/michael_sandbox/Python-3.8.5/python
import re
import pandas as pd
import sqlite3
import numpy as np
import argparse
import glob
import gzip
from scipy import stats

class DatabaseAccess:
    def __init__(self, db_path):
        self.connection = sqlite3.connect(db_path)
        create_chromosome_table_template = """ CREATE TABLE IF NOT EXISTS {table}
        (
        id INTEGER PRIMARY KEY,
        position INTEGER,
        genotype TEXT,
        numerical_genotype TEXT,
        alleles TEXT,
        allele_counts TEXT,
        allele_balances TEXT,
        first_alleles TEXT,
        samples TEXT,
        reference_allele TEXT
        );
        """
        # create a table for each chromosome. the names for tables with be chromosome_1, chromosome_x etc...
        chromosomes = [str(x) for x in range(1, 24)]
        chromosomes.append('x')
        chromosomes.append('y')
        curs = self.connection.cursor()
        for c in chromosomes:
            statement = create_chromosome_table_template.format(table='chromosome_' + c)
            log_file.write(statement + '\n')
            a = curs.execute(statement)

    def print_tables(self):
        curs = self.connection.cursor()
        curs.execute("SELECT name FROM sqlite_master WHERE type='table';")
        print(curs.fetchall())

    def update_db_with_tsv(self, path_to_tsv: str):
        """
        Given a tsv, add its information into the database
        :param path_to_tsv: path to a tsv (typically the output of analyze_sample.py) with these columns:
            chromosome	position	genotype	diff	allele_balance	geno_string	number_of_alleles
        :return: None
        """
        # get the first 2 lines, the second contains the sample name
        sample = None
        with open(path_to_tsv, 'r') as file:
            for i, line in enumerate(file):
                if i == 1:
                    sample = line.split('\t')[1].strip()
                    break

        df = pd.read_csv(path_to_tsv, sep='\t', skiprows=2)
        get_template = """SELECT allele_balances, allele_counts, alleles, first_alleles, samples 
        FROM chromosome_{chrom} WHERE position = {pos} AND genotype = "{genotype}"
        """
        update_template = """ UPDATE chromosome_{chrom}
        SET 
        allele_balances = "{ab}",
        allele_counts = "{ac}",
        alleles = "{alleles}",
        first_alleles = "{ma}",
        samples = "{samps}"
        WHERE
        position = {pos}
        AND
        genotype = "{genotype}"
        """
        insert_template = """INSERT INTO chromosome_{chrom} 
        (position, allele_balances, allele_counts, genotype, first_alleles, alleles, samples, numerical_genotype)
        VALUES ({pos},"{ab}","{ac}","{genotype}","{ma}","{alleles}","{samps}","{ng}")
        """
        # for each chromosome
        for c in set(df['chromosome']):
            if str(c) != '1':
                continue
            else:
                print('Chromosome: ' + str(c))
            sub = df[df['chromosome'] == c]
            sub.index = sub['position']
            # for each position
            curs = self.connection.cursor()
            for p in sub['position']:
                allele_balance = sub.loc[p]['allele_balance']
                first_alleles = sub.loc[p]['first_allele']
                allele_count = sub.loc[p]['allele_counts']
                genotype = sub.loc[p]['genotype']
                alleles = sub.loc[p]['alleles']
                numerical_genotype = sub.loc[p]['numerical_genotype']
                curs.execute(get_template.format(chrom=c, pos=p, genotype=genotype))
                res = curs.fetchall()
                if len(res) == 0:
                    # there is no information at this position and genotype
                    # use the insert template
                    statement = insert_template.format(chrom=c, pos=p, ab=str(allele_balance),
                                                       ac=str(allele_count), genotype=genotype, ma=first_alleles,
                                                       alleles=alleles, samps=sample, ng=numerical_genotype)
                    curs.execute(statement)
                else:
                    # use the update template
                    statement = update_template.format(chrom=c,
                                                       pos=p,
                                                       ab=res[0][0] + ';' + str(allele_balance),
                                                       ac=res[0][1] + ';' + str(allele_count),
                                                       genotype=genotype,
                                                       ma=res[0][2] + ';' + first_alleles,
                                                       alleles=res[0][3] + ';' + alleles,
                                                       samps=res[0][4] + ';' + sample)
                    curs.execute(statement)
            # commit changes after each chromosome
            self.connection.commit()

    def get_mean_and_std_dev(self, chromosome: str, position: int, genotype: str) -> [float, float]:
        """
        Get the mean and standard deviation of allele balance at a certain position
        :param chromosome: string, single character chromosome '1','2' ... 'x', 'y'
        :param position: integer, numerical position on the chromosome
        :param genotype: string, genotype of interest
        :return: list, [allele balance mean, allele balance standard deviation] or
                 [-1,-1] if there is no information about the desired location
        """
        get_template = """SELECT allele_balances, allele_counts, alleles, first_alleles 
        FROM chromosome_{chrom} WHERE position = {pos} AND genotype = "{genotype}"
        """
        curs = self.connection.cursor()
        statement = get_template.format(chrom=chromosome, pos=position, genotype=genotype)
        curs.execute(statement)
        res = curs.fetchall()
        # if there are no results, return [-1, -1]
        if len(res) == 0:
            return [-1, -1]
        # allele balance is the first item returned from the query
        # here I assume there is only 1 entry that matches the query
        ab = res[0][0]
        # split up the string into individual allele balances (';' denotes separation between samples)
        allele_balances = [float(x) for x in ab.split(';')]
        mean = sum(allele_balances) / len(allele_balances)
        std = np.std(allele_balances)
        return [mean, std]

    def query_db(self, samples: list):
        """
        Give an list of samples names, calculate the difference between mean allele balances of the rest of the
        population and the requested samples and calculate a p value of the allele balances of the samples vs rest of
        the population at each site the samples are present at. Creates a file in the current working directory named
        'query_results.tsv'.
        :param samples: list of names of samples that should be compared to the results of the population in the db
        :return: pandas DataFrame with the following columns: chromosome, position, genotype, numerical_genotype,
                 population_size, p_value, avg_allele_balance_difference
        """
        output_dict = {'chromosome': [],
                       'position': [],
                       'genotype': [],
                       'numerical_genotype': [],
                       'population_size': [],
                       'p_value': [],
                       'avg_allele_balance_difference': []}
        chromosomes = [str(x) for x in range(1, 24)]
        chromosomes.append('x')
        chromosomes.append('y')

        curs = self.connection.cursor()
        for c in chromosomes:
            # # remove this after development
            # if c != '1':
            #     return
            get_template = """SELECT position, genotype, allele_balances, numerical_genotype, samples FROM chromosome_{chr} WHERE samples LIKE '%{sample}%' """
            other_sample_template = "AND samples LIKE '%{sample}%'"
            statement = get_template.format(sample=samples[0], chr=c)
            for samp in samples[1:]:
                statement += other_sample_template.format(sample=samp)
            curs.execute(statement)
            res = curs.fetchall()
            # for each row in the db with those samples listed
            for r in res:
                all_samp_names = r[4].split(';')
                all_samp_names = [x.replace('\n','') for x in all_samp_names]
                all_abs = r[2].split(';')
                all_abs = [float(x) for x in all_abs]
                in_samples = [x for i, x in enumerate(all_abs) if all_samp_names[i] in samples]
                # get those in the population
                population = [x for i, x in enumerate(all_abs) if all_samp_names[i] not in samples]
                p_val = stats.ttest_ind(in_samples,population).pvalue
                diff = (sum(population)/len(population)) - (sum(in_samples)/len(in_samples))
                output_dict['chromosome'].append(c)
                output_dict['position'].append(r[0])
                output_dict['genotype'].append(r[1])
                output_dict['numerical_genotype'].append(r[2])
                output_dict['population_size'].append(len(population))
                output_dict['p_value'].append(p_val)
                output_dict['avg_allele_balance_difference'].append(diff)
        df = pd.DataFrame(output_dict)
        df.to_csv('query_results.tsv', sep='\t')
        return df


def collect_and_output_genotypes(input_file, output_file, db_path, coverage_threshold, geno_dict, numerical_geno_dict):
    ab_info = {'chromosome': [],
               'position': [],
               'genotype': [],
               'alleles': [],
               'allele_counts': [],
               'allele_balance': [],
               'num_alleles': [],
               'first_allele': [],
               'reference_allele': [],
               'z_score': [],
               'numerical_genotype': [],
               'strand_bias': []}
    da = DatabaseAccess(db_path)
    count = 0
    geno_count = 0
    geno_count_unknown = 0
    not_enough_count = 0
    with gzip.open(input_file, 'rt') as file:
        for line in file:
            count += 1
            row = line.split('\t')

            ref_allele = row[2]
            seq = row[4]
            # if there if not enough coverage, skip it
            if len(seq) < coverage_threshold:
                continue

            # remove indel information about bases to follow
            seq = re.sub('\\+[0-9]+[ACGTNacgtn]', '', seq)
            seq = re.sub('-[0-9]+[ACGTNacgtn]', '', seq)

            # remove > and < as these represent skips to the reference genome
            seq = seq.replace('>', '').replace('<', '')
            # remove * which denote gap in the read
            seq = seq.replace('*', '')
            # remove $ which denote end of a read
            seq = seq.replace('$', '')
            # remove beginning of read marker and the quality character that follows
            seq = re.sub('\\^.', '', seq)

            # what portion of reads come from the forward stand?
            if len(seq) == 0:
                log_file.write('Error line has no sequence length: ' + line)
                continue
            strand_bias = sum(x == '.' or x.isupper() for x in seq) / len(seq)

            seq = seq.upper()
            seq = seq.replace(',', ref_allele)
            seq = seq.replace('.', ref_allele)

            # count of the occurrence of each allele
            seq_list = [char for char in seq]
            alleles = list(set(seq_list))
            allele_counts = [seq_list.count(x) for x in alleles]

            # sort alleles alphabetically
            alleles, allele_counts = zip(*sorted(zip(alleles, allele_counts)))

            # get the genotype
            genotype = None
            if str(row[0]) not in geno_dict:
                log_file.write('Error chromosome not found in  vcf: ' + str(row[0]))
                continue
            if str(row[1]) in geno_dict[str(row[0])]:
                genotype = geno_dict[str(row[0])][str(row[1])]
                geno_count += 1
            else:
                genotype = row[2] + '/' + row[2]
                geno_count_unknown += 1

            # calc AB as mode count / all count
            first_allele = genotype.split('/')[0]
            try:
                index = alleles.index(first_allele)
            except ValueError:
                log_file.write('Error: genotype not found in alleles\n')
                log_file.write(str(first_allele) + '\n')
                log_file.write(str(row) + '\n\n')
                log_file.write('Error: genotype not found in alleles\n')

            if genotype == 'ERROR':
                log_file.write('Error: genotype from VCF is ERROR: ' + line + '\n')
                continue

            if index >= len(allele_counts):
                log_file.write('Error: The genotype called is not present in the alleles ' + line + '\n')
                continue
            # calculate allele balance
            ab = allele_counts[index] / sum(allele_counts)

            # calc score for each site
            mean_std_arr = da.get_mean_and_std_dev(row[0], int(row[1]), genotype)
            z_score = -1
            if -1 in mean_std_arr:
                not_enough_count += 1
                log_file.write('Not enough info in DB to calculate Z-score\n')
            else:
                # caluclate the z-score, z = (x - mean) / standard deviation
                z_score = (ab - mean_std_arr[0]) / mean_std_arr[1]

            ab_info['chromosome'].append(row[0])
            ab_info['position'].append(row[1])
            ab_info['reference_allele'].append(row[3])
            ab_info['allele_balance'].append(ab)
            ab_info['num_alleles'].append(len(alleles))
            ab_info['first_allele'].append(alleles[index])
            # genotype will be formatted like this: 'A/C' or 'T/T'
            ab_info['genotype'].append(genotype)
            ab_info['alleles'].append('/'.join(alleles))
            ab_info['allele_counts'].append(','.join([str(x) for x in allele_counts]))
            ab_info['z_score'].append(z_score)
            ab_info['strand_bias'].append(strand_bias)
            try:
                ab_info['numerical_genotype'].append(numerical_geno_dict[str(row[0])][str(row[1])])
            except KeyError:
                log_file.write('Error: no numerical genotype found for chromosome ' + str(row[0]) + ' postition ' +
                               str(row[1]) + '\n')
                ab_info['numerical_genotype'].append('Error')

    df = pd.DataFrame(ab_info)

    out_of_range_count = sum(df['z_score'] >= 3) + sum(df['z_score'] <= -3)
    with open(output_file, 'w') as file:
        file.write('# number of sites out of range:\t' + str(out_of_range_count))
    df.to_csv(output_file, mode='a', header=True)

    return df


def get_genotype_dict_from_vcf(vcf: str) -> dict:
    """
    Create a dictionary of the genotype in the given vcf file. Dictionary is formated as:
    a dictionary of dictionaries with chromosome (string) as the first key and position as the second
    The genotype is sorted alphabetically for example 'A/T' or 'C/G' and NEVER 'T/A' or 'G/C'
    :param vcf: path to vcf file
    :return: dictionary of genotypes
    """
    # dictionary of dictionaries with chromosome as the first key and position as the second
    geno_dict = {}
    with gzip.open(vcf, 'rt') as file:
        for line in file:
            # skip headers
            if line[:1] == '#':
                continue
            row = line.split('\t')
            log_file.write(str(row) + '\n')
            ref = row[3]
            alt = row[4]
            chrom = row[0]
            pos = row[1]
            numerical_genotype = row[9].split(':')[0]
            numerical_genotype = numerical_genotype.replace('|', '/')
            genotype = ''
            if numerical_genotype == '1/1':
                genotype = alt[0] + '/' + alt[0]
            elif numerical_genotype == '0/1':
                genotype = ref[0] + '/' + alt[0]
            elif numerical_genotype == '1/2':
                alleles = row[4].split(',')
                alleles = [x[0] for x in alleles]
                # reverse sorting so they get joined in the right order
                alleles.sort(reverse=True)
                genotype = '/'.join(alleles)
            elif numerical_genotype == '0/0':
                genotype = ref[0] + '/' + ref[0]
            else:
                log_file.write('Unknown genotype in get_genotype_dict_from_vcf:\t')
                log_file.write(str(numerical_genotype) + '\n')
                genotype = 'ERROR'
            # if a the chromosome has not been seen before, create a new entry for it
            if chrom not in geno_dict:
                geno_dict[chrom] = {}
            # add the genotype in for the chromosome and position
            geno_dict[chrom][pos] = genotype
    return geno_dict


def get_numerical_genotype_dict_from_vcf(vcf: str) -> dict:
    """
    Create a dictionary of the numerical genotypes in the given vcf file. Dictionary is formated as:
    a dictionary of dictionaries with chromosome (string) as the first key and position as the second
    :param vcf: path to vcf file
    :return: dictionary of genotypes
    """
    # dictionary of dictionaries with chromosome as the first key and position as the second
    geno_dict = {}
    with gzip.open(vcf, 'rt') as file:
        for line in file:
            # skip headers
            if line[:1] == '#':
                continue
            row = line.split('\t')
            log_file.write(str(row) + '\n')
            ref = row[3]
            alt = row[4]
            chrom = row[0]
            pos = row[1]
            numerical_genotype = row[9].split(':')[0]
            numerical_genotype = numerical_genotype.replace('|', '/')
            # if a the chromosome has not been seen before, create a new entry for it
            if chrom not in geno_dict:
                geno_dict[chrom] = {}
            # add the genotype in for the chromosome and position
            geno_dict[chrom][pos] = numerical_genotype
    return geno_dict


def make_fingerprint_report(df, fingerprint, ngd):
    """
    Give a dataframe witb allele balance information, a vcf and .bed file: creates an allele balance report for the
    sites listed in the .bed file :param df: DataFrame formatted like the output of collect_and_output_genotypes
    :param vcf: path to vcf.gz file :param fingerprint: path to .bed file :return: None
    """

    high_quality_report = {'chromosome': [], 'position': [], 'genotype': [], 'numerical_genotype': [],
                           'allele_balance': [], 'z_score': []}
    current_chr = ''
    df['chromosome'] = [str(x) for x in df['chromosome']]
    for line in open(fingerprint, 'r'):
        row = line.split()
        chr = str(row[0].replace('chr', ''))
        pos = int(row[1])
        if current_chr != chr:
            sub = df[df['chromosome'] == chr]
            current_chr = chr
        subsub = sub[sub['position'] == pos]
        if len(subsub) != 0 and str(pos) in ngd[chr]:
            high_quality_report['chromosome'].append(subsub['chromosome'])
            high_quality_report['position'].append(subsub['position'])
            high_quality_report['genotype'].append(subsub['genotype'])
            high_quality_report['numerical_genotype'].append(ngd[chr][str(pos)])
            high_quality_report['allele_balance'].append(subsub['allele_balance'])
            high_quality_report['z_score'].append(subsub['z_score'])
        else:
            high_quality_report['chromosome'].append(subsub['chromosome'])
            high_quality_report['position'].append(subsub['position'])
            high_quality_report['genotype'].append('NOT FOUND')
            high_quality_report['numerical_genotype'].append('NOT FOUND')
            high_quality_report['allele_balance'].append('NOT FOUND')
            high_quality_report['z_score'].append('NOT FOUND')
    df = pd.DataFrame(high_quality_report)
    df.to_csv('fingerprint.tsv', sep='\t')
    return df


def analyze(db_path, pileup_path, vcf_path, sample, fingerprint=None):
    """
    :param db_path: path to database
    :param pileup_path: path to .pileup.gz file
    :param vcf_path: path to .vcf.gz file
    :param sample: name of the sample being analyzed
    :param fingerprint: (optional) path to fingerprinting .bed file
    :return:
            1. (if no fingerprint file is specified) pandas DataFrame with allele balance information
            2. (if fingerprint file is provided) 2 pandas DataFrames. One as previously described and one containing the
               fingerprint profiles
    """
    pileup_file = pileup_path
    vcf_file = vcf_path
    output_file = 'output.tsv'

    gd = get_genotype_dict_from_vcf(vcf_file)
    ngd = get_numerical_genotype_dict_from_vcf(vcf_file)
    df = collect_and_output_genotypes(pileup_file, output_file, db_path, 40, gd, ngd)

    out_of_range_count = sum(df['z_score'] >= 3) + sum(df['z_score'] <= -3)
    with open(output_file, 'w') as file:
        file.write('# number of samples out of range:\t' + str(out_of_range_count) + '\n')
        file.write('# sample name:\t' + sample + '\n')
    df.to_csv(output_file, mode='a', header=True, sep='\t')

    # if no finger printing .bed file is provided go ahead and return
    if fingerprint is None:
        return df
    else:
        fp_df = make_fingerprint_report(df, fingerprint, ngd)
        return df, fp_df


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--name",
                        dest="name",
                        required=False,
                        help="sample name to be used in output files required if using `--type analyze`")

    parser.add_argument("--db",
                        dest="db_path",
                        required=True,
                        help="path to exome or genome sqlite db")

    parser.add_argument("--pileup",
                        dest="pileup",
                        required=False,
                        help="path pileup file")

    parser.add_argument("--vcf",
                        dest="vcf",
                        required=False,
                        help="path vcf file")

    parser.add_argument("--fingerprint",
                        dest="fingerprint",
                        required=False,
                        help="path to a .bed file of locations to be reported in the fingerprint report ("
                             "fingerprint.tsv). This parameter is only used if using `--type analyze` and is optional")

    parser.add_argument("--type",
                        dest="run_type",
                        required=True,
                        help="function to be run: analyze sample `--type analyze`, update db `--type update` or "
                             "query information from db `--type query`")

    parser.add_argument("--inputs",
                        dest="inputs",
                        required=False,
                        help="comma separated list of file paths for `output.tsv` files generated by `--type analyze` "
                             "to add to the data base. This option is only used when using `--type update`")

    parser.add_argument("--force_update",
                        dest="force_update",
                        required=False,
                        help="If this parameter is set to `--force_update true` and you are using `--type update` the "
                             "it will force the db to update it self with the given input files even if they do not "
                             "meet the threshold requirements")

    parser.add_argument("--samples",
                        dest="samples",
                        required=False,
                        help="comma separated list of sample names to be pull allele balance information from the "
                             "database on. Only use if using `--type query`")

    args = parser.parse_args()

    return args


def check_update(files):
    """
    :param files: list of files to check if they meet the requirements to update the db
    :return: true or false (true if the db should be update with these files)
    """
    files_over_thresh = 0
    for file in files:
        # the first line of each file has the number of sites out of range in each file
        for line in open(file, 'r'):
            out_of_range_count = line.split('\t')[1]
            if int(out_of_range_count) > 0:
                files_over_thresh += 1
            break
    return files_over_thresh == 0


def update(db_path, files):
    """
    :param db_path: path to the database
    :param files: list of files to update database with
    """
    da = DatabaseAccess(db_path)
    for file in files:
        print(file)
        da.update_db_with_tsv(file)


def query(db_path, samples):
    """
    :param db_path: path to the database
    :param samples: list of samples to compare to rest of population in db
    :return: pandas DataFrame
    """
    da = DatabaseAccess(db_path)
    return da.query_db(samples)


with open('allele_balance_log.txt', 'a') as log_file:
    if __name__ == "__main__":
        args = vars(get_args())
        if args['run_type'] == 'update':
            da = DatabaseAccess(args['db_path'])
            if args['inputs'] is None:
                log_file.write('Error, no files listed with parameter `--inputs`. See documentation at '
                               'https://github.com/MSBradshaw/AlleleBalance\n')
                quit()
            files = args['inputs'].split(',')
            files = [x for x in files if x != '']
            if check_update(files) or args['force_update'] == 'true':
                update(args['db_path'], files)
            else:
                print('Not updating database, too many samples with an abnormal number of imbalanced alleles. Bad '
                      'batch suspected\n')
                log_file.write('Not updating database, too many samples with an abnormal number of imbalanced '
                               'alleles. Bad batch suspected\n')
        elif args['run_type'] == 'analyze':
            analyze(args['db_path'], args['pileup'], args['vcf'], args['name'], args['fingerprint'])
        elif args['run_type'] == 'query':
            if args['samples'] is None:
                log_file.write('Error: `--samples` parameter is required if using `--type query`')
                quit()
            samps = args['samples'].split(',')
            if len(samps) < 2:
                log_file.write('Error: `--samples` must include at least 2 comma separated items when using `--type '
                               'query`')
                quit()
            query(args['db_path'], samps)
        else:
            print('Invalid `--type`. To analyze a sample use `--type analyze` or to update the db `--type update`')
            log_file.write(
                'Invalid `--type`. To analyze a sample use `--type analyze` or to update the db `--type update`\n')
            quit()
