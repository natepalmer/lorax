# Analyze nanopore combinatorial library data
#
# main.py
# Iterate over paf file
# Filter based on coverage, length, etc
# Use cs string to create sequence
# Classify reference/alternate at each mutation position
# Build count table of library elements
# Calculate summary statistics
#
# To use this script:
# 1. Change the filenames to your target files
#       reference file (in .fasta format)
#       alignment file (in .paf format (see minimap2))
# 2. Change output filename to desired location
#       optionally, change mutation file name
# 3. Select library specification (1 or 2) or create custom
# 4. Run script

from itertools import product as iter_product
from Bio import SeqIO
from random import randint
import re


class PafRecord:
    __slots__ = ["name", "length", "start", "end", "strand", "ref_name", "ref_length", "ref_start", "ref_end",
                 "len_matches", "len_aligned", "mapping_qual", "tags", "n_mutations", "lib_identity", "indel"]

    def __init__(self, text, norm=True):
        fields = text.split(sep="\t")
        if len(fields) < 12:
            raise Exception("Input string is missing required components; see pairwise mapping format specification")
        self.name = fields[0]
        self.length = int(fields[1])
        self.start = int(fields[2])
        self.end = int(fields[3])
        self.strand = fields[4]
        self.ref_name = fields[5]
        self.ref_length = int(fields[6])
        self.ref_start = int(fields[7])
        self.ref_end = int(fields[8])
        self.len_matches = int(fields[9])
        self.len_aligned = int(fields[10])
        self.mapping_qual = int(fields[11])
        self.tags = {x.partition(":")[0]: x.partition(":")[2].partition(":")[2] for x in fields[12:]}
        self.n_mutations = None
        self.lib_identity = None
        self.indel = False
        if norm:
            self.normalize_position()

    def add_alignment_info(self, lib_identity=None, n_mutations=None, replace=False):
        """
        Add or update information about library membership and/or mutation number.

        *lib_identity* is a set containing the library member identifier strings
        *n_mutations* is an integer specifying the number of additional mutations to add to the count
        """
        if n_mutations is not None:
            if self.n_mutations is None or replace:
                self.n_mutations = n_mutations
            else:
                n_mutations += n_mutations
        if lib_identity is not None:
            if self.lib_identity is None or replace:
                self.lib_identity = lib_identity
            else:
                self.lib_identity.update(lib_identity)

    def normalize_position(self):
        if self.ref_start > 0:
            cs_chunks = chunk_paf(self)
            if len(cs_chunks) == 0:
                return
            elif not cs_chunks[0].startswith(":"):
                self.tags["cs"] = ":{}{}".format(str(self.ref_start), self.tags["cs"])
            else:
                cs_chunks[0] = ":{}".format(str(self.ref_start+int(cs_chunks[0][1:])))
                self.tags["cs"] = "".join(cs_chunks)
        return

    def to_csv(self):
        return [self.name, self.length, self.start, self.end, self.ref_start, self.ref_end, self.len_matches,
                self.len_aligned, self.tags["cs"], self.n_mutations,
                str(self.lib_identity).replace("{", "").replace("'", "").replace("}", "")]

    def __len__(self):
        return self.length

    def __str__(self):
        return "PAF record {0} aligned to {1} at pos {2} with quality {3}".format(
            self.name, self.ref_name, str(self.ref_start), str(self.mapping_qual))

    def __repr__(self):
        return "_".join([self.name, self.tags["cs"]])


def chunk_paf(paf_record):
    """ Takes a PafRecord object or a full cs string and returns a list of fields from its cs string """
    if type(paf_record) == str:
        cs = paf_record
    else:
        cs = paf_record.tags["cs"]
    separable_cs = re.sub(r'(?<!\A)([-:*+])', r',\1', cs)
    return separable_cs.split(",")


def build_gapped_alignment(reference, paf_record):
    """Rebuild a gapped alignment from paf_record with coordinates corresponding to reference"""
    reference = reference.upper()
    alignment = "-"*paf_record.ref_start
    current_position = paf_record.ref_start
    cs_fields = chunk_paf(paf_record)
    for field in cs_fields:
        if field[0] == ":":
            match_length = int(field[1:])
            alignment += reference[current_position:current_position+match_length]
            current_position += match_length
        elif field[0] == "*":
            assert reference[current_position] == field[1].upper()
            current_position += 1
            alignment += field[2].upper()
        elif field[0] == "-":
            deletion_length = len(field[1:])
            current_position += deletion_length
            alignment += "-"*deletion_length
        elif field == "+":
            pass
    return alignment


def retrieve_codons(reference, paf_record, site_list):
    """Given a reference sequence and a paf_record object with cs tag,
    return a list of codons in the read at each codon position in site_list specified by DNA position"""
    alignment = build_gapped_alignment(reference, paf_record)
    codon_list = [alignment[x:x+3] for x in site_list]
    return codon_list


def process_codons(ref_fn, paf_fn, site_list, output_fn, report=None):
    reference = str(SeqIO.read(ref_fn, "fasta").seq)
    with open(paf_fn) as paf_file:
        with open(output_fn, mode='w') as outfile:
            for i, record in enumerate(paf_file):
                if report is not None:
                    if i % report == 0:
                        print(f"Processed {i} reads")
                read = PafRecord(record.strip(), norm=False)
                codons = retrieve_codons(reference, read, site_list)
                print(read.name, file=outfile, end=",")
                print(",".join(codons), file=outfile)
    return


def calc_distances(codon, mutations):
    distances = []
    for item in mutations:
        score = 0
        for i, letter in enumerate(codon):
            if letter not in "ATGC":
                pass
            elif letter == item[i]:
                score += 1
            elif letter != item[i]:
                score -= 1
        distances.append(score)
    return distances


def extra_filters(genotype, lib_spec):
    """ Include reads with one gap or one ambiguous site by coding the unknown site randomly """
    if genotype.count('9') == 1:
        idx = genotype.index('9')
        genotype[idx] = randint(0, len(lib_spec[idx]) - 1)
    if genotype.count('8') == 1:
        idx = genotype.index('8')
        genotype[idx] = randint(0, len(lib_spec[idx]) - 1)
    return genotype


def classify_mutations(codon_file, mutation_specs, loose_filt=False):
    numbers = {"total_reads": 0, "reads_assigned": 0, "gaps": 0, "ambiguous": 0}
    element_counts = {item: 0 for item in iter_product(*(range(len(i)) for i in mutation_specs))}

    with open(codon_file) as reads:
        # Remove library elements with double mutant at position 8/9
        #for key in list(element_counts.keys()):
        #    if key[7] == 1 and key[8] == 1:
        #        del element_counts[key]
        for i, read in enumerate(reads):
            genotype = list()
            fields = read.split(",")
            assert len(fields) == len(mutation_specs)+1
            name = fields[0]
            for j, codon in enumerate(fields[1:]):
                if codon in ("", "---", "\n"):
                    mutation_id = '9'  # gap
                else:
                    distances = calc_distances(codon, mutation_specs[j])
                    # Use one of the two following lines
                    if distances.count(max(distances)) == 1:  # One best match
                    #if max(distances) == 3:  # perfect match
                        mutation_id = max(range(len(distances)), key=distances.__getitem__)
                    else:  # Multiple equal scores
                        mutation_id = '8'  # ambiguous
                genotype.append(mutation_id)
            if loose_filt:    # Include reads with one ambiguous or gapped site
                genotype = extra_filters(genotype, mutation_specs)
            try:
                element_counts[tuple(genotype)] += 1
                numbers["reads_assigned"] += 1
            except KeyError:
                if '9' in genotype:
                    numbers["gaps"] += 1
                elif '8' in genotype:
                    numbers["ambiguous"] += 1
                else:
                    # Read not able to be classified
                    print(genotype)
                    raise Exception
            numbers["total_reads"] += 1
    return element_counts, numbers


# Input data files
reference_fn = "sample_data/inputs/C9-r2-ref.fasta"
paf_fn = "sample_data/inputs/Cas9-r2-sample.paf"

# Internally used file with reads and codons - useful for data inspection
mutation_fn = "sample_data/outputs/Cas9-r2-sample_mutations.csv"

# Output file with library element counts
count_fn = "sample_data/outputs/Cas9-r2-sample_counts.csv"

# Loose filtering (keep reads with one ambiguous or gapped site)
loose_filtering = True

# Library specification to use (1 or 2)
lib_spec = 2

Cas9_mutation_positions_r1 = [82, 709, 856, 952, 1102, 1492, 1540, 1846, 1867, 1906,
                           2110, 2179, 2446, 3046, 3733, 3817, 3844, 3880]
Cas9_mutation_positions_r2 = [37,199,481,730,952,973,1333,1471,1903,2089,2419,2635,
                           3046,3115,3283,3415,3436,3658,3859]
cas9_mutations_r1 = [("CCG", "CTC"),
                  ("CTG", "TGC"),
                  ("TAC", "CAG"),
                  ("AGC", "CAT", "TGT"),
                  ("AGC", "TGT"),
                  ("TTT", "ACC"),
                  ("CTG", "ACC", "GGA"),
                  ("CTT", "GGA"),
                  ("CTT", "CAG"),
                  ("TTG", "GAC"),
                  ("TTT", "GCC"),
                  ("CTT", "CCC", "GGA"),
                  ("CTG", "GAC"),
                  ("TAC", "AAG", "GGA"),
                  ("CTC", "GGA"),
                  ("ATA", "CAG"),
                  ("CTC", "GCT", "GAG"),
                  ("TAC", "CAG")
                  ]
cas9_mutations_r2 = [("ACC", "AAT", "GGA"),
                  ("ACC", "CGG"),
                  ("ATG", "CCC"),
                  ("CTG", "GCC", "GGA"),
                  ("AGC", "TGT"),
                  ("TAC", "GAT"),
                  ("ACC", "GAG"),
                  ("TTC", "AAT", "ACA"),
                  ("CGG", "GAA"),
                  ("ATC", "GAC"),
                  ("CAG", "AGA"),
                  ("ATG", "GAC", "CGC"),
                  ("TAC", "GAA"),
                  ("TAC", "GAA"),
                  ("GTG", "AGC"),
                  ("GTG", "AAA"),
                  ("GTG", "GAA"),
                  ("CTG", "AAT"),
                  ("CTG", "AAA")]

if lib_spec == 1:
    cas9_site_list = [x-1 for x in Cas9_mutation_positions_r1]
elif lib_spec == 2:
    cas9_site_list = [x-1 for x in Cas9_mutation_positions_r2]
else:
    raise Exception("Incorrect library specification. Use 1 or 2 or create custom spec")

# This produces a file containing each read (1 per line) and the codon aligned to each mutation site
process_codons(reference_fn, paf_fn, cas9_site_list, mutation_fn, report=10000)

# This will classify the reads created by process_codons into a table of element counts and a list of summary stats
element_counts, numbers = classify_mutations(mutation_fn, cas9_mutations_r2, loose_filt=loose_filtering)

# This will print the data in element counts to file
with open(count_fn, mode="w") as count_file:
    for item in element_counts:
        print("".join([str(x) for x in item])+","+str(element_counts[item]), file=count_file)
print(numbers)
