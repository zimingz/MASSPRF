import subprocess
from pathlib import Path
from argparse import ArgumentParser
from functools import reduce
from multiprocessing import Pool
import itertools
import copy
import os
import gffutils
import vcf
import pandas as pd
from Bio import SeqIO, Seq, SeqRecord, AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MuscleCommandline



'''constants definitions'''
'''tool functions'''

def allele(gt):
    '''OR gate heterozygous/homozygous dominant'''
    return gt in ['0/1','1/1']


''' begin fundamental data structure definitions'''

class HomologyMap(object):

    def __init__(self, name, homologs = 'ALL'):
        homologs = homologs.upper()
        if homologs == 'ALL':
            self.homologs = 'ALL'
        else:
            pass
        self.name = name
        self.group_name_dict = {}
        self.group_genome_dict = {}
        self.group_list = []

    def addGroup(self, groupname):
        self.group_name_dict[groupname] = []
        self.group_list.append(groupname)

    def addToGroup(self, groupname, strainname):
        if not isinstance(strainname,list):
            strainName = [str(strainname)]
        if groupname not in self.group_list:
            self.addGroup(groupname)
        self.group_name_dict[groupname] += strainname

    def buildMapToGenomes(self, genomes):
        for group in self.group_list:
            self.group_genome_dict[group] = [genome for genome in genomes if genome.name in self.group_name_dict[group]]

class Genome(object):
    '''Parent class for ReferenceGenome and VariantGenome

    may change representation to abstract class, not sure how this would work
    '''
    def __init__(self, species, numchromosomes):
        self.species = species
        self.numchromosomes = numchromosomes
        self.CDS = {}
        self.genes_list = []
        self.chromosome_map = {n:'' for n in range(1, numchromosomes+1)}
        self.genes_map = {}
        self.variants = None
    def __repr__(self):
        return repr(self.species + "_" + self._name)

    def __str__(self):
        return str(self.species + "_" + self._name)

    def __iter__(self):
        for genome in self.substrains_list:
            yield self.substrains_map[genome]

    def __getitem__(self,key):
        if isinstance(key,str) and key in self.substrains_list:
            return self.substrains_map[key]
        elif isinstance(key,int) and 0 <= key < len(self.substrains_list):
            return self.substrains_map[self.substrains_list[key]]
        else:
            return self
    @property
    def name(self):
        return str(self)

    def get_chromosome(self,number):
        return self.chromosome_map[number]

    def get_all_chromosomes(self):
        for n in range(1, self.numchromosomes+1):
            yield(self.get_chromosome(n))

    def gene_from_index(self, index):
        try:
            genename = self.genes_list[index]
            return self.genes_map[genename]
        except:
            return False

    def addCDS(self, gene):
        if not isinstance(gene, CodingAnnotation):
            raise AttributeError("passed gene is not of type CodingAnnotation")
        
        if not self.genes_list:
            self.genes_list.append(gene.name)
        else: # insertion sort of genes by length
            index = 0
            curgene = self.gene_from_index(index)
            while len(gene) >= len(curgene):
                index += 1
                curgene = self.gene_from_index(index)
                if curgene == False:
                    break
            self.genes_list.insert(index, gene.name)
        self.genes_map[gene.name] = gene


    def write(self, directory, filename = None):
        if not filename:
            filename = str(self)
        print("writing " + str(self))
        if str(directory):
            #insert event logger here
            f = Path(directory).joinpath(filename)
            f.touch()
            with f.open('w') as writefile:
                SeqIO.write([SeqRecord.SeqRecord(chromo, id=chromo.id, description=chromo.description) for chromo in self.get_all_chromosomes()], writefile, "fasta")
    
    def export(self,directory, filename = None):
        for genome in self:
            genome.write(directory, filename)

class ReferenceGenome(Genome):
    '''ReferenceGenome:
        extends Genome

        attributes:
            species: a blunt grouping name useful
            numchromosomes: the total number of chromosomes in the genome
            chromosome_map: a dictionary mapping each chromosome number to a Chromosome object that holds the sequence
            substrains: dictionary mapping of VariantGenome's

        methods:
            __init__:
                initialize with no substrains, set strain = "reference"
                initialize numchromosomes
                initialize Chromosomes dictionary
            addStrain(strain):
                type(strain) = VariantGenome
                add a variant strain to the substrains dictionary
    '''
    def __init__(self, species, referenceChromosomes):
        #referenceChromosomes is a list of chromosomes w/ sequences
        Genome.__init__(self, species, max(referenceChromosomes))
        
        self._name = "reference"
        self.substrains_map = {str(self):self}
        self.substrains_list = [str(self)]
        for n in range(1, self.numchromosomes+1):
            self.chromosome_map[n] = Chromosome(self, n, referenceChromosomes[n])
    
    def __repr__(self):
        outstr = Genome.__repr__(self) + " mapped to " + repr(self.substrains_list[1:])
        return repr(outstr)

    def addStrain(self, strain):
        if not isinstance(strain, VariantGenome):
            raise AttributeError("error adding variantstrain, improper variantGenome supplied")
        self.substrains_map[str(strain)] = strain
        self.substrains_list.append(str(strain))

    def addCDS(self, gene):
        Genome.addCDS(self,gene)
        for strain in self.substrains_list:
            self.substrains_map[strain].genes_list = self.genes_list
            self.substrains_map[strain].genes_map = self.genes_map

class VariantGenome(Genome):
    '''VariantGenome
    extends: Genome
    
    attributes:
        species: a blunt grouping name useful
        numchromosomes: the total number of chromosomes in the genome
        chromosome_map: a dictionary mapping each chromosome number to a Chromosome object that holds the sequence
        strain: substrains of species grouping
        variants: a list of Variant s, each with a chromosome number and index adjusted position.

    methods:
        __init__: initialized with a ReferenceGenome object, a strain name, and a list of Variant s
            1) initializes chromosomes to ReferenceGenome Chromosomes
            2) iterates over snps (Variants) and inserts the SNPs at corresponding positions in the strings held by
                 the chromosome dictionary
    '''
    def __init__(self, reference, strain, variants):
        Genome.__init__(self, reference.species, reference.numchromosomes)
        self.reference = reference
        self._name = strain
        
        self.substrains_list = [str(self)]
        self.substrains_map = {str(self): self}
        self.relatives = []
        self.variants = variants 
        if isinstance(self.variants, dict): 
            self.chromosome_map = reference.chromosome_map
            # if variants is a list of variant
        else:
            self.variants = None
            for n in range(1, self.numchromosomes+1):
                self.chromosome_map[n] = Chromosome(self, n, variants[n])
        print(str(self) + " genome initialized")

    def get_chromosome(self, number):
        if self.variants:
            chromosome = copy.deepcopy(self.chromosome_map[number])
            chromosome.genome = self
            if self.variants.get(number):
                for snp in self.variants[number]:
                    if snp.is_polymorphic():
                        print("SNP " + snp.gt + " inserted in place of " + snp.ref + " at "
                            + str(snp.coordinate.pos) + " original sequence was " 
                            + chromosome[snp.coordinate.pos-1:snp.coordinate.pos+2])

                        chromosome[snp.coordinate.pos] = snp.gt
                        print("new sequence is " + chromosome[snp.coordinate.pos-1:snp.coordinate.pos+2])
            return chromosome

    def addRelative(self, relative):
        self.relatives.append(relative)

    def addCDS(self, gene):
        self.reference.addCDS(gene)

class Chromosome(Seq.MutableSeq):
    '''Chromosome
    extends Bio.Seq.Seq
    See biopython documentation for full documentation
    (new) attributes:
        reference: parent genome, can be a ReferenceGenome or VariantGenome
        number: chromosome number, cast as int
        sequence: str sequence of chromosome extracted 
    '''
    def __init__(self, reference, number, sequence):
        Seq.MutableSeq.__init__(self, sequence, generic_dna)
        self._genome = reference
        self.number = int(number)
        self.id = str(self.number)
        self._description = str(self._genome) + " chromosome " + str(self.id)

    @property
    def description(self):
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def genome(self):
        return self._genome

    @genome.setter
    def genome(self, genome):
        self._genome = genome
        self.description = str(self._genome) + " chromosome " + str(self.id)

class CodingAnnotation(object):

    '''PLEASE make sure to adopt 0 indexing of coordinates from gff3'''
    def __init__(self, reference, name, coordinates, strand, homolog = None, gene_id = None):
        self.reference = reference
        self.species = reference.species
        self.name = name
        self.coordinates = coordinates
        self.chromosome = coordinates[0].chromosome
        self.strand = strand
        self.homolog = homolog
        if gene_id:
            self.gene_id = gene_id
        else:
            self.gene_id = name
        self.length = len(self)

    def __repr__(self):
        return repr(str(self))

    def __str__(self):
        return str(self.species) + ' ' + str(self.name)

    def __len__(self):
        return reduce(lambda x, y: x+y, map(lambda x: len(x), self.coordinates))

    def getSequence(self, strain):
        curstrain = strain
        chromosome = curstrain.get_chromosome(self.chromosome)
        sequence = ''.join([str(chromosome[coordinate.pos[0]:coordinate.pos[1]]) for coordinate in self.coordinates])
        return CodingSequence(curstrain, self.name, self.strand, sequence, complemented = False, gene_id = self.gene_id, homolog = self.homolog)

class CodingSequence(Seq.Seq):
    """
    docstring for CodingSequence

    Created by coding annotation class

    attributes:
        self.strain - reference to parent strain genome
        self.name - name of coding sequence
        self.strand - string "+" or '-'
        self.sequence - input as a string
        self.complemented - tracks whether has been complemented
    methods:
        reverse_complement: 
            if has been complemented already, returns self;
            if has not been complemented and strand == '-',
            return instance of CodingSequence w/ reverse complement of gene
    """
    def __init__(self, strain, name, strand, sequence, complemented = False, gene_id = None, homolog = None):
        if not isinstance(strain, Genome):
            raise AttributeError("no valid strain linked to coding sequence")
        Seq.Seq.__init__(self, sequence, generic_dna)
        self.strain = strain
        self.name = name
        self.strand = strand
        self.complemented = complemented
        self.homolog = homolog

        if gene_id:
            self.gene_id = gene_id
        else:
            self.gene_id = str(self.strain) + ' ' + self.name 

        if self.strand == '-':
            self = self.reverse_complement()

    def __repr__(self):
        return self.gene_id

    def get_record(self,descriptor = None):
        desc = str(self.strain) + '_' + str(self.name)
        return SeqRecord.SeqRecord(Seq.Seq(), id = desc, description = desc)

    def reverse_complement(self):
        if not self.complemented and self.strand == '-':
            return CodingSequence(self.strain, self.name, self.strand, str(Seq.Seq.reverse_complement(self)), complemented = True)
        else:
            return self

class Variants(object):
    '''Variants, a composite class of Variant leafs
        attributes:
            strains: the strains that have been called
            locations: the Coordinate locations of each variant
            _snps_by_loc: private dictionary of dict[Coordinate] : [Variant,Variant,Variant]
            _snps_by_strain: private dictionary of dict[strain] : dict[chromosome]: [Variant, Variant, Variant]
        methods:
            initialization: provide strains or raise exception;
                initialize all attributes to empty
            addVariant: Given a variant of type Variant, 
                1) add its coordinate to self.locations
                2) add its coordinate to the keys of self._snps_by_loc
                3) append variants to that list
                4) append variants to _snps_by_strain[variant.strain]
            of_strain: 
                arguments: string(strain), string(chromosome)
                returns list of Variants of strain
            on_chromosome:
                returns hits, dict[hits] of variants on a chromosome
            in_range:
                returns variants on a chromosome within a range
    '''
    def __init__(self, strains):
        if not strains:
            raise AttributeError("No strains provided")
        self.strains = strains
        self.locations = []
        self._snps_by_loc = {}
        self._snps_by_strain = {strain:[] for strain in self.strains}

    def addVariant(self, variant):
        if variant.coordinate not in self.locations:
            self.locations.append(variant.coordinate)
        if variant.coordinate not in self._snps_by_loc.keys():
            self._snps_by_loc[variant.coordinate]=[]
        self._snps_by_loc[variant.coordinate].append(variant)
        if not self._snps_by_strain.get(variant.strain):
            self._snps_by_strain[variant.strain] = {}
        if not self._snps_by_strain[variant.strain].get(variant.chromosome):
            self._snps_by_strain[variant.strain][variant.chromosome] = []
        self._snps_by_strain[variant.strain][variant.chromosome].append(variant)

    def of_strain(self, strain, chromosome = False):
        if strain not in self.strains:
            raise ValueError("strain not in Variants")
        else:
            print(strain)
            if chromosome:
                return self._snps_by_strain[strain][chromosome]
            else:
                return self._snps_by_strain[strain]

    def on_chromosome(self, chromosome):
        hits = list(filter(lambda l: l.chromosome == str(chromosome), self.locations))
        variants = {}
        for n in hits:
            variants[n] = self._snps_by_loc[n]
        return (hits, variants)

    def in_range(self, chromosome, x1,x2):
        if x2 <= x1:
            x2, x1 = x1, x2
        hits, variants = self.on_chromosome(chromosome)
        in_range = filter(lambda x: x1 <= x.pos <= x2, hits)

        return (list(in_range), {x:variants[x] for x in in_range})
        
class Variant(object):
    '''Variant

    SNPs
    attributes:
        self.coordinate: type(coordinate) = Coordinate; 0 indexed position on 1-indexed chromosome
        self.strain: string representation of strain
        self.gt: string representation of allele

    methods:
        __init__: 
            type checking of strain, coordinate
            init coordinate, strain, gt
        __repr__:
            internal representation of Variant
            returns repr((self.strain, self.coordinate, self.gt))
    '''
    def __init__(self, coordinate, strain, gt, ref):
        if not strain:
            raise AttributeError("no strains supplied to variants")
        if not isinstance(coordinate, Coordinate):
            raise AttributeError("invalid coordinate supplied")
        self.coordinate = coordinate
        self.chromosome = coordinate.chromosome
        self.strain = strain
        self.gt = str(gt)
        self.ref = str(ref)

    def __repr__(self):
        return repr((self.strain, self.coordinate, self.gt))

    def is_polymorphic(self):
        return self.gt != self.ref

class Coordinate(object):
    '''
    Coordinate

    attributes:
        self.chromosome: integer representation of chromosome. 1 - indexed
        self.start: 0-indexed single identifier position - exchangeable for absolute position in the instance of no range
        self.stop: used if the coordinate exists over a range, as in exons/features.  otherwise, points to self.start
        self._pos: returned by @property self.pos, contains either just self.start or (self.start, self.stop)

    methods:
        __init__(chromosome, start, stop = False:
            initialize with chromosome position - casts chromosome to integer representation if not already
            initialize start/stop positions - cast to int
            default is single position
        __repr__:
            internal representation of Coordinate
            return repr((self.chromosome, self.start))
        @property
        pos:
            return self._pos, content of self._pos dependent on initialization (see above)

    '''
    def __init__(self, chromosome, start, stop = False):
        if not chromosome:
            raise AttributeError("No chromosome supplied for coordinate")
        if not start:
            raise AttributeError("No start position supplied for coordinate")
        self.chromosome = int(chromosome)

        self.start = int(start)
        if not stop:
            self.stop = start
            self._pos = self.start
        else:
            self.stop = int(stop)
            self._pos = (self.start, self.stop)

    def __len__(self):
        if self.stop != self.start:
            return self.stop-self.start
        else:
            return 1

    def __repr__(self):
        return repr((self.chromosome, self.pos))

    @property
    def pos(self):
        if isinstance(self._pos, tuple):
            self._pos = (int(self._pos[0]),int(self._pos[1]))
        else:
            self._pos = int(self._pos)
        return self._pos

'''adaptor builders'''

class CDSAdaptorBuilder(object):
    '''CDSAdaptorBuilderis a factory for coding sequence adaptors classes
        To extend functionality to new file types, write a new class of Adaptor that creates objects of class Gene and accepts **kwargs
        return the new class using an elif statement in __new__ of CDSAdaptorFactory'''
    def __new__(self,cls, **kwargs):
        if cls == 'gff3':
            return GffUtilAdaptor(**kwargs)

class VariantsAdaptorBuilder(object):
    '''Interpret passed variant format file and return properly interfaced variant file'''
    def __new__(self, variants, frmt = 'VCF'):
        frmt = frmt.upper()
        if frmt == 'VCF':
            return PyVCFAdaptor(variants)

class GenomeAdaptorBuilder(object):
    '''similar to other AdaptorBuilders, this one builds genome adaptors based on passed formatting
    add support **kwargs to add variants'''
    def __new__(self, reference, species, frmt = "FASTA", variants = None, **kwargs):
        frmt = frmt.upper()
        if not Path(reference).is_file():
            raise AttributeError("no valid reference genome supplied")
        if not species:
            raise AttributeError("No species supplied")
        if not frmt:
            raise IOError("No parse format supplied")
        if frmt == 'PLACEHOLDER':
            #add new formats for parsing here
            pass
        elif frmt == "FASTA":
            #default to FASTA format for genome
            return FastaGenomeAdaptor(reference, species, variants)
        else:
            pass

class HomologAdaptorBuilder(object):
    '''Adapts input homolog formats to HomologAdaptor'''
    def __new__(self, homologyfile, species, frmt = "CSVB"):
        frmt = frmt.upper()
        if not Path(homologyfile).is_file():
            raise AttributeError("No valid homology file specified")
        if frmt == "CSVA":
            #this file type is useful when you have two distinct divergent species with a CSV file mapping of the homologs
            return CSVaHomologs(homologyfile,species)
        elif frmt == "CSVB":
            return CSVbHomologs(homologyfile,species)

'''adaptors'''

class GffUtilAdaptor(object):
    '''to do:
    Document this class

    1) fix **kwargs checking to account for all scenarios - consider using paper logic diagram first
        - document logic as right now it is a bit obtuse
    2) fix createGene to initalize Gene classes
    3) integrate with Genome?
    '''
    def __init__(self, **kwargs):
        if not kwargs:
            raise IOError("No gffutils database or gff3 supplied")
        self.dbname = ''
        self.gff3 = ''
        if 'gff3' in kwargs:
            self.gff3 = Path(kwargs['gff3'])
            
        if 'db' in kwargs:
            if Path(kwargs['db']).is_file():
                self.dbname = Path(kwargs['db'])
        elif 'db' not in kwargs:
            self.dbname = Path(str(kwargs['gff3']) + 'db')
        if not self.dbname.is_file():
            self.db = gffutils.create_db(str(self.gff3), str(self.dbname))
        self.db = gffutils.FeatureDB(str(self.dbname))

    def getGenes(self, genes = 'ALL'):
        is_gene = lambda x: x.featuretype == 'gene'
        if genes == 'ALL':
            return (feature for feature in self.db.all_features() if is_gene(feature))
        else:
            is_specified = lambda x: x.attributes['Name'][0] in genes
            return (feature for feature in self.db.all_features() if is_gene(feature) and is_specified(feature))

    def createCodingAnnotation(self, feature, reference, **kwargs):
        if not isinstance(feature, gffutils.feature.Feature):
            raise AttributeError("passed feature is not of type gffutils.feature.Feature")
        coordinates = list(map(lambda x: Coordinate(int(x.chrom), x.start-1, x.end), self.db.children(feature, featuretype = 'CDS', order_by = 'start')))
        name = feature.attributes['Name'][0]
        gene_id = feature.attributes['ID'][0]
        strand = feature.strand
        cds = CodingAnnotation(reference, name, coordinates, strand, gene_id = gene_id)
        reference.addCDS(cds)
        return cds


## PLEASE NOTE****** GFF3 is 1-indexed, inclusive on both sides.  Therefore, to translate to python string slicing, the lower coordinate must be -1.  
## if the range on the GFF is 1-10, then the python string indices will be 0-9.  Therefore, str[0:10].  Check this to make sure it is right.
class PyVCFAdaptor(object):
    '''adapts VCF files into Variants class of Variant leafs
    skips all non-snps for now.  Fix!'''
    def __new__(self, variants):
        if not Path(variants).is_file():
            raise IOError("Invalid filename supplied for VCF")
        else:
            self.filename = variants
            self.file = open(self.filename)
            self.vcf = vcf.Reader(self.file)
            
            strains = self.vcf.samples
            variants = Variants(strains)
            for record in self.vcf:
                calls = []
                if record.is_snp:
                    coordinate = Coordinate(record.CHROM, record.POS-1)
                    for strain in variants.strains:
                        calls.append(allele(record.genotype(strain)['GT']))

                    ref = record.REF
                    alt = record.ALT[0] #See above - SNP reading breaks >1 alternate allele

                    alleles = [alt if gt else ref for gt in calls]

                    for a, s in zip(alleles, variants.strains):
                        variants.addVariant(Variant(coordinate, s, a, ref))
            self.file.close()
        return variants

class FastaGenomeAdaptor(object):
    '''FastaGenomeAdaptor adapts biopython SeqIO genome into the genome class
    "species" here is a generic term, eg rice, banana, etc, that is useful for bluntly differentiating objects of type genome

    Consider reworking the external representation of reference genome and its variants
    '''
    def __new__(self, reference, species, variants = None):
        #GenomeAdaptorBuilder already did the type checking of reference to insure it is a filetype
        reference_filename = reference
        reference_file = open(reference)

        parsedchromosomes = {int(n.id):str(n.seq) for n in SeqIO.parse(reference_file, "fasta")}
        referenceGenome = ReferenceGenome(species, parsedchromosomes)
        reference_file.close()

        if not isinstance(variants, Variants) and variants:
            print("variants is not a Variants type, checking for multi-genome import")
            if isinstance(variants, list) and isinstance(variants[0], str):
                for file in variants:
                    if Path(file).is_file():
                        with open(file) as variantfile:
                            parsedchromosomes = {int(n.id):str(n.seq) for n in SeqIO.parse(reference_file, "fasta")}
                            referenceGenome.addStrain(VariantGenome(referenceGenome, species, parsedchromosomes))
        else:
            if variants:
                for strain in variants.strains:
                    referenceGenome.addStrain(VariantGenome(referenceGenome, strain, variants.of_strain(strain)))

        return referenceGenome
        #if variants are passed, we are going to build a bunch of new genome objects
        #the genomes should build & insert their variants themselves. 

class CSVHomologs(object):
    def __init__(self, homologyfile, name):
        self.csvfilename = homologyfile
        self.csvfilehandle = open(self.csvfilename)
        self.csv_df = pd.read_csv(self.csvfilehandle)

    def close(self):
        self.csvfilehandle.close()

class CSVaHomologs(CSVHomologs):
    '''accepts CSV homology files where row = homolog, column = species'''
    def __init__(self):
        pass

class CSVbHomologs(CSVHomologs):
    '''this file format has two groupings of species/strains, and the species are using the same reference genome'''
    def __new__(self, homologyfile, name):
        CSVHomologs.__init__(self, homologyfile, name)
        homologymap = HomologyMap(name, homologs = 'ALL')
        for column in self.csv_df:
            homologymap.addGroup(column)
            homologymap.addToGroup(column,[name + '_'+ strain for strain in self.csv_df[column][pd.notnull(self.csv_df[column])]])
        return homologymap

class SeqRecordPairing(object):
    def __init__(self, name, inseqs):
        self.group1 = inseqs[0]
        self.group2 = inseqs[1]
        self.name = name
        self.group1_identifiers = [ele.id for ele in self.group1]
        self.group2_identifiers = [ele.id for ele in self.group2]
        

    def edit_groups(self, edited_groups):
        self.group1 = [ele for ele in edited_groups if ele.id in self.group1_identifiers]
        self.group2= [ele for ele in edited_groups if ele.id in self.group2_identifiers]

    @property
    def unnested_groups(self):
        return [ele for ele in itertools.chain(self.group1, self.group2)]

    def write(self,outdir,txt_to_append=''):
        fileouts = [outdir.joinpath("GROUP1_%s_%s.txt" % (self.name, txt_to_append)), outdir.joinpath("GROUP2_%s_%s.txt" % (self.name, txt_to_append))]
        for out, file in zip([self.group1, self.group2], fileouts):
            with file.open('w') as writefile:
                SeqIO.write(out, writefile, "fasta")
        return fileouts
''' begin algorithm objects '''

class Aligner(object):
    '''Aligner: Requires MUSCLE, Biopython
    given an unnested list of seqrecord, align on a subprocess and return the alignment
    '''
    def __new__(self, inseqs):
        muscle_cline = MuscleCommandline(clwstrict=True)   
        self.child = subprocess.Popen(str(muscle_cline), stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True, shell=True, preexec_fn=os.setsid)
        SeqIO.write(inseqs, self.child.stdin, "fasta")
        self.child.stdin.close()
        alignment = AlignIO.read(self.child.stdout,"clustal")

        return alignment

class Trimmer(object):
    '''Trimmer:
    Given an alignment, trim all '-' and return sequences in their native structure
    '''
    def __new__(self, inseqs):
        inserts = []
        codons = []
        sequences = inseqs
        # this could probably be converted to np.arrays and use logical indexing to delete gaps for (?) more efficiency
        for index, sequence in enumerate(sequences): # go through each sequence in the alignment
            sequence = str(sequence.seq)
            codons.append([sequence[n:n+3] for n in range (0, len(sequence), 3)])
            inserts.append([c for c, x in enumerate(codons[index]) if '-' in x])
             # check the sequence for gaps ('-'), append a list of the indices of gaps for each alignment
        inserts=set([item for sublist in inserts for item in sublist]) # iterate through the nested list and condense into a set (removes duplicates)
        inserts = list(inserts) #flip it back to a list for sorting and indexing
        for i in range(len(codons)): # go back through each alignment
            record = codons[i] # set to current alignment
            for z in sorted(inserts,reverse=True): # index of gaps, in reverse so that the indices do not change
                del record[z]

            codons[i] = ''.join(record)
            sequences[i].seq=Seq.Seq(codons[i])
            
        return sequences        

class DirectoryTree(object):
    """builds an ouput directory tree"""
    __tree_dict = {"genomedir": "genomes", 
                    "pre_align": "pre_alignment", 
                    "alignments": "alignments",
                    "trimmed": "trimmed",
                    "csv": "csv"}
    def __init__(self, rootdir):
        
        self._rootdir = Path(rootdir)
        self.outdir = self.mksubdir(self._rootdir, "out")
        for variable, human_readable in self.__tree_dict.items():
            setattr(self, variable, self.mksubdir(self.outdir,human_readable))   


    def mksubdir(self, root, directory):
        directory = root.joinpath(str(directory))
        if not directory.is_dir():
            directory.mkdir()
        return directory

class Program(object):
    """docstring for program"""
    def __init__(self, 
                genomes, 
                homologymap,
                species, 
                annotations = None, 
                annotationsdb = None, 
                rootdir = '.', 
                variants = None, 
                ignore_small_sample = False):
        if not Path(genomes).is_file():
            raise AttributeError("No genomes supplied for input")
        if not Path(homologymap).is_file():
            raise AttributeError("No homology map supplied for polymorphism-divergence mapping")
        if not annotations and not annotationsdb:
            raise AttributeError("No annotationsdb supplied and no annotations map for import")
        if not variants and len(genomes) <= 2 and not ignore_small_sample:
            raise AttributeError("<= 2 genomes supplied and no variants to build, massprf analysis is incomplete")
        self.directories = DirectoryTree(str(rootdir))
        
        if variants:
            self.variants = VariantsAdaptorBuilder(str(variants))
        else:
            self.variants = None
        self.homologymap = HomologAdaptorBuilder(str(homologymap),species)
        self.genomes = GenomeAdaptorBuilder(str(genomes),str(species), variants = self.variants)
        self.homologymap.buildMapToGenomes(self.genomes)


        annotationsargs = {key:value for key,value in {'db':annotationsdb,'gff3':annotations}.items() if value}
        self.annotations = CDSAdaptorBuilder('gff3', **annotationsargs)


    def export_genomes(self):
        self.genomes.export(self.directories.genomedir)

def build_genes_parallel(inputs):
    cds = inputs[0]
    annotations = inputs[1]
    homologymap = inputs[2]
    genomes = inputs[3]
    directories = inputs[4]


    grouping = []

    for group in homologymap.group_list:
        outbuffer = []
        for genome in homologymap.group_genome_dict[group]:
            outbuffer.append(genome.genes_map[cds.name].getSequence(genome))
        outbuffer = [SeqRecord.SeqRecord(Seq.Seq(str(ele)), id=str(ele.strain)+'_'+str(ele.name)) for ele in outbuffer]
        grouping.append(outbuffer)

    grouping = SeqRecordPairing(cds.name, grouping)
    print(cds.name, ' generated \n')
    pre_alignedfiles = grouping.write(directories.pre_align, 'pre_alignment')
    alignment = Aligner(grouping.unnested_groups)
    print(cds.name, ' aligned \n')
    grouping.edit_groups(alignment)
    alignedfiles = grouping.write(directories.alignments, 'aligned')
    trimmed = Trimmer(grouping.unnested_groups)
    print(cds.name, ' trimmed \n')
    grouping.edit_groups(trimmed)
    trimmedfiles = grouping.write(directories.trimmed, 'trimmed')
    
    return [cds.name, 
            len(cds), 
            pre_alignedfiles[0], 
            pre_alignedfiles[1],
            alignedfiles[0],
            alignedfiles[1],
            trimmedfiles[0],
            trimmedfiles[1]]

def run(genes = 'ALL'):
    program = Program("../ricemassprf/referencegenome_12chro.fa",
                        '../ricemassprf/junruiricemap.csv',
                        'rice',
                        annotationsdb = '../ricemassprf/ricefeatureDB',
                        variants = '../ricemassprf/partfilerice')
    gene_annotations = map(lambda ele: program.annotations.createCodingAnnotation(ele, program.genomes), [gene for gene in program.annotations.getGenes(genes)])
    inputs = [[gene, 
                program.annotations,
                program.homologymap,
                program.genomes,
                program.directories] for gene in gene_annotations]

    #proc_pool = Pool(2)

    #output = proc_pool.map(build_genes_parallel, inputs)
    output = list(map(build_genes_parallel, inputs))
    print(output)

def test():
    #tree = DirectoryTree('.')
    #variants = VariantsAdaptorBui(lder("../ricemassprf/mass_rice")
    #genomes = GenomeAdaptorBuilder("../ricemassprf/referencegenome_12chro.fa", "rice", variants = variants)
    #genomes.export()
    #return genomes
    '''annotation_args = {'db':'../ricemassprf/ricefeatureDB'}
    annotations = CDSAdaptorBuilder('gff3',**annotation_args)
    target_genes = ['OS05T0522500-01']
    gene = annotations.getGenes(target_genes)
    output_annotation = []
    output_sequences = []
    for c in gene:
        output_annotation.append(annotations.createCodingAnnotation(c, genomes))
    for g in output_annotation:
        output_sequences.append(g.getSequence(genomes.name))

    return genomes, output_annotation, output_sequences'''


'''skeleton code for chromosome tester'''
'''
genome = main.test()
refchr = genome.get_chromosome(3)
var = genome[1]
varchr = var.get_chromosome(3)
refchr != varchr
'''
'''command line
if __name__ == '__main__':
    parser = ArgumentParser(description = "Pipeline for processing genome data for downstream MASSPRF analysis")
    parser.add_argument('-dir', dest = 'rootdir', )'''