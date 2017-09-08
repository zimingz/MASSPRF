'''This is an example of massprf_pipeline to illustrate how to use different classes'''
from pathlib import Path
from Bio import SeqIO, Seq, SeqRecord, AlignIO
import pandas as pd
import massprf_pipeline as pipeline
from multiprocessing import Pool
from argparse import ArgumentParser

class Program(object):
    """Initialize all relevant databases and objects for gene collation"""
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
        #initialize the output directories
        self.directories = pipeline.DirectoryTree(str(rootdir))
        #check for variants
        if variants:
            self.variants = pipeline.VariantsAdaptorBuilder(str(variants))
        else:
            self.variants = None
        #initialize groupings
        self.homologymap = pipeline.HomologAdaptorBuilder(str(homologymap),species)
        #initialize genomes w/ variants
        self.genomes = pipeline.GenomeAdaptorBuilder(str(genomes),str(species), variants = self.variants)
        #map the genomes to the groupings
        self.homologymap.buildMapToGenomes(self.genomes)

        #initialize annotations
        annotationsargs = {key:value for key,value in {'db':annotationsdb,'gff3':annotations}.items() if value}
        self.annotations = pipeline.CDSAdaptorBuilder('gff3', **annotationsargs)

    def export_genomes(self):
        self.genomes.export(self.directories.genomedir)

def build_genes_parallel(inputs):
    '''an example gene mapping function that builds groups of genes based on input and aligns/trims them
    returns a list of relevant files'''
    #input parsing
    cds = inputs[0]
    annotations = inputs[1]
    homologymap = inputs[2]
    genomes = inputs[3]
    directories = inputs[4]

    #grouping is a list of lists, each sublist is a group of seqrecords
    grouping = []

    for group in homologymap.group_list:
        #for each group, build a buffer of SeqRecord(genes)
        outbuffer = []
        for genome in homologymap.group_genome_dict[group]:
            outbuffer.append(genome.genes_map[cds.name].getSequence(genome))
        outbuffer = [SeqRecord.SeqRecord(Seq.Seq(str(ele)), id=str(ele.strain)+'_'+str(ele.name)) for ele in outbuffer]
        #id is used because the alignment function preserves ID keys
        #write to grouping
        grouping.append(outbuffer)

    grouping = pipeline.SeqRecordPairing(cds.name, grouping) #init a pairing so we can keep track of the genes and write them
    print(cds.name, ' generated \n')
    pre_alignedfiles = grouping.write(directories.pre_align, 'pre_alignment') #write the pre aligned genes
    alignment = piepline.Aligner(grouping.unnested_groups) #build an alignment
    print(cds.name, ' aligned \n')
    grouping.edit_groups(alignment) #sync the alignment back into its groupings
    alignedfiles = grouping.write(directories.alignments, 'aligned') #write the grouped alignments
    trimmed = pipeline.Trimmer(grouping.unnested_groups) #initialize a trimming
    print(cds.name, ' trimmed \n')
    grouping.edit_groups(trimmed) #sync trimming back into grouping
    trimmedfiles = grouping.write(directories.trimmed, 'trimmed') #write the trimming
    
    return [cds.name, 
            len(cds), 
            str(pre_alignedfiles[0]), 
            str(pre_alignedfiles[1]),
            str(alignedfiles[0]),
            str(alignedfiles[1]),
            str(trimmedfiles[0]),
            str(trimmedfiles[1])]



if __name__ == "__main__":
    parser = ArgumentParser(description='example massprf-pipeline usage for a single genome with a one-to-one homolog mapping')
    parser.add_argument('-ref', dest='reference_genome', required=True, type=str, nargs=1, help="reference genome file")
    parser.add_argument('-hom', dest='homology_map', required=True, type=str, nargs=1, help="CSV homology map in format of example.csv")
    parser.add_argument('-spc', dest='species', required=True, type=str, nargs=1, help= 'placeholder species name')
    parser.add_argument('-anno', dest='annotations', required=False, type=str, nargs=1, help='annotations file if db not used')
    parser.add_argument('-annodb', dest='annotations_db', required=False, type=str, nargs=1, help='annotations database if initialized')
    parser.add_argument('-cores', dest='cores', required = True, type=str, nargs=1, help='number of cores desired to use')
    parser.add_argument('-vcf', dest='vcf', required = True, type=str, nargs=1, help = 'vcf file')

    args = parser.parse_args()

    if not (args.annotations or args.annotations_db):
        raise Exception("No annotations or database provided")
    elif args.annotations and not args.annotations_db:
        program = Program(args.reference_genome[0], 
                            args.homology_map[0],
                            args.species[0],
                            annotations = args.annotations[0],
                            variants = args.vcf[0])

    elif args.annotations_db and not args.annotations:
        program = Program(args.reference_genome[0], 
                            args.homology_map[0],
                            args.species[0],
                            annotationsdb = args.annotations_db[0],
                            variants = args.vcf[0])
    gene_annotations = map(lambda ele: program.annotations.createCodingAnnotation(ele, program.genomes), [gene for gene in program.annotations.getGenes()])
    print("program initialized")
    inputs = [[gene, 
                program.annotations,
                program.homologymap,
                program.genomes,
                program.directories] for gene in gene_annotations]
    if cores <= 2:
        output = list(map(build_genes_parallel, inputs))
    else:
        proc_pool = Pool(cores/2) #divide by 2 because alignment uses a separate subprocess
        output = proc_pool.map(build_genes_parallel, inputs)

    output = sorted(output, key= lambda ele: ele[1]) #sort by gene length

    output_df = pd.DataFrame(data = output, columns = ['gene','length','group1_pre_aligned','group2_pre_aligned','group1_aligned','group2_aligned','group1_trimmed','group2_trimmed'])

    output_df.to_csv(program.directories.csv.joinpath('output.csv').open('w'))
