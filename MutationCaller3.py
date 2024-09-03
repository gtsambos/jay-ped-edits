import math
import multiprocessing as mp
import functools
from collections import defaultdict

import vcf
import sys
import numpy
from numpy import exp
from numpy import log1p
import concurrent.futures

#TODO: Handle Phase and IBD uncertainty better
#TODO: Finish testing internal statistical behaviors
#TODO: Repair parallelization (This requires that I make a forked version of the pyvcf package that uses the correct settings for multiple iterators, which the current cannot do)
#TODO: Finish integrating two separate IBD
ibd_error_rate_flat = 0.01

def sum_log_prob(a, b):
    if a > b:
        return a + log1p(exp(b-a))
    else:
        return b + log1p(exp(a-b))


#Input:
#ibdKnown;  if IBDdata is available. Currently used as always true without a secondary input file containing doubt pairings.
#phased; Phasing status
#baseLikelihood; Used score for IBD confidence. Combination of experimentally determined and from two possible IBD segment files of origin

#Output:
def ibd_probability(ibdKnown, phased, distance, baseLikelihood):
    if(ibdKnown):
        return baseLikelihood
    elif(phased):
        return 0.5*baseLikelihood
    else:
        return math.pow(0.5,distance)*baseLikelihood

#Input:
#focal_allele; value of allele currently assumed to be the true De Novo mutation
#node_of_read_origin; node structure for the current read
#ibd_likelihood; likelihood score for the focal sample haplotype of origin for the allele and haplotype of the read to be IBD
#descendant; boolean for if the read of origin comes from the focal sample haplotype or a descendant of the focal sample.
#Output;
def read_likelihood(focal_allele, node_of_read_origin, ibd_likelihood, descendant):
    if(descendant and (focal_allele==node_of_read_origin.allele)):
        return (1.0-node_of_read_origin.hidden_var_error)*ibd_likelihood+node_of_read_origin.hidden_var_error*ibd_likelihood
    elif(descendant and (focal_allele!=node_of_read_origin.allele)):
        return (1.0-node_of_read_origin.hidden_var_error)*(1-ibd_likelihood)+(node_of_read_origin.hidden_var_error)*ibd_likelihood
    elif((not descendant) and (focal_allele==node_of_read_origin.allele)):
        return (1.0-node_of_read_origin.false_var_error)*(1-ibd_likelihood)+(node_of_read_origin.false_var_error)*ibd_likelihood
    else:#if(not descendent and (focal_allele!=node_of_read_origin.allele)
        return (1-node_of_read_origin.false_var_error)

#TODO: possibly set up cache memoization to take advantage of repeated calls. Needs to use non-recursive so that function calls can be separated from relative distance between nodes, making their output identical from repeated calls
#Input:
#Output
def gatherLogLikelihood(viewed_ids, current_id, original_allele, original_id, original_haplotype, site, distance, descendent, outer_pedigree, inner_pedigree, seg_structure, ibd_error_rate):
    #EndCase, prevents redundent analysis
    if (current_id == -1 or current_id in viewed_ids):
        return 0.0
    viewed_ids.add(current_id)
    log_likelihood = 0.0
    #Iterates over possible haplotypes for current sample, accounts for read evidence in each depending on IBD beliefs
    for read_haplotype in inner_pedigree.map_to_inner_nodes[current_id]:
        inner_node = inner_pedigree.map_to_inner_nodes[current_id][read_haplotype]
        baseLikelihood = ibd_error_rate
        if(seg_structure.is_ibd(original_id, original_haplotype, current_id, read_haplotype, site)):
            baseLikelihood=1.0-baseLikelihood
        ibd_likelihood = ibd_probability(True, inner_node.phased, distance, baseLikelihood)
        log_likelihood += inner_node.depth * math.log(read_likelihood(original_allele, inner_node, ibd_likelihood, descendent))
    for parent in outer_pedigree[current_id].parents:#set descendent status to False when going into a parent
        log_likelihood += gatherLogLikelihood(viewed_ids, parent, original_allele, original_id, original_haplotype, site, distance + 1, False, outer_pedigree, inner_pedigree, seg_structure, ibd_error_rate)
    for child in outer_pedigree[current_id].children:
        log_likelihood += gatherLogLikelihood(viewed_ids, child, original_allele, original_id, original_haplotype, site, distance + 1, descendent, outer_pedigree, inner_pedigree, seg_structure, ibd_error_rate)
    return log_likelihood#TODO Need to also return IBD presence in decendents to shore up IBD confidence concerns


#In its final version, this is meant to be used by parallel processes. Currently, it is not.
#Input:
#contig_iterator; an iterator aimed at a range of sites within a tabix-indexed vcf. When parallelized, is one of many aimed at the same file
#outer_tree; a tree containing sample names and parent/child relationships. When parallelized, all threads are looking at the same shared tree object.
#seg_first_order_filename; filename for the segments ordered by their startpoint. Is not an actual opened file yet.
#seg_last_order_filename; filename for the segments ordered by their endpoint. Is not an actual opened file yet.
#q; the listener for parallel printing
#Output:
#None, but prints likelihood values per variant per sample, ordered in the same order as the input vcf, to a file. When parallelized, as a result of no control over which sites are finished first, the output is not ordered by site, and must be sorted in post-processing
def gather_likelihoods_for_variants(contig_iterator, outer_tree, seg_first_order_filename, seg_last_order_filename, q):
    with open(seg_first_order_filename) as first_ordered_file:
        with open(seg_last_order_filename) as last_ordered_file:
            structure = SegSiteStructure(first_ordered_file, last_ordered_file)#Starts a segment tracker for an efficient storage of IBD near the current site.
            site_record = contig_iterator.next()
            while(site_record):
                inner_tree = InnerTree(site_record)#builds a site-specific inner tree structure from the genotype calls at this site
                site_record = contig_iterator.next()
                likelihoods_to_print = str(site_record.site)
                for sample in site_record.samples:

                    viewed_ids = set()
                    haplotype = 0
                    original_node = inner_tree.map_to_inner_nodes[sample][haplotype]
                    hap_likelihood_0 = gatherLogLikelihood(viewed_ids, sample, original_node.allele, sample, haplotype, site_record.site, 0, True, outer_tree, inner_tree, structure, ibd_error_rate_flat)
                    viewed_ids = set()
                    haplotype = 1
                    original_node = inner_tree.map_to_inner_nodes[sample][haplotype]
                    hap_likelihood_1 = gatherLogLikelihood(viewed_ids, sample, original_node.allele[haplotype], sample,
                                                           haplotype, site_record.site, 0, True, outer_tree, inner_tree,
                                                           structure, ibd_error_rate_flat)

                    likelihoods_to_print += "\t{:.2f}".format(sum_log_prob(hap_likelihood_0,hap_likelihood_1))
            print(likelihoods_to_print)
            #q.put(likelihoods_to_print)

#Contains identifying information to be able to add or remove IBD segments from the current site
class Seg:
    #TODO: Introduce segment endpoint error that buffers edges of overlap checks
    #Populate segment information from input file line
    def __init__(self, seg_file_line):
        segTokens = seg_file_line.split()
        self.name1 = segTokens[0]
        self.haplotype1 = int(segTokens[1])
        self.name2 = segTokens[2]
        self.haplotype2 = int(segTokens[3])
        #Base coordinates of segment endpoints
        self.start = int(segTokens[4])
        self.end = int(segTokens[5])

    #Booleans to check if this segment overlaps a given site
    def overlap(self, site):
        return self.start < site and self.end > site

    def check_starts_before(self, site):
        return self.start <= site

    def check_starts_after(self, site):
        return self.start > site

    def check_ends_before(self, site):
        return self.end <= site

    def check_ends_after(self, site):
        return self.end > site

    #Used for accessing if two haplotypes/samples are IBD
    def get_pairing_id(self):
        return (self.name1,self.haplotype1,self.name2,self.haplotype2)

    def get_reverse_pairing_id(self):#NOTE: there may be a better way to handle segment orientation consistency in hashing. There are time/memory pros/cons.
        return (self.name2,self.haplotype2,self.name1,self.haplotype1)


#Contains the connected datastructures involved in IBD segment traversal and tracking. Adds and removes segments that overlap with the site being investigated
class SegSiteStructure:
    def __init__(self, first_ordered_file, last_ordered_file):
        self.first_ordered_file = first_ordered_file#used to efficiently remove segments from analyses
        self.last_ordered_file = last_ordered_file#Used to efficiently add segments to analyses
        self.seg_map = {} #Used to access IBD status of pairs of samples
        self.doubt_map = set() # TODO: set up input to account for uncertainty in IBD

    def add_segment(self, seg):
            segID = seg.get_pairing_id()
            if(not segID in self.seg_map):
                #segList.add(seg)
                self.seg_map[segID]=seg
                self.seg_map[seg.get_reverse_pairing_id()]=seg
    def remove_segment(self, seg):
            segID = seg.get_pairing_id()
            if(segID in self.seg_map):
                #segList.add(seg)
                del self.seg_map[segID]
                del self.seg_map[seg.get_reverse_pairing_id()]

    #Allows querying a line from a file and then putting it back in its previous interal position.
    def peek_line_last(self):
        pos = self.last_ordered_file.tell()
        line = self.last_ordered_file.readline()
        self.last_ordered_file.seek(pos)
        return line
    def peek_line_first(self):
        pos = self.first_ordered_file.tell()
        line = self.first_ordered_file.readline()
        self.first_ordered_file.seek(pos)
        return line

    #Checks what segments overlap with current cite, adds and removes tracked segments as needed
    def updateSegs(self, site):
        last_ordered_line = self.peek_line_last()
        # Moves to first potentially overlapping segment, and attempts removal of any segments that no longer overlap after site update.
        while (last_ordered_line):
            last_ordered_seg = Seg(last_ordered_line)
            if(last_ordered_seg.check_ends_before(site)):
                self.remove_segment(last_ordered_seg)
                self.last_ordered_file.readline()
            else:
                break
            last_ordered_line = self.peek_line_last()
        first_ordered_line = self.peek_line_first()
        #Moves to first potentially overlapping segment.
        while(first_ordered_line):
            first_ordered_seg = Seg(first_ordered_line)
            if(first_ordered_seg.check_starts_before(site)):
                if(first_ordered_seg.overlap(site)):
                    self.add_segment(first_ordered_seg)
                first_ordered_line = self.first_ordered_file.readline()
            else:
                break
            first_ordered_line = first_ordered_line = self.peek_line_first()


    #Checks if the current site is registered as IBD for the pair of samples and haplotypes, or if the samples/haplotypes refer to the same haplotype.
    def is_ibd(self, sample1, haplotype1, sample2, haplotype2, site):
        return ((sample1,haplotype1) == (sample2,haplotype2)) or (((sample1, haplotype1, sample2, haplotype2) in self.seg_map) and self.seg_map[(sample1, haplotype1, sample2, haplotype2)].overlap(site))







#Structure for representing the overall tree structure. Not unique to any site
class OuterNode:
    def __init__(self, name, mother, father):
        self.name = name
        self.parents = set()
        self.children = set()


    def add_child(self, child_name):
        self.children.add(child_name)

    def add_parent(self, parent_name):
        self.parents.add(parent_name)

#Site, sample, and haplotype specific node class with depth and allele data. Remade for each site analyzed.
#Represents a haplotype for a specific sample, and its allele controls what allele reads associated with it show for matching purposes
#For homozygous sites, there's a hypothetical node (hap -1) made for the other allele with zero coverage and no given IBD, and the two real haplotypes are given half read evidence.
class InnerNode:
    def __init__(self, allele, read_depth, false_var_error, hidden_var_error):

        self.depth = read_depth
        self.allele = allele
        self.false_var_error = false_var_error
        self.hidden_var_error = hidden_var_error
        self.phase = True




#Knows haplotype node alternatives, maps homozygous alleles to the same InnerNode with half read depth
class InnerTree:
    #Input:
    #data_entry; site record from the vcf iterator.
    def __init__(self, data_entry):
        self.map_to_inner_nodes = defaultdict(lambda: {}) #map of sample name to maps of haplotypes to inner nodes
        self.ordered_samples = []  # array for consistent ordered sample access based on the input file order
        for sample in data_entry.samples:
            self.call_to_nodes(data_entry.samples[sample])

    #CURRENTLY ASSUMES BIALLELIC
    #loads a node for the inner tree using a call datapoint from the input record
    def call_to_nodes(self, input_call):
        name = input_call.sample
        print(input_call.data)
        depth_values = input_call.data.AD
        pl_values = []
        self.ordered_samples.append(name)
        if getattr(input_call.data, 'PL', None) is not None:
            pl_values_log=input_call.data.PL
            pl_values = [math.pow(10, (pl_values_log[0])/-10.0), math.pow(10, (pl_values_log[1])/-10.0), math.pow(10, (pl_values_log[2])/-10.0)]
        else:
            gq_value = input_call.data.GQ
            pl_values=[gq_value, gq_value, gq_value]
            pl_values[input_call.gt_type]=0#The GQ value is, in practice, the difference between the two most likely PL values, so the correct genotype should have a PL of zero, and the others should have GQ or higher.
        if(input_call.gt_type==0):
            node = InnerNode(0, depth_values[0]/2, pl_values[1], pl_values[2])
            self.map_to_inner_nodes[name][0] = node
            self.map_to_inner_nodes[name][1] = node
        elif(input_call.gt_type==2):
            node = InnerNode(1, depth_values[1]/2, pl_values[1], pl_values[0])
            self.map_to_inner_nodes[name][0] = node
            self.map_to_inner_nodes[name][1] = node
        elif(input_call.gt_type==1):
            node1 = InnerNode(0, depth_values[0], pl_values[0], pl_values[2])
            node2 = InnerNode(1, depth_values[1], pl_values[2], pl_values[0])
            self.map_to_inner_nodes[name][0] = node1
            self.map_to_inner_nodes[name][1] = node2












#Tree structure consistent across all individual site runs. Built from pedigreed input
class OuterTree:
    def __init__(self, treeFileName):
        self.samples = {}#Map for fast sample access
        self.load_tree(treeFileName)

    #takes name of a tsv file to use as a pedigree. Still needs a method for taking in IBD uncertainty
    #Expects: sampleID motherID fatherID sex family
    #Does NOT check consistency for sex in input.
    def load_tree(self, treeFileName):
        with open(treeFileName) as treeFile:
            self.sampleCount = 0
            self.trueSampleCount = 0
            for line in treeFile:
                tokens = line.split()
                self.samples[tokens[0]] = OuterNode(tokens[0],tokens[1],tokens[2])


#input:
#vcf_reader; iterator for the entire input
#pool_count; number of processes available simultaneously

#output:
#list of iterators, n per contig, where n is number of pools, so one contig can be processed in parallel at a time.
def vcf_to_range_list(vcf_reader, pool_count):
    ranges = []
    for contig in vcf_reader.contigs:
        try:
            length = vcf_reader.contigs[contig][1]
            step = length/pool_count
            start = 0
            while(start<length):
                ranges.append((contig, start, start+step))
                start+=step
        except:#Avoids an issue with asking for fetching an empty contig region, for which tabix throws an error rather than giving an empty object.
            continue
    return ranges


#input:
#vcf_reader; iterator for the entire input
#pool_count; number of processes available simultaneously

#output:
#list of iterators, assuming one given contig with a string name, where n is number of pools.
def single_contig_vcf_to_range_list(vcf_reader, pool_count, contig):
    ranges = []
    length = vcf_reader.contigs[contig][1]
    step = length / pool_count
    start = 0
    while (start < length):
        ranges.append((contig, start, start + step))
        start += step
    return ranges

#Parallel process printing function used by one of the pool processes
def listener(q, fn):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'w') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                f.write('killed')
                break
            f.write(str(m) + '\n')
            f.flush()

def get_call(iterator_range):
    vcf_reader = vcf.Reader(filename="TestInput/183.interval_list_shailee_run_1.vcf.gz")
    iterator = vcf_reader.fetch(iterator_range[0],iterator_range[1],iterator_range[2])
    record = next(iterator, None)
    if(record == None):
        return -1
    else:
        return record.start

def main():
    #Sets a parallelization strategy that avoids duplicating things not passed to individual processes, and is robust to multiple POSIX systems.
    #DOES use one thread to manage processes, so one core is going to be devoted towards handling the other spawned processes


    mutation_candidate_vcf_name = sys.argv[1]
    tree_file_name = sys.argv[2]
    seg_first_order_filename = sys.argv[3]
    seg_last_order_filename = sys.argv[4]

    output_filename = sys.argv[6]

    pool_jobs = int(sys.argv[5])# mp.cpu_count() can theoretically get available cores

    outer_tree = OuterTree(tree_file_name)

    vcf_reader = vcf.Reader(filename=mutation_candidate_vcf_name)
    #gather_likelihoods_for_variants(vcf_reader, outer_tree, seg_first_order_filename, seg_last_order_filename, q)


    #TODO: The parallel process I was attempting here does not work, will replace with a placeholder that performs the same task shortly, then will repair the parallel version.
    #When working properly, splits up contigs into subregions, assigns iterators to each, then within each region, iterates site by site, using the segment updates to store overlapping subsets of IBD segments.
    #Then iterates over each sample, and performs the recursive likelihood function gatherLogLikelihood. For the current statistical model, it doesn't use reads from outside the subtree

    # mp.set_start_method('spawn')
    # manager = mp.Manager()
    # pool = mp.Pool(processes=pool_jobs)
    # q = manager.Queue()  # Used for multiprocess-safe printing to an output file.
    # watcher = pool.apply_async(listener, (q, output_filename))
    # vcf_reader = vcf.Reader(filename=mutation_candidate_vcf_name)

    # iterators_for_pools = []

    #q.put('kill')
    #pool.close()
    #pool.join()

    vcfReader = vcf.Reader(filename="TestInput/183.interval_list_shailee_run_1.vcf.gz")
    # iterators = vcf_to_iterator_list(vcfReader, 3)
    iterators = single_contig_vcf_to_range_list(vcfReader, pool_jobs, "Chr18")
    with concurrent.futures.ProcessPoolExecutor(max_workers=pool_jobs) as executor:
        for callsite in zip(executor.map(get_call, iterators)):
            print(callsite)
if __name__ == "__main__":
    main()
