import unittest
import itertools
import collections
import multiprocessing as mp
import functools
from collections import defaultdict

import vcf
import sys
import numpy
from MutationCaller3 import *

class VCF_Input_Tests(unittest.TestCase):
    def notest_create_inner_node_from_call(self):
        vcfReader = vcf.Reader(filename="TestInput/183.interval_list_shailee_run_1.vcf.gz")
        variant1 = next(vcfReader)
        call = variant1.samples[0]
        print(variant1)
        inner_tree = InnerTree()
        inner_tree.call_to_nodes(call)#Should be a homozygous ref call
        self.assertEqual(len(inner_tree.map_to_inner_nodes),1)
        self.assertEqual(len(inner_tree.map_to_inner_nodes["APA1"]), 2)
        self.assertEqual(inner_tree.map_to_inner_nodes["APA1"][0],inner_tree.map_to_inner_nodes["APA1"][1])
        self.assertEqual(inner_tree.map_to_inner_nodes["APA1"][0].allele, 0)

    def nost_vcf_input_size_check(self):
        vcfReader = vcf.Reader(filename="TestInput/183.interval_list_shailee_run_1.vcf.gz")
        #elem_count = itertools.count()
        #collections.deque(zip(vcfReader, elem_count))
        #print('Iterator has {} elements'.format(next(elem_count)))
        self.assertEqual(True, True)  # add assertion here


class Seg_Input_Tests(unittest.TestCase):
    def test_segment_updates_1(self):
        with open("TestInput/segsTest1.seg.last_order") as last_ordered_file:
            with open("TestInput/segsTest1.seg.first_order") as first_ordered_file:
                structure = SegSiteStructure(first_ordered_file, last_ordered_file)
                structure.updateSegs(5)
                self.assertTrue(structure.is_ibd("sample1", 0, "sample3", 0,5))
                self.assertFalse(structure.is_ibd("sample1", 1, "sample3", 0,5))

    def test_segment_update_and_remove_1(self):
        with open("TestInput/segsTest1.seg.last_order") as last_ordered_file:
            with open("TestInput/segsTest1.seg.first_order") as first_ordered_file:
                structure = SegSiteStructure(first_ordered_file, last_ordered_file)
                structure.updateSegs(5)
                self.assertTrue(structure.is_ibd("sample1", 0, "sample2", 0,5))
                self.assertTrue(structure.is_ibd("sample1", 0, "sample3", 0,5))
                self.assertFalse(structure.is_ibd("sample2", 0, "sample3", 0,5))
                structure.updateSegs(999999)
                self.assertTrue(structure.is_ibd("sample1", 0, "sample2", 0, 999999))
                self.assertFalse(structure.is_ibd("sample1", 0, "sample3", 0, 999999))
                self.assertTrue(structure.is_ibd("sample2", 0, "sample3", 0, 999999))

class Pedigree_Input_Tests(unittest.TestCase):
    def test_tree_input(self):
        outer_tree = OuterTree("TestInput/testPedigree1.ped")
        self.assertEquals(len(outer_tree.samples),3)
        self.assertTrue(True)


class Parallelism_Tests(unittest.TestCase):
    def get_call(self, iterator):
        call = iterator.next().samples[0]
        return call.sample.site.pos

    def test_pooled_vcf_iterators(self):
        mp.set_start_method('spawn')
        vcfReader = vcf.Reader(filename="TestInput/183.interval_list_shailee_run_1.vcf.gz")
        #iterators = vcf_to_iterator_list(vcfReader, 3)
        iterators = single_contig_vcf_to_iterator(vcfReader, 3, "Chr18")
        with mp.Pool(processes=3) as pool:
            multiple_results = [pool.apply_async(self.get_call, (iterators[i])) for i in range(3)]
            self.assertTrue(multiple_results[0].get()<multiple_results[1].get())
            self.assertTrue(multiple_results[1].get()<multiple_results[2].get())



if __name__ == '__main__':
    unittest.main()
