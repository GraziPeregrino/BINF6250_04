# Introduction
A De Bruijn Graph visualizes the *de novo* genome assembly of short-read sequencing to construct longer consensus sequences (contigs). To construct the contigs, a Eulerian walk (path) is used to connect the nodes (prefixes and suffixs of the k-1 mer) with edges (the overlaps). 

This project is broken down into three challenges: using toy/ sample data, using the yeast genome and attempting the mouse genome. In this code, we implement object-oriented programming by building ```class DeBruijn Graph``` to manipulate sequencing reads to conduct genome assembly by using Eulerian paths. 

# Pseudocode
Put pseudocode in this box:

```
class DeBruijnGraph

def __init__(self, reads, k):

  define variables reads and k

    reads = DNA sequences as a string

    k = k-mer as a number

def add_edge(self, left, right):

  pair prefixes with suffixes to create an edge

  check for duplicates with an "if" statement

    if not in the current dictionary, append the next k-1 mer

def remove_edge(self, left, right):

  opposite but works in tangent with add_edge() function

  removes edges directed edges from graph with a guard check: checks for the edge and the final stop

def build_graph_from_reads(self, reads, k):

  calls the add_edge function to run many times

  loop the sequencing reads

    loop the reads in the length of the k-mer

    slice the k-mer, split the k-mer and call add_edge() function

def eulerian_walk(self, node, graph, seed=None):

  select a random edge with seed, remove it and recurse

  stack = the Eulerian path being walked from the list of nodes

  loop stack

    start at a node

    continue matching prefixes and suffixes

    keep track of the nodes and edges

    return the string of sequences (tour) 

  tour = the completed sequence from the Eulerian path

def assemble_contigs(self, seed=None):

  construct the contigs from the Eulerian walk

  ~ need to think more about this ~

def tour_to_sequence(self, tour):

  ~ need to think more ~ 

def get_assembly_stats(self, contigs):

  create a list of the lengths of the contigs sorted from smallest to largest 

  number of contigs = length of each of the newly assembled contigs

  total length = the sum of the lengths of the contig

  largest contig = the biggest in length of the list (should be the last number so -1)

  shortest contig = smallest in length (should be the first number in the list so 0) 

  mean length = the average of all the lengths of the contig

  n50 = length of the contig at which the running total, adding longest

def write_fasta(self, contigs, filename):

  open and write a new FASTA file

  create a header for each with the contig number and the length

  export the file 

```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
