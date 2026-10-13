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

## Other members

### Aamna
I actually found this project a lot easier to understand conceptually than some of the previous ones, which is surprising because genome assembly in general is a pretty complicated process. Breaking sequencing reads into overlapping k-mers and using those overlaps to build a De Bruijn graph made sense to me. We discussed De Bruijn graphs earlier in the program too which helped because I had a biological application I was already somewhat familiar with so I understood what we were trying to do before getting into the code.

I particularly liked working through the add_edge()/remove_edge() and build_graph_from_reads() functions. Once I understood how each k-mer was split into a prefix and suffix of length k-1, the graph construction clicked and was simple. The toy example helped a lot because I could actually follow how the reads were connected.

The Eulerian walk was probably the most challenging part for me. I understood the idea of traversing every edge exactly once, but recursion took longer to figure out. I had to go through the code a few times to understand why we append nodes after visiting their outgoing edges and why the tour needs to be reversed. I also struggled a little with assemble_contigs(), especially understanding how to choose a starting node and what happens when there are multiple possible paths.

One thing I liked was comparing different approaches to the same problem. I followed the recursive Eulerian walk from Langmead's lecture, while Graziano used an iterative approach with a stack. I also hadn't really thought about what happens when a graph branches and there isn't one clear path to follow. Initially, I was mainly thinking about connecting overlapping sequences, but looking at non-branching paths helped me understand why we can't always assume which direction a sequence should continue.

Overall, I enjoyed this project. I still want to get more comfortable with recursion and understanding how different functions fit together, but I felt more confident with the underlying algorithm this time. I think this project helped me connect the coding side of bioinformatics with the actual biological problem we were trying to solve.

# Generative AI Appendix
As per the syllabus
