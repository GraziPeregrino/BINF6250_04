# Introduction
After short-read sequencing, the output is millions of short fragments (reads) of about 100-300 bp. To reconstruct the genome, we use de novo assembly because it does not require a reference genome. Comparing and reading each read against every other read to find overlaps is computationally long and difficult at this scale. The De Bruijn graphs can mitigate some of this overload by splitting each read into overlapping substrings of length k (k-mers).

In a De Bruijn graph, each k-mer is a directed edge. Its first k−1 bases (the prefix) and last k−1 bases (the suffix) are the nodes it connects. Two k-mers that overlap by k−1 bases share a node, so overlaps between reads appear as connected edges. Assembly then becomes a path problem. A Eulerian path, which traverses every edge exactly once, uses every k-mer exactly once, and the sequence of nodes along it spells out a contig, a contiguous assembled sequence.

This project implements that algorithm in Python using object-oriented programming. The ```class DeBruijnGraph``` stores the graph as an adjacent list, builds it from reads, performs the Eulerian walk, converts paths to sequences, and reports assembly statistics. We assume error-free reads, which removes the need for error correction. We take into account two features of real data: reads come from both DNA strands, so reverse-complement duplicate contigs must be removed, and repeats and coverage gaps mean no single Eulerian path exists, so contigs break at those points.

The project has three stages of increasing scale: (1) a small toy example to verify correctness, (2) the S. cerevisiae (yeast) genome (~12 Mb, 2.4 million simulated 150 bp reads at 30× coverage), and (3) the mouse genome, which tests how the approach scales in time and memory. Assemblies are evaluated by number of contigs, total length, longest contig, and N50.

```text
CLASS DeBruijnGraph
    FUNCTION __init__(reads, k):
        CREATE empty graph
        STORE k-mer size
        CALL build_graph_from_reads(reads, k)

    FUNCTION add_edge(left, right):
        IF edge left → right does not already exist:
            ADD right as an outgoing neighbor of left

    FUNCTION remove_edge(left, right):
        IF edge left → right exists:
            REMOVE the edge

    FUNCTION build_graph_from_reads(reads, k):
        FOR each read:
            FOR each k-mer in the read:
                left  = prefix of length k - 1
                right = suffix of length k - 1
                CALL add_edge(left, right)

FUNCTION eulerian_walk(node, graph, seed=None):
    IF seed is provided:
        SET the random seed
    CREATE tour containing node
    SET current ← node
    WHILE current exists in graph AND current has outgoing edges:
        RANDOMLY SELECT one outgoing edge from current
        SET next_node ← destination of selected edge
        REMOVE the selected edge from graph[current]
        IF graph[current] has no remaining outgoing edges:
            REMOVE current from graph
        APPEND next_node to tour
        SET current ← next_node
    RETURN tour

FUNCTION assemble_contigs(seed=None):
    IF seed is provided:
        SET the random seed
    CREATE work as a copy of the graph
        COPY every node and its outgoing edges
        INCLUDE only nodes that have outgoing edges
    CREATE an empty dictionary called indegree
    FOR each node in work:
        FOR each destination node connected from node:
            INCREASE indegree[destination] by 1
    CREATE an empty list called start_nodes
    CREATE an empty list called other_nodes
    FOR each node in work:
        SET outgoing_degree ← number of outgoing edges from node
        SET incoming_degree ← indegree[node]
        IF outgoing_degree > incoming_degree:
            APPEND node to start_nodes
        ELSE:
            APPEND node to other_nodes
    CREATE an empty list called contigs
    CREATE ordered_nodes by combining:
        start_nodes followed by other_nodes
    FOR each start node in ordered_nodes:
        WHILE start still exists in work:
            SET tour ← eulerian_walk(start, work)
            IF tour contains more than one node:
                SET sequence ← tour_to_sequence(tour)
                APPEND sequence to contigs
    RETURN contigs

    FUNCTION tour_to_sequence(tour):
        IF tour is empty:
            RETURN empty sequence
        START sequence with the first (k - 1)-mer
        FOR each remaining node:
            APPEND the last nucleotide of the node
        RETURN assembled sequence

    FUNCTION get_assembly_stats(contigs):
        CALCULATE length of each contig
        SORT lengths from longest to shortest
        CALCULATE:
            number of contigs
            total assembly length
            longest contig
            shortest contig
            mean contig length
        CALCULATE N50:
            ADD contig lengths from longest to shortest
            STOP when cumulative length reaches half
            of the total assembly length
        RETURN assembly statistics

    FUNCTION write_fasta(contigs, filename):
        OPEN output FASTA file
        FOR each contig:
            WRITE contig header
            WRITE contig sequence
        CLOSE file
```
## Output 

Please review the output files named genome_assembly_stats.txt.

## Successes

There were major successes while working on this project. 

*Scalability of the De Bruijn graph algorithm*: running from a small toy dataset to more than 2.4. million yeast sequencing reads. Using the small dataset is the best approach to validate the code and the graph construction. 

*Implementation of the program using Object Oriented Programming*: Use class to separate responsibilities into methods for the graph construction was an improvement to the code. 

*Performance optimization*: While running the code I called the Add_edge into the build_graph_from_reads function and notice that the code would run correct with the toy dataset, but once we scale it, we notice that the time that the code was running increased. So I removed the add_edge calling and added a conditional instatement inline to the function and that improved the code. 

*Output Result*: We were successfully able to generate an output for the toy dataset as well as the yeast genome. 

## Struggles

We would like to bring few challenges that we faced as we developed the code.  

*Memory usage- an initial implementation stored forward and reverse graph*: Maintaining both adjacency structures required a large amount of memory. The final implementation stores the forward adjacency graph and only the incoming degree of each node, which contains the information required by the config assembly algorithm without storing every incoming adjacency list.

*Balancing the algorithm described in the assignment with the requirements of a large real dataset*: The Eulerian walk works well for demonstrating the required graph-traversal algorithm and for validating small examples. For the complete yeast graph, the final implementation constructs maximal non-branching paths instead. This required understanding why a node with one incoming edge and one outgoing edge represents an unambiguous continuation of a `contig` and why branching nodes should terminate it.

*Timing and Optimization*: An interesting challenge that came up throughout this project was the timing of the code within different IDEs. Using Pycharm on the Yeast genome took about an hour to complete; however, when using VSCode, the code took 5-8 minutes to complete. Optimization was a difficult task to complete as well. We tried our best to optimize the time, but with using different IDEs it became difficult to consider if the code is running at a faster time. We would love to get insight on this struggle and see if this happened to other groups as well. 

## Personal Reflections
### Group Leader

### Kailey Solorzano

The overall project changed my logical thinking. I received a brief introduction to DeBruijn Graphs in BINF6310, but nothing as extensive as this course. The concept was difficult to learn at first (I had way too many tabs open), but after learning how to reconstruct it by hand provided the basis for my pseudocode. I started by using a sentence like "My cat is black" and would go through each step. This made the pseudocode easier to write and build the foundation to how I would write the code. 

As we worked through the code, it was difficult to create certain functions because I knew what I wanted to do, but could not figure out how to create it. Learning about optimization sped up the timing of the code and created cleaner code. When collaborating together on Teams, we came to a consensus when we had similar or overlapping code. Aamna and I shared a lot of code and ideas. Graziano also came up with the idea of helper functions that would optimize the code; however we opted not to put them in because it would extended the run time in Pycharm. I think ideally these would be supplemental functions if we had additional time to perfect them and if I was using a different IDE. 

For future projects, I believe I want to create a similar workflow as I think I am able to grasp the biological concept first, which facilitates how I should approach the code. Another thing to consider here is having patience for larger datasets. Larger datasets take more computational power and will take longer to produce an output. 

### Other members

### Aamna Sheikh

I actually found this project a lot easier to understand conceptually than some of the previous ones, which is surprising because genome assembly in general is a pretty complicated process. Breaking sequencing reads into overlapping k-mers and using those overlaps to build a De Bruijn graph made sense to me. We discussed De Bruijn graphs earlier in the program too which helped because I had a biological application I was already somewhat familiar with so I understood what we were trying to do before getting into the code.

I particularly liked working through the add_edge()/remove_edge() and build_graph_from_reads() functions. Once I understood how each k-mer was split into a prefix and suffix of length k-1, the graph construction clicked and was simple. The toy example helped a lot because I could actually follow how the reads were connected.

The Eulerian walk was probably the most challenging part for me. I understood the idea of traversing every edge exactly once, but recursion took longer to figure out. I had to go through the code a few times to understand why we append nodes after visiting their outgoing edges and why the tour needs to be reversed. I also struggled a little with assemble_contigs(), especially understanding how to choose a starting node and what happens when there are multiple possible paths.

One thing I liked was comparing different approaches to the same problem. I followed the recursive Eulerian walk from Langmead's lecture, while Graziano used an iterative approach with a stack. I also hadn't really thought about what happens when a graph branches and there isn't one clear path to follow. Initially, I was mainly thinking about connecting overlapping sequences, but looking at non-branching paths helped me understand why we can't always assume which direction a sequence should continue.

Overall, I enjoyed this project. I still want to get more comfortable with recursion and understanding how different functions fit together, but I felt more confident with the underlying algorithm this time. I think this project helped me connect the coding side of bioinformatics with the actual biological problem we were trying to solve.

### Graziano Peregrino 

`Issue  - Performance and Scalability`
One of the biggest challenges while working on this project was the program runtime, the toy dataset and the yeast dataset helped me to understand that the scalability can be an issue overtime and it is always good to think how the code will support as the data increases. At first the code took up to 76 and the stats documentation was not correct as output. 

`What I Learned`
I learned about the performance problems on the code and how important it is to optimize the code to support the data increase. The graph construction was the main bottleneck. Repeatedly calling `add_edge()` innermost loop was significantly slower than performing the equivalent graph update directly.

`Next Steps`
As next steps I would like to explore more memory efficient data structures, especially when working with larger genomes such as the mouse genome.

`Issue - Object-Oriented Programming and Code Organization`
During development, my `assemble_contigs()` function became difficult to understand because it contained helper functions inside the main function and several different responsibilities were grouped together.

`What I Learned`
I learned how object-oriented programming can help separate a large algorithm into smaller methods with specific responsibilities. I moved operations such as checking whether a node is one-in/one-out, finding canonical DNA sequences, building non-branching paths, and handling cycles into private methods of the `DeBruijnGraph` class. I also learned how methods such as `@staticmethod` can be used when an operation belongs logically to a class but does not require information from a particular object. I also improved my understanding of object-oriented programming.

`Next Steps`
My next step is to continue learning how bioinformatics software represents very large graphs more efficiently and how production genome assemblers handle sequencing errors, repeated regions, and much larger genomes.


## Generative AI Appendix

The appendix entry must contain:
1. Description of which generative AI was used and its version. Claude - Opus 5

2. The entire prompt that was used to generate the content.
* Can you help me to fix my pseudocode to a latex format and the github? 
* Can you guide me on the process to implement the @staticmethod?
* Can you help to understand the De Bruijn Graph algorithm and the use of the algorithm in bioinformatics?

3. An explanation of how it was used (e.g., to generate ideas).
* It was used as guidance to understand the algorithm and how it is implemented in bioinformatics
* It was also used to understand the output and why the first implementation took longer to complete and what is the best approach to achieve the goal.

4. A justification for why generative AI was used. 
* Used as further explanation on the content and how to use the De Bruijn algorithm in bioinformatics 
