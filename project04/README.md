# Introduction
A De Bruijn Graph visualizes the *de novo* genome assembly of short-read sequencing to construct longer consensus sequences (contigs). To construct the contigs, a Eulerian walk (path) is used to connect the nodes (prefixes and suffixs of the k-1 mer) with edges (the overlaps). 

This project is broken down into three challenges: using toy/ sample data, using the yeast genome and attempting the mouse genome. In this code, we implement object-oriented programming by building ```class DeBruijn Graph``` to manipulate sequencing reads to conduct genome assembly by using Eulerian paths. 

```text
# ---------- Setup ----------

INIT(reads, k):
    graph     = {}    # node -> outgoing neighbors
    in_degree = {}    # node -> incoming edge count
    BUILD_GRAPH(reads, k)

# ---------- Utilities ----------

REVERSE_COMPLEMENT(seq):
    RETURN reverse(complement(seq))          # A<->T, C<->G

CANONICAL(seq):
    RETURN min(seq, REVERSE_COMPLEMENT(seq))

ONE_IN_ONE_OUT(node):
    RETURN in_degree[node] == 1 AND |graph[node]| == 1

PATH_TO_SEQUENCE(path):
    RETURN path[0] + last base of each remaining node

# ---------- Graph construction ----------

ADD_EDGE(u, v):
    IF v not in graph[u]:
        add v to graph[u]
        in_degree[v] += 1

REMOVE_EDGE(u, v):
    remove v from graph[u]
    in_degree[v] -= 1
    IF in_degree[v] == 0: delete in_degree[v]

BUILD_GRAPH(reads, k):
    FOR each read:
        FOR each segment in split(read, "N"):
            IF |segment| < k: CONTINUE
            FOR each k-mer in segment:
                ADD_EDGE(prefix(k-mer), suffix(k-mer))

# ---------- Contig assembly ----------

NON_BRANCHING_PATH(u, v, visited):
    path = [u, v]
    WHILE ONE_IN_ONE_OUT(v):
        mark v visited
        v = only neighbor of v
        append v to path
    RETURN path

CYCLE(start, visited):
    path = [start], node = start
    REPEAT:
        mark node visited
        node = only neighbor of node
        append node to path
    UNTIL node == start
    RETURN path

ASSEMBLE_CONTIGS():
    contigs, emitted, visited = [], {}, {}

    FOR each node NOT ONE_IN_ONE_OUT:                 # non-branching paths
        FOR each neighbor:
            EMIT(NON_BRANCHING_PATH(node, neighbor, visited))

    FOR each unvisited node that is ONE_IN_ONE_OUT:   # isolated cycles
        EMIT(CYCLE(node, visited))

    RETURN contigs

EMIT(path):
    seq = PATH_TO_SEQUENCE(path)
    IF CANONICAL(seq) not in emitted:
        add CANONICAL(seq) to emitted
        append seq to contigs

# ---------- Eulerian walk (Hierholzer) ----------

EULERIAN_WALK(start, seed):
    stack, tour = [start], []
    WHILE stack:
        u = top(stack)
        IF graph[u] not empty:
            v = random neighbor of u (seeded)
            REMOVE_EDGE(u, v)
            push v
        ELSE:
            pop u, append u to tour
    RETURN reverse(tour)

# ---------- Output ----------

ASSEMBLY_STATS(contigs):
    lengths = sorted lengths, descending
    report count, total, max, min, mean
    N50 = first length where running sum >= total / 2

WRITE_FASTA(contigs, file):
    FOR i, seq in contigs:
        write ">contig_i" and seq
```
## Output 

Please review the output files named genome_assembly_stats.txt and the FASTA file attached to the Github. 

## Successes

There were major successes while working on this project. 

*Scalability of the De Bruijn graph algorithm*: running from a small toy dataset to more than 2.4. million yeast sequencing reads. Using the small dataset is the best approach to validate the code and the graph construction. 

*Implementation of the program using Object Oriented Programming*: Use class to separate responsibilities into methods for the graph construction was an improvement to the code. 

*Performance optimization*: While running the code I called the Add_edge into the build_graph_from_reads function and notice that the code would run correct with the toy dataset, but once we scale it, we notice that the time that the code was running increased. So I removed the add_edge calling and added a conditional instatement inline to the function and that improved the code. 

*Output Result*: We were successfully able to generate an output for the toy dataset as well as the yeast genome. 

## Struggles

We would like to bring few challenges that we faced as we developed the code. 

*```add_edge()``` function*: The code implementation was showing to add the ```Add_edge``` function which was great while working on the toy dataset, but once we increase the dataset the best approach was to add the inline edge insertion instead of calling ```add_edge```. 

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
I learned about the performance problems on the code and how important it is to optimize the code to support the data increase. The graph construction was the main bottleneck. Repeatedly calling `add_edge()` innermost loop was significantly slower than performing the equivalent graph update directly. I also learned that replacing the complete reverse graph with an `in_degree` dictionary reduced the amount of information that needed to be stored and improved performance.

`Next Steps`
As next steps I would like to explore more memory efficient data structures, especially when working with larger genomes such as the mouse genome.

`Issue - Object-Oriented Programming and Code Organization`
During development, my `assemble_contigs()` function became difficult to understand because it contained helper functions inside the main function and several different responsibilities were grouped together.

`What I Learned`
I learned how object-oriented programming can help separate a large algorithm into smaller methods with specific responsibilities. I moved operations such as checking whether a node is one-in/one-out, finding canonical DNA sequences, building non-branching paths, and handling cycles into private methods of the `DeBruijnGraph` class. I also learned how methods such as `@staticmethod` can be used when an operation belongs logically to a class but does not require information from a particular object. I also improved my understanding of object-oriented programming. Refactoring helper operations into private methods such as `_is_one_in_one_out()`, `_canonical_sequence()`, `_build_non_branching_path()`, and `_build_cycle()` made the class easier for me to understand and explain. 

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
