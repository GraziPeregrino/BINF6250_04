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
## Successes

There were major successes while working on this project. I would like to mention the first one as the scalability of the De Bruijn graph algorithm, running from a small toy dataset to more than 2.4. million yeast sequencing reads. Using the small dataset is the best approach to validate the code and the graph construction. 
I also want to mentioned how the implementation of the program using Object Oriented Programming. Use class to separate responsibilities into methods for the graph construction was an improvement to the code. 
The performance optimization is another point I want to bring. While running the code I called the Add_edge into the build_graph_from_reads function and notice that the code would run correct with the toy dataset, but once we scale it, we notice that the time that the code was running increased. So I removed the add_edge calling and added a conditional instatement inline to the function and that improved the code. 
The statistics result was showing great results and that for sure was a improvement result from the first code we wrote that took a long time and did not provide a good result.

## Struggles

I want to bring few challenges that we faced as we developed the code. The code implementation was showing to add the Add_edge function which was great while working on the toy dataset, but once we increase the dataset the best approach was to add the inline edge insertion instead of calling the add_edge. 
Another struggle we faced was the memory usage, an initial implementation stored forward and reverse graph. Maintaining both adjacency structures required a large amount of memory. The final implementation stores the forward adjacency graph and only the incoming degree of each node, which contains the information required by the config assembly algorithm without storing every incoming adjacency list.
A third challenge was balancing the algorithm described in the assignment with the requirements of a large real dataset. The Eulerian walk works well for demonstrating the required graph-traversal algorithm and for validating small examples. For the complete yeast graph, the final implementation constructs maximal non-branching paths instead. This required understanding why a node with one incoming edge and one outgoing edge represents an unambiguous continuation of a `contig` and why branching nodes should terminate it.


## Personal Reflections
### Group Leader

Group leader's reflection on the project

### Other member - Graziano Peregrino 
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



### Other members' reflections on the project
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
