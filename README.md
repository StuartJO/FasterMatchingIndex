# FasterMatchingIndex

This is code for computing the matching index much faster. This is of particular relevance for people using generative network models

## What is the matching index?

The matching index, as typically used, is a measure that quantifies the similarity in two nodes connectivity profiles*. The index is calculated on an adjacency matrix $A$, and is mathematically defined as:

(1) $$M_{ij} = \frac{|\Gamma_{i}-{j}\cap \Gamma_{j}-{i}|}{|\Gamma_{i}-{j}\cup \Gamma_{j}-{i}|}$$

where $\Gamma_{i}-{j}$ is the set of neighbours $N$ of node $i$ excluding node $j$ (if it is at all connected to node $j$). 

<sub>* See the [section below](#what-did-you-mean-when-you-said-the-matching-index-may-not-measure-what-I-think-it-does) discussing how the matching index may not measure what we usually think it does...</sub>

## So what is the problem?

The matching index has often been found to give the best performance among topological generative network models (for example, see [here](https://doi.org/10.1016/j.neuroimage.2015.09.041) [here](https://www.science.org/doi/10.1126/sciadv.abm6127) [here](https://doi.org/10.1016/j.neuroimage.2020.117510) [here](https://doi.org/10.1038/s41467-021-24430-z) and [here](https://doi.org/10.1002/dev.22405)), meaning that it is the model people are most interested in running.

_HOWEVER_

It takes agesssssssssssssssssssssss to run (in my experience I could run ever single other topological model in the time it takes to run the matching model)

This likely makes you sad

But what if I told you there was a way to make it (somewhat) better

## I'm interested now, why can it be made faster?

To understand why it can be made faster, we need to think about how it is actually calculated. To do so, we need to talk more maths (sorry).

We can also calculate the matching index as

(2) $$M_{ij} = \frac{2(N_{ij}-A_{ij})}{k_{i}+k_{j}-2A_{ij}}$$

which is just the number of neighbours $i$ and $j$ share multiplied by two (whilst excluding themselves as neighbours of the other, thats what the $-A_{ij}$ is for), then dvided by the summed degree $k$ of nodes $i$ and $j$ (whilst ignoring any connection that may exist between nodes $i$ and $j$).

When written this way it is trivial* to see how the calculation could easily be done progromatically. If you look at the original code provided in the [BCT](https://sites.google.com/site/bctnet/), you'll notice it is actually calculating it the second way and not the first. However it is looping over all the nodes to calculate it. We can actually forgo any loops when calculating this measure resulting in a considerable speed up in processing speed

<sub>* I'vew always wanted to say this haha. The reason it is trivial is because we can take advatage of matrix operations to compute this value (as the number of neighbours can be easily computed by taking the square of the matrix, and degree can very easy be obtained by just taking the sum)</sub>

## So how much faster is it?

First, lets compare calculating the matching index on networks with different numbers of nodes:


You can see that the old way takes significantly longer when computing the index than the new way, and this time only increases as the network gets larger! So the benefits of using the new code gets considerably better as the network gets larger

## Yes that's neat Stuart, but what about its use in generative network models?

An important thing to note about generative network models is they _iteratively_ add edges. As I alluded to before, the new code largely benefits because it can calculate everything in one hit instead of needing to loop over all the nodes. In the old generative model node, at each iteration the loop to calculate the matching index only runs over nodes who will be affected by the newly added edge (i.e., the loop very likely doesn't need to be run over all nodes). 



You'll notice as the number of edges increases, the improvement of the new codes lessens (but it is still better). The reason for this is simple. When the generative model starts, it initalises a topology matrix. In other words, at very beginning it needs to run the matching index for the entire network. As I demonstrated earlier, this is wehere the new code will result in fairly massive improvements. If we remove this factor however (by initalising the topology matrix the same way in the new and old code), we still see a benefit!



If we directly compare the timings of the code as they progress over iterations of the model (timings have been normalised in the plot below), we can see that as the old code progresses it gets slower, this is because there will be more neighbours to iterate over at later steps (because the network is forming)


But in summary, no matter how you cut it, the new code is faster

## I want to see this with my own eyes

Easy! Just run the script I wrote to demonstrate this in MATLAB

```
matchingSpeedTest.m
```

Note it does use the Parellel Computing Toolbox, just to speed things up a bit. It takes a bit over 30 minutes to run everything (using the Toolbox) on a i7 6700k FYI 

## What did you mean when you said the matching index may not measure what I think it does?

The matching index is commonly considered a normalised measure of the overlap of two nodes' neighbourhoods. Conceptually, we would understand this as meaning the number of shared neighbours divided by the total unique neighbours, as the original mathematical definition shown in Equation 1 suggests.

Yet consider the network below:

![An example network to illustrate how the matching index is calculated. The matching index is being calculated between the two blue nodes (labelled i and j). They share three neighbours (red nodes). Their combined neighbourhood size is 10 (red nodes plus grey nodes). They also share six connections to neighbouring nodes (red edges), and have 13 total connections (red edges plus black edges). Note that by convention, direct connections between the nodes of interest (in this case the blue nodes) are discounted](./images/MatchingDemo.svg)

Node $i$ and $j$ share three neighbours (red nodes). The combined total of (unique) neighbours $i$ and $j$ have is 10 (red and grey nodes), so we would assume the matching index of these nodes is $\frac{3}{10}$.

We would _technically_ be incorrect however, or rather we would have a different answer to what the code (both old and new) provides.

As mentioned above, the matching index is similartity in the _connectivity profiles_ of two nodes. This means the matching index is actually calculated as the number of connections a pair of nodes have to same neighbours, over the total number of connections those nodes have. So in the example above, as node $i$ and $j$ share three neighbours they have six connections in common (red edges). They have 13 connections in total (red edges plus the black edges, not the connection between them is excluded by convention), so the matching index is $\frac{6}{13}$. I would say that this definition isn't exactly consistent with Equation 1 (in my opinion, although I am not overly across set theory), but it is exactly how Equation 2 is done. However I would argue that this definition (which I shall call the connectivity profiles definition) isn't the most intuitive. We can change Equation 2 to be more consistent with the intuitive conceptualisation (which I shall call the normalised overlapping neighbourhood definition) by doing the following

(3) $$M_{ij} = \frac{N_{ij}-A_{ij}}{k_{i}+k_{j}-2A_{ij}-N_{ij}}$$

## Do these differing conceptualisations affect anything?

No not really. They will give different results as I showed above, but so long as the same calculation is being used throughout the paper, it should be ok. May affect how you discuss and interpret this measure though. The code I provided does the connectivity profiles definition by default, but does allow for the normalised overlapping neighbourhood definition to be done as well (note this is only done for the "matching.m" function, all the generative modelling functions at current can only use the connectivity profiiles formulation)

## I would like to incorporate these new ways of computing the matching index into my own code, is there an easy way to do this?

Good news! I have written code which allows you to do this! The inputs and outputs should be very similar to what the BCT/Betzel implementation used (and is consistent with the code I wrote for my [paper](https://www.science.org/doi/10.1126/sciadv.abm6127)) 

## Wow this is so great, you're amazing Stuart

Thank you

## Who should I cite for this implementation?

Please just reference this GitHub, there is no paper attached to this new code (however I really won't object to you [citing my work in this space at all](https://www.science.org/doi/10.1126/sciadv.abm6127)).  

## Is there a Python version?

Ha no. I leave it as an exercise to the reader to figure that one out.