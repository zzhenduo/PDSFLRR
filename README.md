# PDSFLRR
The code is from the paper:PDSFLRR:Low-Rank Representation with Projection Distance and Sparsity Constraints for Clustering Single-Cell RNA Sequencing Data.
This is abstract: The development of single-cell RNA sequenc-
ing (scRNA-seq) technology has provided essential re-
sources for identifying cellular heterogeneity and diver-
sity. Clustering is one of the key techniques for analyzing
scRNA-seq data. However, due to the high dimensionality,
sparsity, and dropout issues of such data, existing cluster-
ing methods still need improvement. This paper proposes
a novel clustering method, PDSFLRR, for scRNA-seq data.
Specifically, to overcome the limitation of traditional low-
rank representation (LRR) models that fix the dictionary to
the data itself, PDSFLRR employs a projection-based strat-
egy for dictionary updating. Additionally, it incorporates
projection distance constraints to enhance the model’s
ability to learn local structures and ensure the validity of
the dictionary matrix. Recognizing the sparse nature of
scRNA-seq data, PDSFLRR introduces a sparse constraint.
Moreover, it replaces the nuclear norm with the Frobe-
nius norm to avoid extensive singular value decomposition
(SVD) computations. The Alternating Direction Method of
Multipliers algorithm (ADMM) is used to solve the objective
function. Finally, we compared the proposed PDSFLRR
method with nine popular methods across twelve datasets.
Experimental results demonstrate the superiority of PDS-
FLRR.
