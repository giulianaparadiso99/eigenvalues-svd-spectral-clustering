# Eigenvalue Methods, Spectral Clustering and PCA

Homework projects for **Computational Linear Algebra for Large Scale Problems** – MSc in Mathematical Engineering, Politecnico di Torino (A.Y. 2024/2025).

**Authors:** Lucia Ghezzi, Giuliana Paradiso

The repository collects three homework projects. The first two are implemented in MATLAB with self-written numerical routines; the third is a Python notebook.

| Project | Topics | Language | Report |
|---|---|---|---|
| [Methods for computing eigenvalues](#1-methods-for-computing-eigenvalues) | Power and inverse power methods, shifting, deflation, QR, Lanczos, SVD | MATLAB | [HW3_Ghezzi_Paradiso.pdf](eigenvalue-methods/HW3_Ghezzi_Paradiso.pdf) |
| [Spectral clustering](#2-spectral-clustering) | k-nearest neighbour similarity graphs, graph Laplacian, eigenvector-based clustering | MATLAB | [HW_SC_Ghezzi_Paradiso.pdf](spectral-clustering/HW_SC_Ghezzi_Paradiso.pdf) |
| [PCA and k-means](#3-pca-and-k-means-clustering) | Data cleaning, scaling, principal component analysis, k-means, cluster evaluation | Python | [HWpca_Ghezzi_Paradiso.pdf](pca/HWpca_Ghezzi_Paradiso.pdf) |

## 1. Methods for computing eigenvalues

Implementation and comparison of iterative methods for computing eigenvalues and eigenvectors of square matrices, and singular values and vectors of rectangular matrices.

**Methods**

- **Power method** for the dominant eigenpair and **inverse power method** for the smallest one.
- **Shifting** (for symmetric matrices) and **deflation** with Householder reflections, to compute the subsequent eigenpairs.
- **QR method**: reduction to upper Hessenberg form, QR iterations, and shifted QR.
- **Lanczos method**: reduction of symmetric matrices to tridiagonal form.
- **Singular value decomposition** through bidiagonalisation.

**Evaluation.** Each method is tested on random matrices with prescribed properties. For each method the report gives convergence, correctness of the eigenvalues (tolerance 1e-7) and eigenvectors (tolerance 1e-4), execution time, number of iterations, experimental order of convergence and residual ‖Av − λv‖.

**Main results**

- All methods compute the eigenvalues accurately; the eigenvectors are less accurate.
- The power and inverse power methods are highly accurate for a single eigenpair. When they are iterated to obtain all the eigenpairs, the accuracy of the eigenvectors decreases, especially with deflation.
- Lanczos is the most accurate method for symmetric positive definite matrices, but it applies only to that class. QR applies to a wider range of matrices, at the cost of less accurate eigenvectors.

## 2. Spectral clustering

Self-implemented spectral clustering on two two-dimensional datasets, `Circle.mat` and `Spiral.mat`. The Spiral dataset also contains the true cluster labels, which are used for validation.

**Pipeline**

1. Gaussian similarity matrix and **k-nearest neighbour similarity graph** (k = 10, 20, 40), with symmetrised adjacency matrix W.
2. Degree matrix D, **Laplacian** L = D − W and normalised symmetric Laplacian L_sym, stored as sparse matrices.
3. Number of connected components of the graph, and choice of the number of clusters M from the smallest eigenvalues of L.
4. Computation of the M smallest eigenpairs of L with a **self-implemented function combining the inverse power method and deflation**.
5. k-means on the rows of the eigenvector matrix, and assignment of the original points to the clusters.

**Comparisons**

- L versus L_sym, for different values of k.
- The self-implemented algorithm versus MATLAB's built-in spectral clustering, k-means and k-medoids.
- Different numbers of clusters (M = 2 to 12).

**Main results**

- M = 3 is the most suitable number of clusters for both datasets.
- The best results are obtained with k = 10, where the self-implemented method matches MATLAB's spectral clustering and reproduces the true clusters of the Spiral dataset. Accuracy decreases as k grows.
- Spectral clustering separates the non-convex structures of both datasets much better than centroid-based methods (k-means, k-medoids).

## 3. PCA and k-means clustering

Exploratory analysis of a survey dataset on the interests, habits, personality traits and demographics of young people (`responses_hw.csv`, 672 respondents). The features analysed are a subset selected from the student IDs.

**Workflow**

1. **Data preparation:** analysis of missing values, mode imputation for the features and removal of the rows with missing labels, label and one-hot encoding of categorical variables.
2. **Scaling and variance:** comparison of the raw data with StandardScaler and MinMaxScaler, and analysis of the cumulative explained variance.
3. **Dimensionality reduction:** PCA on the scaled datasets, interpretation of the principal components, and score graphs coloured by demographic labels.
4. **Clustering:** k-means in the space of the principal components, with interpretation of the centroids.
5. **Cluster evaluation:** external evaluation against the labels (homogeneity, completeness, adjusted Rand index) and internal evaluation with silhouette scores.

**Main results**

- The first five principal components explain less than 33% of the variance, so the reduction loses a significant amount of information.
- The score graphs show no clear structure, except for a partial separation by gender.
- k-means finds poorly separated clusters (low silhouette scores), which confirms that the dataset has no evident group structure.
- The report discusses possible improvements: more principal components, or a more focused selection of features.

## Repository structure

```
.
├── README.md
├── eigenvalue-methods/
│   ├── HW3_Ghezzi_Paradiso.pdf
│   └── MATLAB scripts and functions
├── spectral-clustering/
│   ├── HW_SC_Ghezzi_Paradiso.pdf
│   └── MATLAB scripts and functions
└── pca/
    ├── HWpca_Ghezzi_Paradiso.pdf
    └── HWpca_Ghezzi_Paradiso.ipynb
```

## Requirements

- **MATLAB** for the eigenvalue and spectral clustering projects (the Statistics and Machine Learning Toolbox is needed for `kmeans`, `kmedoids` and `spectralcluster`).
- **Python 3** with `numpy`, `pandas`, `matplotlib`, `scikit-learn` and `pyyaml` for the PCA notebook.

The datasets (`Circle.mat`, `Spiral.mat`, `responses_hw.csv`) were provided by the course and are not included in this repository.
