# cengenAnnotate: Automated CeNGEN-Based Annotation of Seurat Clusters

Automates annotation of C. elegans single-cell RNA-seq clusters by
scoring Seurat cluster marker genes against CeNGEN reference expression
data. Computes a dot-plot-equivalent coherence score for each candidate
cell type from unthresholded expression and percent-expressing values,
and classifies clusters as confidently annotatable or in need of manual
review, replacing the manual workflow of pasting marker genes into the
CeNGEN website and visually reading the resulting dot plot. Reference
datasets are discovered and fetched at runtime from a 'pins' board so
new CeNGEN dataset releases (developmental stage and sex variants)
become available without package changes.

## See also

Useful links:

- <https://github.com/cengenproject/cengenAnnotate>

- Report bugs at
  <https://github.com/cengenproject/cengenAnnotate/issues>

## Author

**Maintainer**: Marc Hammarlund <marc.hammarlund@yale.edu>
