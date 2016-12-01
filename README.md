# Supplemental material for papers on model propriety of restricted Boltzmann machines and deep learning models

#### Quick Links:

* Shiny [Applets](https://andeek.shinyapps.io/rbms/) illustrating the degeneracy and instability results for small restricted Boltzmann machine models.
* The [degeneracy](https://github.com/andeek/rbm/tree/master/degeneracy) folder contains code for creating the small model explorations in Section 4 of the RBM paper.
* The [4x4 model](https://github.com/andeek/rbm/tree/master/4x4%20model) folder contains code for fitting the 4 visible, 4 hidden model in Section 5 of the RBM paper.
* The [presentations](https://github.com/andeek/rbm/tree/master/presentations) folder contains slide decks and posters that have been presented on this material.
* The [writing](https://github.com/andeek/rbm/tree/master/writing) folder contains fully reproducible versions of both papers.

## ["On the propiety of restricted Boltzmann machines"](https://github.com/andeek/rbm/blob/master/writing/draft.Rmd)

**Authors:** [Andee Kaplan](mailto:ajkaplan@iastate.edu?subject=RBM paper), Daniel Nordman, and Stephen Vardeman  

**Abstract:** A restricted Boltzmann machine (RBM) is an undirected graphical model constructed for discrete or continuous random variables, with two layers, one hidden and one visible, and no conditional dependency within a layer. In recent years, RBMs have risen to prominence due to their connection to deep learning. By treating a hidden layer of one RBM as the visible layer in a second RBM, a deep architecture can be created. RBMs are thought to thereby have the ability to encode very complex and rich structures in data, making them attractive for supervised learning. However, the generative behavior of RBMs is largely unexplored. In this paper, we discuss the relationship between RBM parameter specification in the binary case and the tendency to undesirable model properties such as degeneracy, instability and uninterpretability. We also describe the difficulties that arise in likelihood-based and Bayes fitting of such (highly flexible) models, especially as Gibbs sampling (quasi-Bayes) methods are often advocated for the RBM model structure.

## ["A note on the instability and degeneracy of deep learning models"](https://github.com/andeek/rbm/blob/master/writing/note.Rmd)

**Authors:** [Andee Kaplan](mailto:ajkaplan@iastate.edu?subject=Instability paper), Daniel Nordman, and Stephen Vardeman  

**Abstract:**   A probability model exhibits instability if small changes in a data outcome result in large, and often unanticipated, changes in probability. For correlated data structures found in several application areas, there is increasing interest in predicting/identifying instability. We consider the problem of quantifying instability for general probability models defined on sequences of observations, where each sequence of length N has a finite number of possible outcomes. (A sequence of probability models results indexed by N that accommodates data of expanding dimension.) Model instability is formally shown to occur when a certain log-probability ratio under such models grows faster than N. In this case, a one component change in the data sequence can shift probability by orders of magnitude. Also, as a measure of instability becomes more extreme, the resulting probability models are shown  to tend to degeneracy, placing all their probability on arbitrarily small portions of the sample space. These results on instability apply to large classes of models commonly used in random graphs, network analysis, and machine learning contexts.
