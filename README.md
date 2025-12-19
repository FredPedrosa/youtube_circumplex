# Before the Circumplex: The Hierarchical Structure of Affect

This repository contains the data, analysis scripts, and supplementary materials for the study **"The Hierarchical Structure of Affect"**. The research explores the existence of a general factor of Affective Intensity/Salience, contrasting two domains: natural language semantics (via word embeddings) and psychometric self-reports (via PANAS).

## 📝 Abstract

Two-dimensional models of Valence and Arousal dominate the literature on the structure of affect, despite the persistent emergence of a general factor typically treated as a methodological artifact. This study proposes a reinterpretation of this factor as a fundamental dimension of **Affective Intensity/Salience**, organizing affect into a hierarchical structure. 

**Keywords:** Affect, Psychometrics, Psychological Testing, Artificial Intelligence, Word Embeddings.

We investigated this proposal in two domains: 
1. The **semantic structure** of affect, by analyzing word embeddings from a vast linguistic corpus (Study 1).
2. The **psychometric structure** of self-report, using data from the Positive and Negative Affect Schedule - PANAS (Study 2).

A robust analytical approach, including Principal Component Analysis (PCA), Partial Least Squares Structural Equation Modeling (PLS-SEM), and Confirmatory Factor Analysis (CFA), was employed.

## 📂 Repository Structure

### 🛠 Scripts & Notebooks
*   `youtube_scraping.ipynb`: Python notebook used for scraping YouTube comments, text preprocessing, and generating word embeddings.
*   `Circumplex_Suplementar.Rmd`: Main R script containing the statistical pipeline (PCA, Parallel Analysis, PLS-SEM, CFA, and convergence plots).

### 📊 Datasets
*   `embeddings_circumplex.csv`: Matrix of mean contextual embeddings for the 45 affective keywords (Study 1).
*   `dataset_final_classificado.csv`: The final dataset used for the psychometric analyses in Study 2.
*   `comentarios_classificados_final.csv`: Processed comments with sentiment classification.
*   `ANEW_BR.xlsx` & `ANEW2025.xlsx`: Affective Norms for Portuguese (lexical validation).
*   `comentarios_youtube_balanceado.csv` & `robust_comments`: Raw and balanced datasets from the YouTube corpus.

### 📜 Supplementary Documents
*   `Acircomplexmodelofaffect.pdf`: Russell's paper on the Circumplex model.
*   `CircumplexClean.pdf`: Syntax in Brazilian portuguese.

## 🚀 How to Replicate

### Prerequisites
To run the R analysis, you will need the following packages:
`pacman`, `tidyverse`, `psych`, `lavaan`, `seminr`, `ggplot2`, `patchwork`, `ggrepel`.

### Steps
1. Clone the repository.
2. Open `Circumplex_Suplementar.R` in RStudio.
3. Ensure all data files (`.csv`) are in the same directory.
4. Run the script to generate the models and figures.

ps.: The PANAS data is not public avaiable. 

## ✍️ Authors

*   **Frederico G. Pedrosa** - *Federal University of Minas Gerais (UFMG)* - [ORCID](https://orcid.org/0000-0002-0682-0734)

---

## 🎓 How to Cite this Code

If you use the code or analysis pipeline from this repository, please cite it as follows:

**APA Style:**
Pedrosa, F. G. (2025). *The Hierarchical Structure of Affect: Data and Analysis Code* [GitHub repository]. https://github.com/FredPedrosa/youtube_circumplex

**BibTeX:**
```bibtex
@manual{pedrosa2025affect,
  title  = {The Hierarchical Structure of Affect: Data and Analysis Code},
  author = {Pedrosa, Frederico G.},
  year   = {2025},
  url    = {https://github.com/FredPedrosa/youtube_circumplex},
  note   = {GitHub repository}
}
```
