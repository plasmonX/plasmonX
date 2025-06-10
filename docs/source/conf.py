# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'plasmonX'
copyright = '2025 The Authors'
author = 'Tommaso Giovannini, Luca Bonatti, Stefano Corni, Pablo Grobas Illobre, Piero Lafiosca, Luca Nicoli, Chiara Cappelli'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

# -- General configuration ------------------------------------------------

# -- General configuration ------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',  # Documentazione automatica del codice
    'sphinx.ext.napoleon',  # Supporto per Google-style e NumPy-style docstrings
    'sphinx.ext.viewcode',  # Per aggiungere link al codice sorgente
]

# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'  # Tema di Read the Docs

# Aggiungi un logo personalizzato (opzionale)
html_logo = 'images/logo_name.png'

# Opzioni per il tema sphinx_rtd_theme
html_theme_options = {
    'logo_only': True,
    'navigation_depth': 3,
    'collapse_navigation': False,
}

# -- Options for HTMLHelp output ------------------------------------------
htmlhelp_basename = 'MyProjectDoc'

# -- Options for LaTeX output ---------------------------------------------
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '11pt',
    'inputenc': '',  # disattiva inputenc (per modern TeX)
'utf8extra': '',
    'maketitle': r'''
\begin{titlepage}
    \centering
    \vspace*{3cm}

    \includegraphics[width=0.4\textwidth]{logo_latex.png}\par
    \vspace{1cm}

    {\Huge\bfseries plasmonX\par}
    \vspace{1cm}

    {\Large Tommaso Giovannini, Luca Bonatti, Stefano Corni, Pablo Grobas Illobre,
     Piero Lafiosca, Luca Nicoli, Chiara Cappelli\par}
    \vspace{2cm}

    {\large \today\par}
    \vfill

    \thispagestyle{empty}
\end{titlepage}

\setcounter{tocdepth}{2}
\tableofcontents
\clearpage
''',    
    'preamble': r'''
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{hyperref}
\usepackage{titlesec}
\usepackage{helvet}
\usepackage[utf8]{inputenc}
\usepackage{newunicodechar}
\newunicodechar{ω}{$\omega$}
\newunicodechar{μ}{$\mu$}
\newunicodechar{ε}{$\varepsilon$}
\newunicodechar{🎸}{}
\newunicodechar{💀}{}
\newunicodechar{🚪}{}
\renewcommand{\familydefault}{\sfdefault}
\titleformat{\section}{\Large\bfseries}{\thesection}{1em}{}
\titleformat{\subsection}{\large\bfseries}{\thesubsection}{1em}{}
''',
    'tableofcontents': '',
    'figure_align': 'H',
}

latex_additional_files = ['images/logo_latex.png']

latex_documents = [
    ('index', 'plasmonX.tex', 'plasmonX',
     'Tommaso Giovannini et al.', 'manual'),
]

# -- Options for manual page output -------------------------------------
man_pages = [
    ('index', 'MyProject', 'My Project Documentation',
     ['Author Name'], 1)
]

# -- Options for Texinfo output -----------------------------------------
texinfo_documents = [
    ('index', 'MyProject', 'My Project Documentation',
     'Author Name', 'MyProject', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
    'style.css',
]

