---
project: MSOLV-fortran
summary: MSOLV-fortran -- A modern Fortran API for Sparse Direct Solvers, as part of <br/><br/> ![MUSST_img](media/MUSST_long_petit.png)
author:  Arthur Francisco - Noël Brunetière
email:   arthur.francisco@univ-poitiers.fr
website: https://www.pprime.fr/?q=fr/recherche-scientifique/d3/mecanique-des-interfaces-lubrifiees
github: https://github.com/Arthur-Francisco
project_github: https://github.com/TRIBO-Pprime/MSOLV
include:    ./inc
src_dir:    ./src
exclude: hsl_common90.f90
         hsl_ddeps90.f90
         hsl_ma48d.f90
output_dir: ./docs
media_dir:  ./img
favicon:    ./img/logo.png
docmark:        !
docmark_alt:    #
predocmark:     >
predocmark_alt: <
display: public
         protected
         private
source: true
graph:  true
search: true
sort: src
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc(anchorlink=False)
               markdown.extensions.smarty(smart_quotes=False)
---

-----------------
{!README.md!}

