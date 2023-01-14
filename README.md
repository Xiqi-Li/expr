A code pool for RNA expression analysis in Woodman Lab.

## Install

Without SSH key: 

    devtools::install_git(
        url='http://gitlab.mdanderson.edu/XLi23/expr.git',
    )

With SSH key set up:

    devtools::install_git(
        url='git@gitlab.mdanderson.edu:XLi23/expr.git',
        quiet=FALSE
    )
 
  - [How to set up SSH key?](https://docs.gitlab.com/ee/user/ssh.html)

## Usage
View package vignette with `browseVignettes("expr")` in Rstudio console, or [here](http://127.0.0.1:26015/library/expr/doc/expr_vignette.html)

## Development
- This package is documented with [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html). Before pushing a commit, 'Cmd + Shift + B' to generate document, then 'Cmd + Shift + B' to build the package. Keep new functions going!
- Package tests are as good as how you write it. Please kindly report issues under "issue" tab every time an error is encountered.
