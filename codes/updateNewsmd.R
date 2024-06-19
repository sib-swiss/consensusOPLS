### LICENCE.md
library(usethis)
use_gpl3_license()

### cran-comments.md
use_cran_comments()

### NEWS.md
library(newsmd)

my_news <- newsmd(text=c(paste0("## version 1.0.0"),
                         "", "---", "", "### NEWS.md setup", "",
                         "- added NEWS.md creation with [newsmd]", ""))
my_news$write()
