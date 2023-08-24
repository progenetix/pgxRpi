context("Retrieve frequency")
require(httr)

test_that("retrieve frequencies with pgxseg by filter",{
#    url_e1 <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxseg&filters=43"
    url <- "http://www.progenetix.org/services/intervalFrequencies/?filters=NCIT:C4323,pgx:icdom-85003&output=pgxfreq"
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result  <- read.table(url, header=T, sep="\t"))
    expect_gt(nrow(result),0)
    expect_equal(length(unique(result[,1])),2)
    })


test_that("retrieve frequencies with pgxmatrix by filter",{
    #url_e1 <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxmatrix&filters=NCIT:C1234"
    url <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxmatrix&filters=NCIT:C4323,pgx:icdom-85003"
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result  <- read.table(url, header=T, sep="\t"))
    expect_equal(dim(result),c(2,6213))
})

