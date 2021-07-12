context("Retrieve frequency")


test_that("retrieve frequencies with pgxseg by filter",{
#    url_e1 <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxseg&filters=43"
    url <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxseg&filters=NCIT:C4323,NCIT:C3697"
    cat(paste("\n trying:",url,"\n"))
    table  <- read.table(url, header=T, sep="\t", na="NA")
    expect_gt(nrow(table),0)
    })

test_that("retrieve frequencies with pgxseg by id",{
    url <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxseg&id=NCIT:C4323"
    cat(paste("\n trying:",url,"\n"))
    table  <- read.table(url, header=T, sep="\t", na="NA")
    expect_gt(nrow(table ),0)
})


test_that("retrieve frequencies with pgxmatrix by filter",{
    #url_e1 <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxmatrix&filters=NCIT:C1234"
    url <- "http://www.progenetix.org/services/intervalFrequencies/?output=pgxmatrix&filters=NCIT:C4323,NCIT:C3697"
    cat(paste("\n trying:",url,"\n"))
    table  <- read.table(url, header=T, sep="\t", na="NA")
    expect_gt(nrow(table),0)
})

