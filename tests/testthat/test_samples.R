context("Retrieve samples")
library(httr)

test_that("retrieve samples with group id",{
    
    #url_e1 <- "http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=NCIT:C1234&output=table"
    #url_e2 <- "http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=N&output=table"
    url <- "http://progenetix.org/beacon/biosamples/?filters=NCIT:C3697&output=table"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_error(read.table(url, skip = 1, stringsAsFactors = FALSE, sep = "\t",fill=TRUE),NA)
    

})

test_that("retrieve samples with biosample id",{

    url <- "http://progenetix.org/beacon/biosamples/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972&output=datatable"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_error(read.table(url, skip = 1, stringsAsFactors = FALSE, sep = "\t",fill=TRUE),NA)

})