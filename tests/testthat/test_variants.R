context("Retrieve variants")
require(httr)
require(dplyr)
#url_e1 <- "http://progenetix.org/cgi/bycon/beaconServer/variants.py?biosampleIds=1"
#url_e2 <- "https://progenetix.org/cgi/bycon/beaconServer/variants.py?biosampleIds=pgxbs-kftvh98e"
url <- "http://progenetix.org/beacon/variants/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972"
url_2 <- "http://progenetix.org/beacon/biosamples/pgxbs-kftvh94d/g_variants"
test_that("retrieve variants with JSON",{
        cat(paste("\n trying:",url_2,"\n"))
        result <-  content(GET(url_2))
        expect_equal(result$responseSummary$exists,TRUE)
        table <- lapply(result$response$resultSets[[1]]$results,unlist)
        table <- as.data.frame(bind_rows(table))
        expect_gt(nrow(table),0)
})


test_that("retrieve variants with pgxseg",{
    url <- paste0(url,"&output=pgxseg")
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=T))
    expect_gt(nrow(result),0)
    expect_gt(ncol(result),0)
})
