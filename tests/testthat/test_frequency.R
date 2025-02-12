context("Retrieve frequency")
require(httr)

test_that("retrieve frequencies by filter",{
    url <- "http://www.progenetix.org/services/intervalFrequencies/?filters=NCIT:C4323,pgx:icdom-85003"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(length(result$response$results),2)
    expect_equal(length(result$response$results[[1]]$intervalFrequencies),3106)
})


# test_that("retrieve frequencies by filter in cellz",{
#     url <- "https://cancercelllines.org/services/intervalFrequencies/?filters=NCIT:C5228,pgx:icdom-87203"
#     cat(paste("\n trying:",url,"\n"))
#     r <- GET(url)
#     expect_equal(http_type(r), "text/plain")
#     expect_no_error(result  <- read.table(url, header=T, sep="\t"))
#     expect_equal(dim(result),c(2,6213))
# })