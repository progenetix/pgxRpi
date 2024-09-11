context("Retrieve samples/individuals")
require(httr)

# "http://progenetix.org/services/sampletable/?filters=NCIT:C3697&responseEntityPathId=biosamples&datasetIds=progenetix"
# "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g&responseEntityPathId=biosamples&datasetIds=progenetix"
# "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvjji1&responseEntityPathId=individuals&datasetIds=progenetix"
# "http://progenetix.org/services/sampletable/?individualIds=pgxind-kftx3565,pgxind-kftx5g4v&responseEntityPathId=individuals&datasetIds=progenetix"
# "http://progenetix.org/services/sampletable/?filters=NCIT:C7707&responseEntityPathId=biosamples&datasetIds=cellz"

test_that("retrieve samples with group id",{  
    url <- "http://progenetix.org/beacon/biosamples?filters=NCIT:C3697"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_gt(result$responseSummary$numTotalResults,0)
})

test_that("retrieve samples with biosample id",{
    url <- "http://progenetix.org/beacon/biosamples?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_gt(result$responseSummary$numTotalResults,0)
})

test_that("retrieve individuals with biosample id",{
    url <- "http://progenetix.org/beacon/individuals?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvjji1"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_gt(result$responseSummary$numTotalResults,0)
})

test_that("retrieve individuals with individual id",{
    url <- "http://progenetix.org/beacon/individuals?individualIds=pgxind-kftx3565,pgxind-kftx5g4v"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_gt(result$responseSummary$numTotalResults,0)
})

test_that("retrieve limited individuals with group id",{
    url <- "http://progenetix.org/beacon/individuals?filters=NCIT:C3512&limit=5"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_equal(length(result$response$resultSets[[1]]$results),5)
})

test_that("retrieve samples with group id in cellz",{
    url <- "https://cancercelllines.org/beacon/biosamples?filters=NCIT:C7707"
    cat(paste("\n trying:",url,"\n"))
    result <-  content(GET(url))
    expect_equal(result$responseSummary$exists,TRUE)
    expect_gt(result$responseSummary$numTotalResults,0)
})

