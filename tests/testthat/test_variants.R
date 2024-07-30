context("Retrieve variants")
require(httr)
require(dplyr)
url <- "http://progenetix.org/beacon/biosamples/pgxbs-kftvh94d/g_variants?datasetIds=progenetix"
url_2 <- "https://progenetix.org/services/pgxsegvariants/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972&datasetIds=progenetix"

url_3 <- "https://progenetix.org/services/samplematrix/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972&datasetIds=progenetix"
url_4 <- "https://progenetix.org/services/samplematrix/?individualIds=pgxind-kftx3565,pgxind-kftx5g4v&datasetIds=progenetix"
url_5 <- "https://progenetix.org/services/samplematrix/?filters=pgx:icdom-88503&datasetIds=progenetix"

url_6 <- "https://progenetix.org/beacon/analyses/?output=cnvstats&biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972&datasetIds=progenetix"
url_7 <- "https://progenetix.org/beacon/analyses/?output=cnvstats&individualIds=pgxind-kftx3565,pgxind-kftx5g4v&datasetIds=progenetix"
url_8 <- "https://progenetix.org/beacon/analyses/?output=cnvstats&filters=pgx:icdom-88503&datasetIds=progenetix"
url_9 <- "http://progenetix.org/beacon/biosamples/cellzbs-0821E6df/g_variants?datasetIds=cellz"
test_that("retrieve variants with JSON",{
        cat(paste("\n trying:",url,"\n"))
        result <-  content(GET(url))
        expect_equal(result$responseSummary$exists,TRUE)
        table <- lapply(result$response$resultSets[[1]]$results,unlist)
        table <- as.data.frame(bind_rows(table))
        expect_gt(nrow(table),0)
})


test_that("retrieve variants with pgxseg",{
    cat(paste("\n trying:",url_2,"\n"))
    r <- GET(url_2)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result <- read.table(url_2, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=T))
    expect_gt(nrow(result),0)
    expect_gt(ncol(result),5)
})

test_that("retrieve pgxmatrix variant with biosample id",{
    cat(paste("\n trying:",url_3,"\n"))
    r <- GET(url_3)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result <- read.table(url_3, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=T))
    expect_gt(nrow(result),0)
    expect_gt(ncol(result),6212)
})

test_that("retrieve pgxmatrix variant with individual id",{
    cat(paste("\n trying:",url_4,"\n"))
    r <- GET(url_4)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result <- read.table(url_4, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=T))
    expect_gt(nrow(result),0)
    expect_gt(ncol(result),6212)
})

test_that("retrieve pgxmatrix variant with filters",{
    cat(paste("\n trying:",url_5,"\n"))
    r <- GET(url_5)
    expect_equal(http_type(r), "text/plain")
    expect_no_error(result <- read.table(url_5, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header=T))
    expect_gt(nrow(result),0)
    expect_gt(ncol(result),6212)
})

test_that("retrieve fraction variant with biosample id",{
    cat(paste("\n trying:",url_6,"\n")) 
    result <-  content(GET(url_6))
    expect_equal(result$responseSummary$exists,TRUE)
    table <- lapply(result$response$resultSets[[1]]$results,unlist)
    table <- as.data.frame(bind_rows(table))
    expect_gt(nrow(table),0)
})

test_that("retrieve fraction variant with individual id",{
    cat(paste("\n trying:",url_7,"\n"))
    result <-  content(GET(url_7))
    expect_equal(result$responseSummary$exists,TRUE)
    table <- lapply(result$response$resultSets[[1]]$results,unlist)
    table <- as.data.frame(bind_rows(table))
    expect_gt(nrow(table),0)
})

test_that("retrieve fraction variant with filters",{
    cat(paste("\n trying:",url_8,"\n"))
    result <-  content(GET(url_8))
    expect_equal(result$responseSummary$exists,TRUE)
    table <- lapply(result$response$resultSets[[1]]$results,unlist)
    table <- as.data.frame(bind_rows(table))
    expect_gt(nrow(table),0)
})

test_that("retrieve variants with JSON in cellz",{
        cat(paste("\n trying:",url_9,"\n"))
        result <-  content(GET(url_9))
        expect_equal(result$responseSummary$exists,TRUE)
        table <- lapply(result$response$resultSets[[1]]$results,unlist)
        table <- as.data.frame(bind_rows(table))
        expect_gt(nrow(table),0)
})

