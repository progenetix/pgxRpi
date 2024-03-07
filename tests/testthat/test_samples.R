context("Retrieve samples/individuals")
require(httr)

test_that("retrieve samples with group id",{
    
    url <- "http://progenetix.org/services/sampletable/?filters=NCIT:C3697&datasetIds=progenetix"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve samples with biosample id",{

    url <- "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g&datasetIds=progenetix"
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve individuals with biosample id",{
  url <- "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvjji1&requestEntityPathId=individual&datasetIds=progenetix"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve individuals with individual id",{
  url <- "http://progenetix.org/services/sampletable/?individualIds=pgxind-kftx3565,pgxind-kftx5g4v&requestEntityPathId=individual&datasetIds=progenetix"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve limited individuals with group id",{
  url <- "http://progenetix.org/services/sampletable/?filters=NCIT:C3512&requestEntityPathId=individual&limit=1000&datasetIds=progenetix"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_equal(nrow(result),1000)
})

test_that("retrieve samples with group id in cellz",{
    
    url <- "http://progenetix.org/services/sampletable/?filters=NCIT:C7707&datasetIds=cellz"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

