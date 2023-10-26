context("Retrieve samples/individuals")
require(httr)

test_that("retrieve samples with group id",{
    
    url <- "http://progenetix.org/services/sampletable/?filters=NCIT:C3697"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve samples with biosample id",{

    url <- "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g"
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve individuals with biosample id",{
  url <- "http://progenetix.org/services/sampletable/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvjji1&responseEntityId=individual"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve individuals with individual id",{
  url <- "http://progenetix.org/services/sampletable/?individualIds=pgxind-kftx3565,pgxind-kftx5g4v&responseEntityId=individual"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve limited individuals with group id",{
  url <- "http://progenetix.org/services/sampletable/?filters=NCIT:C3512&responseEntityId=individual&limit=1000"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_equal(nrow(result),1000)
})
