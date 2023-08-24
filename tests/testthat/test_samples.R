context("Retrieve samples/individuals")
require(httr)

test_that("retrieve samples with group id",{
    
    #url_e1 <- "http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=NCIT:C1234&output=table"
    #url_e2 <- "http://progenetix.org/cgi/bycon/beaconServer/biosamples.py?filters=N&output=table"
    url <- "http://progenetix.org/beacon/biosamples/?filters=NCIT:C3697&output=table"

    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve samples with biosample id",{

    url <- "http://progenetix.org/beacon/biosamples/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972&output=datatable"
    cat(paste("\n trying:",url,"\n"))
    r <- GET(url)
    expect_equal(http_type(r), "text/tsv")
    expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE,header = T))
    expect_gt(nrow(result),0)
})

test_that("retrieve individuals with biosample id",{
  url <- "https://progenetix.org/beacon/individuals/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvjji1&output=datatable"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve individuals with individual id",{
  url <- "https://progenetix.org/beacon/individuals/?individualIds=pgxind-kftx3565,pgxind-kftx5g4v&output=datatable"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_gt(nrow(result),0)
})

test_that("retrieve limited individuals with group id",{
  url <- "https://progenetix.org/beacon/individuals/?filters=NCIT:C3512&limit=1000&output=datatable"
  cat(paste("\n trying:",url,"\n"))
  r <- GET(url)
  expect_equal(http_type(r), "text/tsv")
  expect_no_error(result <- read.table(url, stringsAsFactors = FALSE, sep = "\t",fill=TRUE, header=T))
  expect_equal(nrow(result),1000)
})
