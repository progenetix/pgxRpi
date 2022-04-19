context("Retrieve variants")

#url_e1 <- "http://progenetix.org/cgi/bycon/beaconServer/variants.py?biosampleIds=1"
#url_e2 <- "https://progenetix.org/cgi/bycon/beaconServer/variants.py?biosampleIds=pgxbs-kftvh98e"
url <- "http://progenetix.org/beacon/variants/?biosampleIds=pgxbs-kftvh94d,pgxbs-kftvh94g,pgxbs-kftvh972"

test_that("retrieve variants with JSON",{
        cat(paste("\n trying:",url,"\n"))
        table <- rjson::fromJSON(file = url)
        table <- lapply(table$response$resultSets[[1]]$results,unlist)
        table <- as.data.frame(dplyr::bind_rows(table))
        expect_gt(nrow(table),0)
})


test_that("retrieve variants with pgxseg",{
    url <- paste0(url,"&output=pgxseg")
    cat(paste("\n trying:",url,"\n"))
    table <- read.table(url, header = T, sep="\t")
    expect_gt(nrow(table),0)
    expect_gt(ncol(table),1)
})
