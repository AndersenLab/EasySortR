library(testthat)
library(COPASutils)
library(dplyr)
library(kernlab)

#Test that plate reading and summarization fuctions are 
context("Plate reading")

test_that("Reading plates is not functioning", {
    #Load in a data frame on which to test the summarization function
    testPlate <- readPlate_worms("../testPlate.txt")
    #Test for the right number of and labels for columns 
    expect_identical(colnames(testPlate), c("row", "col", "sort", "TOF", "EXT", "time", "green", "yellow", "red", "norm.EXT", "norm.green", "norm.yellow", "norm.red", "object", "call50", "stage"))
    #Check that the number of rows read in is correct
    expect_equal(nrow(testPlate), 1863)
    #Check that the classes of each column are correct
    expect_identical(as.vector(sapply(testPlate, class)), c("factor", "factor", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "factor", "character"))
    # Test to ensure that normalization to TOF was done correctly (checking multiple rows and all normalizaed traits)
    expect_equal(testPlate[1, 10], testPlate[1, 5]/testPlate[1, 4])
    expect_equal(testPlate[2, 11], testPlate[2, 7]/testPlate[2, 4])
    expect_equal(testPlate[3, 12], testPlate[3, 8]/testPlate[3, 4])
    expect_equal(testPlate[4, 13], testPlate[4, 9]/testPlate[4, 4])
})

context("Plate summarization")

test_that("Summarizing plates is not functioning", {
    testPlate <- readPlate_worms("../testPlate.txt")
    summarizedPlate <- summarizePlate_worms(testPlate)
    
    #Summarized plates should have 96 rows, one for each well
    expect_equal(nrow(summarizedPlate), 96)
    
    # The mean of all of the values in the A3 well should equal that of the value for A3 in the summarized plate data frame
    testPlate_worms_A3 <- testPlate %>% filter(call50 == "object", row == "A", col == 3) %>% select(TOF)
    expect_equal(as.numeric(summarizedPlate[3, "mean.TOF"]), mean(testPlate_worms_A3$TOF))
    
})

