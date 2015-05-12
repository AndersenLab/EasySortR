library(testthat)
library(COPASutils)
library(dplyr)
library(kernlab)
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
})

context("Plate summarization")

test_that("Summarizing plates is not functioning", {
    testPlate <- readPlate_worms("../testPlate.txt")
    summarizedPlate <- summarizePlate_worms(testPlate)
    
})

