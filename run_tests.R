#!/usr/bin/env Rscript
library('methods')
library('testthat')
try(test_dir('tests', reporter = 'Summary'))
