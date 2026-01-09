#!/usr/bin/Rscript

library("webshot")

# Set up variable to control command line arguments
args <- commandArgs(TRUE)
# Set up variable for filename
input_filename <- args[1]

png_height <- 3500
png_width <- 6000

filename_png <- gsub(".html", ".png", input_filename)
webshot(input_filename, filename_png, vwidth = png_width, vheight = png_height)
