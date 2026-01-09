#!/usr/bin/Rscript

# necessary libraries
library("networkD3")
library("tibble")
library("htmlwidgets", warn.conflicts = FALSE)

# Set up variable to control command line arguments
args <- commandArgs(TRUE)
# Set up variable for filename
input_filename <- args[1]

data_links = read.csv(input_filename, stringsAsFactors = FALSE)
png_height <- 3500
png_width <- 6000
font_size1 <- 50
nodepadding_size <- 30

taxonomy <- as.character(data_links$label)
nodes <- data.frame("name" = taxonomy)
# this part is necessary for the OnRender function to display names on the left side
nodes <- tibble(name = nodes$name, target = grepl(':', name))
nodes <- as.data.frame(nodes)

source1 <- as.numeric(data_links$link1[2:nrow(data_links)])
target1 <- as.numeric(data_links$link2[2:nrow(data_links)])
value1 <- data_links$value[2:nrow(data_links)]

links = as.data.frame(matrix(c(source1, target1, value1), ncol = 3))
names(links) = c("source", "target", "value")

# function to make sankey plot
sankey = sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "source", Target = "target",
                     Value = "value", NodeID = "name",
                     fontSize= font_size1, sinksRight = FALSE, width = png_width, height = png_height,
                     nodePadding = nodepadding_size)

sankey$x$nodes$target <- nodes$target
sankey$x$nodes$target[1] <- FALSE

sankey <- onRender(sankey,
                 '
function(el) {
d3.select(el)
  .selectAll(".node text")
  .filter(d => d.target)
  .attr("x", -16)
  .attr("text-anchor", "end");
}
'
)

filename_html <- gsub(".csv", ".html", input_filename)
# save the widget
saveNetwork(sankey, filename_html)


