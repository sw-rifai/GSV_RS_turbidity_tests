library(ggplot2)

# Functions ---------------------------------------------------------------
fn_theme <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          # legend.position = c(0.5,0.025),
          # legend.justification = c(0.5,0.025),
          # legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE)
}

make_fullsize <- function() structure("", class = "fullsizebar")
ggplot_add.fullsizebar <- function(obj, g, name = "fullsizebar") {
  h <- ggplotGrob(g)$heights
  panel <- which(grid::unitType(h) == "null")
  panel_height <- unit(1, "npc") - sum(h[-panel])
  panel_height <- panel_height*0.9
  g + 
    guides(fill = guide_colorbar(barheight = panel_height,
                                 title.position = "top")) +
    theme(legend.title = element_text(angle = 0, hjust = 0.5,
                                      size = rel(2)))
}
