theme_bd <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_classic() %+replace%    #replace elements we want to change
    
    theme(
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 14),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 12),                #font size
      
      legend.title = element_text(
        family = font, 
        size = 14), 
      
      legend.text = element_text(
        family = font, 
        size = 12
      )

    )
}
