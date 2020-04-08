# Link function for a long genome
# Provided by Henry RCircos after reporting bug
# Modified April 2016

track.num = 1
by.chromosome = FALSE

RCircos.Long.Genome.Link.Plot <- function(link.data, track.num, 
                                          by.chromosome=FALSE)
{
  #   Arguments validate
  #   =================================================================
  #
  if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }
  if(by.chromosome!=TRUE && by.chromosome!=FALSE)
  { stop("by.chromosome must be either TRUE or FALSE.\n"); }
  
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  
  #   Check chromosome names, Start, and End positions
  #   =================================================================
  #
  #link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type="link"); # commented april
  
  #   Plot position for link track.
  #   =================================================================
  #
  one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
  start <- RCircos.Par$track.in.start - (track.num-1)*one.track;
  base.positions <- RCircos.Pos*start;
  
  data.points <- matrix(rep(0, nrow(link.data)*2), ncol=2);
  for(a.link in 1:nrow(link.data))
  {
    data.points[a.link, 1] <- RCircos.Data.Point(
      link.data[a.link, 1], link.data[a.link, 2]);
    data.points[a.link, 2] <- RCircos.Data.Point(
      link.data[a.link, 4], link.data[a.link, 5]);
    
    if(data.points[a.link, 1] < 1) data.points[a.link, 1] <- 1;
    if(data.points[a.link, 2] < 1) data.points[a.link, 2] <- 1;
  }
  
  #   Get link line colors for each pair of locations
  #   ============================================================
  #
  link.colors <- RCircos.Get.Link.Colors(link.data, by.chromosome, genomic.columns = 3); # added genomic col arg 
  
  #    Draw link lines for each pair of locations
  #   ===========================================
  #
  for(a.link in 1:nrow(data.points))
  {  
    point.one <- data.points[a.link, 1];
    point.two <- data.points[a.link, 2];
    if(point.one > point.two)
    { 
      point.one <- data.points[a.link, 2];
      point.two <- data.points[a.link, 1];
    }
    
    P0 <- as.numeric(base.positions[point.one,]);
    P2 <- as.numeric(base.positions[point.two,]);
    links <- RCircos.Link.Line(P0, P2); 
    lines(links$pos.x, links$pos.y, type="l", 
          col=link.colors[a.link] ); 
  }
}