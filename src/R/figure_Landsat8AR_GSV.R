aq <- vect("data/aquaculture/AquacultureLeasesAndLicences_GDA2020.shp")
aq


fl <- list.files("/media/sami/drive_a/landsat_marine/AMC/")



path <- "/media/sami/drive_a/landsat_marine/AMC/LC080980842013042702T1-SC20231013063220/"

fn_plotRGB <- function(path){
  fl <- list.files(path, full.names = T)
  p_b <- fl[str_detect(fl, "AR_BAND2")]
  p_g <- fl[str_detect(fl, "AR_BAND3")]
  p_r <- fl[str_detect(fl, "AR_BAND4")]
  r <- rast(c(p_b, p_g, p_r))
  
  terra::RGB(r) <- c(3,2,1)
  plotRGB(r,stretch='lin')
}

fn_plotRGB(path)
