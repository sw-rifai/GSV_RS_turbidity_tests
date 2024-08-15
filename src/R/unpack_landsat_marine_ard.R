pacman::p_load(data.table, tidyverse, future.apply)

base <- "/media/sami/drive_a/landsat_marine/espa-srifai@gmail.com-10132023-002201-116/"
target <- "/media/sami/drive_a/landsat_marine/AMC"
fp_list <- list.files(base)


fn_unpack <- function(fp){
  dest_dir <- file.path(target,str_remove(fp,".tar.gz"))
  src_path <- file.path(base,fp)
  untar(tarfile = src_path,
        exdir = dest_dir)
  # R.utils::gunzip(filename=file.path(base,fp),
  #                 pathname = dest_dir,
  #                 remove =F)
  print(paste(fp," unpacked"))
}

plan(multisession)

future_lapply(fp_list, FUN = fn_unpack)

plan(sequential)
gc(reset=T,full=T)
