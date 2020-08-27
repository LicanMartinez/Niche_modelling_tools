#' clim_data()

#' Crops the spatial polygon of interest on global raster data;
#'   compute the 19 annual bioclimatic variables for that polygon;
#'     (intermediate files are exported to local memory here)
#'   average these variables over defined periods of time.
#'     (final files are exported to local memory here).

#' Libraries needed:
#' c(raster, dismo, rgdal, dplyr)


#' @param ex = vector of geographical limits (W, E, S, N) 
#' @param cuts = vector of years to produce the cuts for period averaging. c(1990) separates (earliest_year, 1989) and (1990, latest_year) periods
#' @param out_crop_dir = directory where to export geographically cropped bioclimatic data
#' @param out_average_dir = directory where to export averaged bioclimatic data
#' @param year = vector giving listing the years to consider (from 1961 to 2018 by default)
#' @param print_stack_per_list = (FALSE by default) if TRUE, prints final list of stacks
#' @export

clim_data<-function(data_dir = './raw climatic data/', 
                    ex, cuts, year = 1961:2018, 
                    out_crop_dir = './cropped climatic data/', 
                    out_average_dir = './averaged and cropped climatic data/', 
                    print_stack_per_list=F){
  
  ## binning years in periods for cropped files tagging
  per_delim<-c(min(year), cuts, max(year)) # defining limits and cuts
  per_names<-paste0('period', 1:(length(cuts)+1)) %>% as.character() # defining period names
  per<-cut(year, breaks = per_delim, labels = per_names, 
           include.lowest = T, rigth=T) # factor giving the corresponding period
  
  ## Geographical cropping and bioclimatic variables calculation 
  cropbox <- as(extent(ex), 'SpatialPolygons') # defining geographical limits 
  months <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  for (i in 1:length(year)) {
    
    # import and stack raw climatic data of corresponding year i
    prec <- data_dir %>% 
      paste0("wc2.1_2.5m_prec_", year[i], "-", months, ".tif") %>% 
      stack()
    tmax <- data_dir %>% 
      paste0("wc2.1_2.5m_tmax_", year[i], "-", months, ".tif") %>% 
      stack()
    tmin <- data_dir %>% 
      paste0("wc2.1_2.5m_tmin_", year[i], "-", months, ".tif") %>% 
      stack()
    
    # crop
    crpd_prec <- crop(x = prec, y = cropbox)
    crpd_tmax <- crop(x = tmax, y = cropbox)
    crpd_tmin <- crop(x = tmin, y = cropbox)
    
    # Compute bioclimatic variables
    biovars <- biovars(crpd_prec, crpd_tmin, crpd_tmax)
    layer_filenames <- paste0(out_crop_dir, "bio", seq(1,19,1), 
                              "_", year[i], '_', per[i])
    for(j in 1:19) { # export each bioclimatic raster j of the corresponding year i
      single_band <- raster(biovars, layer = j)
      writeRaster(single_band, layer_filenames[j], format = "GTiff")
      
      # print cropping status
      crop_percentage<-100*(19*i+j)/(length(year)*19)
      print(paste('Cropping ', year[i], ' bio ', j, 
                  '(', crop_percentage, '%)', sep='')) 
    }
    # deleting objects from memory
    rm(prec, tmax, tmin, crpd_prec, crpd_tmin, crpd_tmax, biovars)
  }
  
  ## Averaging by period
  all_files<-list.files(out_crop_dir, full.names = T) # create list of paths to files
  bio_names<-paste0('bio', 1:19) # names of bioclimatic layers
  
  stacks_per_list<-list() # empty list to fill
  
  for (i in 1:length(per_names)) { # Big loop goes along periods
    per.i<-per_names[i]
    files_per<-all_files[grep(per.i, all_files)] # select paths of i period files
    av_bios<-list() # to fill with 19 mean bios corresponding to the i period
    
    for (j in 1:19) { # Small loop goes along bioclimatic variables
      b.j<-bio_names[j]
      files_bio<-files_per[grep(b.j, files_per)] # select paths of j bio files
      bio.j_stack<-stack(files_bio) # import and stack raster files
      bio.j_av<-calc(bio.j_stack, fun = mean) # average them
      av_bios[[j]]<-bio.j_av # introduce to previously created list
      
      # Export
      writeRaster(bio.j_av, filename = paste0(out_average_dir, per.i, 
                                              '_', b.j), format = "GTiff")
      # print averaging status
      av_percentage<-100*(19*i+j)/(length(per_names)*19) # percentage calculation
      print(paste('Averaging ', per.i, ' ', b.j, 
                  ' (', av_percentage, '%)', sep='')) # printing
    }
    rm(bio.j_stack, bio.j_av) # Memory cleaning
    
    per.i_stack<-stack(av_bios) # stack the averaged bios of the i period
    names(per.i_stack)<-paste(bio_names, per.i, sep='_') # proper naming 
    
    stacks_per_list[[i]]<-per.i_stack # introduce to previously created list
  }
  if(print_stack_per_list){
    stacks_per_list # print final list if print_stack_per_list == TRUE
  }
}






