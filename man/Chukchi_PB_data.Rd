\docType{data}
\name{Chukchi_PB_data}
\alias{Chukchi_PB_data}
\title{Data for implementing analysis of 2016 Chukchi Sea polar bear aerial surveys}
\description{
  This dataset includes the following objects
  \itemize{
    \item{\code{Area_table}} A data.frame providing the flight, grid cell surveyed, area surveyed and day of survey for each grid cell surveyed in US aerial surveys.  Note that the
       grid cell here is for a smaller grid just for the US portion of the study area.  The 'Cell_lookup_table' is needed to convert these to the larger U.S.-Russia grid used for analysis.
    \item{\code{bears}} An sf spatial object giving the location and time of each ``on effort" polar bear in US survey flights
    \item{\code{bears_off}} An sf spatial object giving the location and time of each polar bear detected by a human observer but out of the range of infrared sensors in US survey flights 
    \item{\code{Cell_lookup_table}} FA vector providing the grid cell ID of the larger U.S. + Russia grid for each cell of a smaller U.S. only grid
    \item{\code{Grid_list}} A list with 55 elements, each of which is an sf object providing the survey grid and associated covariates (e.g. sea ice, distance from mainland, time from snow melt onset, etc.) for each day of the survey (i.e. the first object is the grid and covariates for day one of the survey, the second is for the second day of the survey, etc.)
    \item{\code{I_US}} A binary vector indicating whether or not a grid cell is mostly in US waters
    \item{\code{pb_dist}}  A data.frame providing location and covariates associated with each Russian polar bear detection.    This includes horizontal distance from the aircraft needed for distance samplign analysis
    \item{\code{pb_ru}} A data.frame providing effort data for Russian survey flights.  Each row corresponds to a surveyed grid cell, and includes information on the date, length of the transect through a given grid cell, as well as counts of polar bears (adults and cubs) obtained in each cell.
    \item{\code{PBSG}} An sf spatial object showing the the Polar Bear Specialist Group Chukotka population boundary
    \item{\code{Regehr_study_area}} An sf spatial object showing the boundary of the study area used in the Regehr et al. (2018) physical capture study 
    \item{\code{tracks}} An sf spatial points object, where each point (row) gives a photograph (from U.S. surveys) examined for presence or absence of polar bear tracks
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

