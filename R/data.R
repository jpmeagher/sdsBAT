#'Bat Echolocation Calls
#'
#'Bat echolocation calls were recorded across north and central Mexico from June
#'to November 2012 and from February to May 2013. Live trapped bats were
#'measured and identified to species level using field keys and recorded in
#'flight. Echolocation calls were recorded with a Pettersson 1000x bat detector
#'(Pettersson Elektronik AB, Uppsala, Sweden)  set to record calls manually in
#'realtime, full spectrum at 500 kHz. Each recording consisted of multiple calls
#'from a single individual bat. As many individual calls as possible up to a
#'maximum of 100 calls per species are included. In total the dataset consists
#'of 22 species from five families, 449 individual bats and 1816 individual
#'echolocation call recordings. See source material for species level details.
#'

#'
#'@format  Datafame of 449 observations of 36 variables: \describe{
#'  \item{species}{Key given by the first two letters of the genus and species
#'  for each bat} \item{sex}{Male/Female indicator} \item{chirp1, chirp2 , ... ,
#'  chirp34}{Individual echolocation chirp/call recordings for each bat} }
#'
#'@source V. Stathopoulos, V. Zamora-Gutierrez, K. E. Jones, and M. Girolami,
#'  Bat echolocation call identifcation for biodiversity monitoring: a
#'  probabilistic approach, Journal of the Royal Statistical Society: Series C
#'  (Applied Statistics) (2017).
#'  \url{http://www.engage-project.org/publications/}
'calls'



#' Phylogenetic Tree for selected Bat species
#'
#' Tree encoding phylogenetic relationships between the species of bat recorded
#' in the 'calls' dataset. Branch lengths set such that 1 unit = 1 million
#' years.
#'
#' @format List of type 'phylo'. Requires package 'ape'.
#'
#' @source A. Collen. The evolution of echolocation in bats: a comparative
#'   approach. PhD thesis, UCL (University College London) (2012).
'phylogeny'

#' Spectral Components
#'
#' Approximate Functional Principal Components of the Bat Echolocation Call
#' Energy Spectral Densities obtained by a resampling of call datasets balanced
#' by both species and individual and grouping components by peak frequency.
'spectral_fpca'

#'Spectral Mean of Bat Echolocation Calls
#'
#'Approximate mean Energy Spectral density for Bat Echolocation calls generated
#'by a resampling of call datasets balanced by both species and individual.
'spectral_mean'
