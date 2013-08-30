\name{population}
\alias{population}
\title{ individual based population model }
\description{
  \code{population} is an individual based model. It uses demographic parameter of a population to calculate annual effective sizes and effective population sizes (per generation). Heterozygosity over time is observed and used to calculate genetic estimates. Demographic estimates are based on a demographic formula.
}
\usage{
population(max.repeat, max.time, Na, N.loci, N.allele, n.recruit.fem, 
skew.recruit.fem, skew.recruit.mal, surv.ad)
}
\arguments{
  \item{max.repeat}{ number of repeats (must be at least 2, but at least 10 are required that results are not based on one or two examples but are means over several simulation runs) }
  \item{max.time}{ in years; for how long simulations should run (good results require several hundred years) }
  \item{Na}{ number of adults in the actual population; will be used in the simulation as well as in demographic calculations }
  \item{N.loci}{ number of loci in the simulation (at least 100 for SNPs, at least 5-10 for microsattelites, but you can also use your actual data) }
  \item{N.allele}{ number of alleles in the simulation (2 for SNPs, usually at least 20 for microsatellies, but you can also use your actual data) }
  \item{n.recruit.fem}{ number of offspring per female per year }
  \item{skew.recruit.fem}{ female variance in reproductive success (if unknown, use 0) }
  \item{skew.recruit.mal}{ male variance in reproductive success (if unknown, use 0) }
  \item{surv.ad}{ adult (female) survival rate }
}
\details{
Plese note that the function population will run simulations to calculate effective sizes and in parallel, will calculate effective sizes using the formulae described in Grimm et al. (submitted to Molecular Ecology Resources). If your population does not show variance in male or female reproductive success, you should use the results calculated (referred to as [calc]). On the contrary, the simulation will account for variance in male or female reproductive success. Under these circumstances, the calculations might be defective and simulation results should be used (referred to as [simulation]). The output will show you important parameters as well as effective sizes with Ny being the annual effective size and Ne being the effective population size (per def. i.e. per generation). For more interpretations of the results, please refer to Grimm et al. (submitted to Molecular Ecology Resources).
}
\value{
The output is a summary of the most important parameters. These parameters include the slope of heterozygosity loss per year and per generation (note the conditions of the model: no migration or mutation) calculated by means of the individual based model according to the given demographic parameters. A mean generation length is calculated. Effective sizes are indicated if they were a result of the individual based model [simulation] or calculated by a derived demographic formula [calc]. Ny refers to annual effective sizes and Ne to effective population size (which is, per def., per generation). A graph indicating heterozygosity loss over years (left column) respectively generations (right column) is provided in the graph window and refers to means in the simulation of the individual based model. The lower graphs are log transformed and a regression line is shown in red. For further definitions and interpretations and a detailed description of the example, please refer to Grimm et al. (submitted to Molecular Ecology Resources).
}
\references{
Grimm, A., B. Gruber, M. Hoehn, K. Enders, and K. Henle. submitted to Molecular Ecology Resources. Assessing effective population size in vertebrate populations under different degrees of generational overlap.
}
\author{ Annegret Grimm }
\section{Warning }{Error in rowMeans(Het_pop) : 'x' must be an array of at least two dimensions --> number of repeats must be at least 2} 
\examples{
# image you will run your population 3 times over 300 years,
# your population contains 50 adults,
# you are looking at 10 loci holding 20 alleles on average,
# each female produces 2 offspring per year and you have 
# no idea about female or male reproductive skew,
# and the adult survival rate is 0.7, then call your population as follows:
population(3,300,50,10,20,2,0,0,0.7)
}
\keyword{ effective population size }
\keyword{ annual effective size }
\keyword{ individual based model }
\keyword{ heterozygosity }
\keyword{ demography }
