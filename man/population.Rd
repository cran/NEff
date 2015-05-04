\name{population}
\alias{population}
\title{ individual based population model }
\description{
  \code{population} calls the individual based population model. Demographic parameter of a population are used as input values to calculate annual effective sizes and effective population sizes (per generation). Heterozygosity over time is observed and used to calculate effective sizes. The population can be adapted to every life history parameter combination or genetic variation.
}
\usage{
population (max.repeat, max.time, Na, sex.ratio=FALSE, N.loci=10, 
N.allele=20, mut.rate=FALSE, one.step=0.9,n.recruit.fem, skew.recruit.fem=0, 
skew.recruit.mal=0, delay=1,surv.juv=FALSE, surv.ad.fem, 
surv.ad.mal=surv.ad.fem, const=TRUE, quiet=TRUE)
}
\arguments{
  \item{const}{ logical; if TRUE, the population size is assumed to be constant over time (default); if FALSE; the population size is changing according to respective survival and recruitment rates }
  \item{delay}{ delay until sexual maturity is reached in years (=age at sexual maturity); default: 1 (individuals read sexual maturity at the age of 1 year) }
  \item{max.repeat}{ number of repeats (at least 10 are required that results are not based on one or two examples but are means over several simulation runs) }
  \item{max.time}{ in years; for how long simulations should run (good statisitcal results require several hundred years) }
  \item{mut.rate}{ proportion of mutations occurring per generation; can be logical or numeric; default: FALSE (no mutation occurs); if set to TRUE, a standard mutation rate of 5e-3 is assumed; if numeric, any mutation rate can be declared }
  \item{Na}{ number of adults (and subadults) in the actual population; will be used as starting point of the simulation }
  \item{N.allele}{ number of possible alleles per locus in the simulation; default: 20 }
  \item{N.loci}{ number of loci in the simulation; default: 10 }
  \item{n.recruit.fem}{ number of offspring per female per year }
  \item{one.step}{ just valid if mut.rate = TRUE or numerical; defines the proportion of single-step mutations in a two-phase mutation model; default: 0.9 (0.9 singe-step, 0.1 multi-step mutations) }
  \item{quiet}{ logical; just necessary if no repeat of the simulation is used (max.repeat=1); if TRUE, calculation will be stopped, if FALSE, the warning of no repetition will be ignored and results without repetition will be presented }
  \item{sex.ratio}{ proportion of females in the population; default: 0.5 (1:1 sex ratio) }
  \item{skew.recruit.fem}{ female variance in reproductive success (i.e., variance around mean number of offspring per female per year); default: 0 (no variation between females) }
  \item{skew.recruit.mal}{ male variance in reproductive success (i.e., variance around mean number of offspring per male per year); default: 0 (no variation between males) }
  \item{surv.ad.fem}{ adult female survival rate }
  \item{surv.ad.mal}{ adult male survival rate; default: equal to adult female survival rate }
  \item{surv.juv}{ juvenile survival rate; not applicable to the default constant population size, but mandatory to non-constant population sizes }
}
\details{
Plese note that the function population() will run simulations to calculate effective sizes and in parallel, will calculate effective sizes using the formulae described in Grimm et al. (submitted to Methods in Ecology and Evolution). The demographic formula just applies to constant population sizes, no variance in reproductive success, an equal 1:1 sex ratio and no mutation. If at least one of these conditions are not met, the demographic formula should be used carefully and one should refer to calculations based on simulations. A good indicator is the logarithmic heterozygosity graph: If the regression line fits the data, the dempgraphic formula might be applicable. For further information, please refer to Grimm et al. (submitted to Method in Ecology and Evolution).
}
\value{
The output is list of the following items: a summary of the most important parameters; heterozygosity values over time; generation length calculated each year (annual mean age of reproducing females); population size over time. Time series are annual mean values over all repeats. Parameters of the summery include the slope of heterozygosity loss per year and per generation calculated by means of the individual based model according to the given demographic parameters. Further, it includes a mean generation length. Effective sizes are indicated by the way of estimation: a result of the individual based model [simulation] or calculated by a derived demographic formula [calc]. See above for restrictions of the formula. Ny refers to annual effective sizes and Ne to effective population size (which is, per def., per generation). All values are followed by a 95% confidence interval. A graph indicating population sizes and heterozygosity loss over years (left column) respectively generations (right column) is provided in the graph window and refers to annual means in the simulation of the individual based model. The lower heterozygosity graphs are log transformed and a regression line is shown in red. 
}
\references{
Grimm, A., B. Gruber, M. Hoehn, K. Enders, and K. Henle. submitted to Molecular Ecology Resources. A model-derived short-term estimation method of effective population size for species with overlapping generations.
}
\author{ Annegret Grimm }
\examples{
# image you will run your population 3 times over 50 years 
# (note that more years and more repeats will lead to more precise results),
# your population contains 20 adults with a default 1:1 sex ratio,
# each female produces 2 offspring per year without variation in reproductive success
# and the adult survival rate is 0.7, then call your population as follows:
population(max.repeat=3, max.time=50, Na=20, n.recruit.fem=2, surv.ad.fem=0.7)

# you can adapt further life history parameters
# like a sex ratio of F:M = 1:0.67 (60% females)
# a male adult survival rate of 0.6
# non-constant population size and juvenile survival rate of 0.39
# and individuals reach sexual maturity at the age of 2
population(max.repeat=3, max.time=50, Na=20, sex.ratio=0.6, n.recruit.fem=2, delay=2,
surv.ad.fem=0.7, surv.ad.mal=0.6, surv.juv=0.39, const=FALSE)

# or you can include other genetic values like
population(max.repeat=3, max.time=50, Na=20, n.recruit.fem=2, surv.ad.fem=0.7, 
N.loci=15, N.allele=25)
}
\keyword{ effective population size }
\keyword{ annual effective size }
\keyword{ individual based model }
\keyword{ heterozygosity }
\keyword{ demography }
\keyword{ life history }
\keyword{ population development }
\keyword{ population growth }
