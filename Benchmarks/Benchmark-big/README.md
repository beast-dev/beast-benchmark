## benchmark-big-1.xml

Source data: SARS-CoV-2_B.1.1.7_UK.fasta
Taxa: 976
Sites: 29409
Site patterns: 2078

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators: 
  - Scale on rate with weight 3.0
  - Scale on popSize with weight 3.0
  - SubtreeLeap on tree with weight 976.0
  - fixedHeightSPR on tree with weight 97.6
  - upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states

## benchmark-big-2.xml

EBOV_1610_genomes.fasta	
Taxa: 1610	
Sites: 18996	
Site patterns: 7927
Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
- Scale on rate with weight 3.0
- Scale on popSize with weight 3.0
- SubtreeLeap on tree with weight 1610.0
- fixedHeightSPR on tree with weight 161.0
- upDown on rate and popSize with weight 3.0

MCMC: 1,000,000 states, logging every 1000 states

## benchmark-big-3.xml

Source data: bdbv_drc_2026-10-16_n899 (Bundibugyo virus, DRC; sparse mutations, many ambiguous sites)
Taxa: 899
Sites: 18903
Site patterns: 5132 (min); for the MAT variant, 2405 variable sites, with sites that are constant
apart from ambiguities in at most 25% of sequences counted as constant

Model: HKY, strict clock, coalescent constant population tree prior (replacing the source's
Bayesian SkyGrid, as in the other big benchmarks; the starting tree is simulated with popSize 1.0)
Priors: kappa ~ LogNormal(1, 1.25), frequencies ~ Dirichlet(1), rate ~ CTMC reference prior,
popSize ~ Gamma(0.001, 1000)
Operators:
- Scale on kappa with weight 1.0
- DeltaExchange on frequencies with weight 1.0
- Scale on rate with weight 3.0
- upDown on node heights and rate with weight 3.0
- SubtreeLeap on tree with weight 899.0
- fixedHeightSPR on tree with weight 89.9
- Scale on popSize with weight 3.0

Variants (model label):
- min: Felsenstein likelihood (BEAGLE)
- mat: MAT likelihood with its default branch history model (integrated hits: multiple hits summed
  analytically, which gives the same posterior as the Felsenstein likelihood)

The MAT variant resamples the mutation history of every site with the SPR operators (by default with the
local block), and updates it with a Gibbs operator on the history (mutationHistoryGibbsOperator, weight
89.9, which the local block needs for a correct posterior). It needs the MAT branch of BEAST
(v10.6.0-mat2 or later).

MCMC: 2,000,000 states, logging every 2,000 states
