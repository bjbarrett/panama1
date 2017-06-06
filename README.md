# panama1
This folder includes rstan code, graphing code, data simulations, and data for Barrett, McElreath, and Perry (2017): Payoff-biased social learning underlies the diffusion of novel extractive foraging traditions in a wild primate in Proceedings of the Royal Society-B.

1) Raw data used for model can be found in "panama_data_14days.csv" file. 

2) Global model with age effects can be found in the "PN_social_global_age.stan" file. For researchers interested in a simpler learning model, this model can be simplified to develop a frequency-dependent learning, payoff-bias, or model bias only model.

3) To fit global model look at "EWA model fits.r." It also contains code to recreate several figures.

4) mono_indexing.csv contains information for recreating several figures.

5) code for figure s3 (individual predictions) is in "figure s3 individual prediction code git.r" 

6) Code for figures 2, s4, and s5 is in "Panama raw data graphs_git.r"

7) EWA_simulations_git.R contains data simulation code. 

8) PN_social_combo.stan is model for simulation.
