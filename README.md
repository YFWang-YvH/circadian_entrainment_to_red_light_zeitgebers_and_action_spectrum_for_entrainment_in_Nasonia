# circadian_entrainment_to_red_light_zeitgebers_and_action_spectrum_for_entrainment_in_Nasonia

## Project information 

This repository contains scripts used to analyze the behavioural data, as publiched in the Journal of Comparative Physiology A - Circadian entrainment to red-light Zeitgebers and action spectrum for entrainment in the jewel wasp Nasonia vitripennis
([https://doi.org/10.1098/rspb.2022.2319]).

In this project, we examed the effect of wavelength and T-cycle on entrainment and analyzed the circadian action spectrum of entrainment using generalized linear model and dose response analysis. We characterized the strength of red light as a potential Zeitgeber and established an action spectrum for light entrainment by light-dark cycles. 

Author:
Yifan Wang (YF Wang), ORCID ID: 0000-0002-6541-7435

## 1. **[lightentrainmentmanuscript.R](https://github.com/YFWang-YvH/circadian_entrainment_to_red_light_zeitgebers_and_action_spectrum_for_entrainment_in_Nasonia/blob/main/lightentrainmentmanuscript.R)**

1. Analysis of the effect of wavelengths and T-cycle on entrainment with a binomial generalized linear model (with free running animals defined as 0 and entrained animals defined as 1)
2. Analysis of the effect of wavelengths and T-cycles on phase angle of entrainment with a linear effect model and posthoc tests for pairwise comparison between wavelengths
3. Analysis of the circadian action spectrum of entrainment with a binomial generalized linear model to assess the effect of wavelength and light intensities on entrainment
4. Dose response (non-linear) analysis for constructing the circadian action spectrum and identying effective dose from the dose-response curves as an indication of Nasonia's spectral sensitivity
5. Visualization of data and model output for manuscript

---

## 2. **[supplementary.R]([https://github.com/YFWang-YvH/circadian_entrainment_of_Nasonia_by_antagonistic_interactions_of_multiple_spectral_inputs/blob/main/gam3dplot.m](https://github.com/YFWang-YvH/circadian_entrainment_to_red_light_zeitgebers_and_action_spectrum_for_entrainment_in_Nasonia/blob/main/supplementary.R))**

1. Additional visualization of periodograms, centre of gravity plots, and actograms as supplementary data to the manuscript
