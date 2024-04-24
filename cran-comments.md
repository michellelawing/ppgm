## Resubmission

This is a resubmission. In this version I have:

* In DESCRIPTION - removed space after 'doi:' and written title in quotes

* In ppgm.Rd, replaced T and F with TRUE and FALSE

* Added \value to ppgmMESS and plotTraitGramMultiPhylo

* Exported plotAnimatePPGMMultiPhylo

* Changed dontrun{} to donttest{} for all functions except plotAnimatedPPGM and plotAnimatedPPGMMultiPhylo, as both these functions require additional software - ImageMagick or GraphicsMagick, and for plotGeoRates and plotGeoRatesCon as these require write permission in the directory

* added path so examples for plotting functions write to tempdir()

* added on.exit code to plotAnimatedPPGMMultiPhylo.R and plotGeoRates.R to ensure does not change user's par

* changed ppgm and ppgmConsensus outputs from aicmin and aicmax to AICcmin and AICcmax for clarity to users

## R CMD check results

0 errors | 0 warnings | 0 notes
