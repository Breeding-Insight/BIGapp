# BIGapp 1.9.0

## Improvements

* polyRAD dosage calling: the overdispersion parameter is now estimated much faster.
  The rough genotype estimates are computed once, and the optimum is found with a
  wide coarse grid (2 to 1000) followed by an integer refinement, in place of
  repeated `TestOverdispersion()` calls over slowly expanding ranges. The search is
  capped at 1000; if the optimum reaches the cap, a warning is given and 1000 is used
  (read counts close to binomial). Previously, high-depth data could make the
  search run indefinitely.
* Documented the parameters and return value of `polyRAD_dosage_call()`.

## Changes to VCF output

* `polyRAD2vcf()`: the polyRAD genotype is now always written to the `GT` FORMAT
  field, as the first field. If the input VCF already has `GT` or `UD` fields, they
  are kept and renamed to `GTO` and `UDO` (original values), in both the header and
  the FORMAT column. Previously the polyRAD genotype was written to `GTP` when a
  `GT` field existed.

## Bug fixes

* `polyRAD2vcf()`: the header declared `GT` while the data rows used `GTP` when the
  input VCF had a `GT` field.
* `polyRAD2vcf()`: input VCFs that already had a `UD` field (e.g. updog output)
  produced a duplicated `UD` key.
* `polyRAD2vcf()`: the FORMAT string is taken from each row of the input VCF instead
  of copying the first row's FORMAT to every marker.
* Fixed the Docker image build.

# BIGapp 1.8.0

* updated minor edits by @alex-sandercock in #150
* include main updates by @alex-sandercock in #152
* madc2vcf updates by @Cristianetaniguti in #154
* Merge GS-Amat support to development by @alex-sandercock in #153
* Implement UI and madc2vcf updates with GS-Amat support by @alex-sandercock in #155

# BIGapp 1.6.0

* Improved user warnings and addressed minor bug fixes

# BIGapp 1.5.1

This release includes bug fixes and new feature updates:

* Genotype files are now validated for proper formatting when uploaded
* The correct Ref and Alt bases are retained in the VCF through the genotype processing steps
* The list of removed samples is provided in the VCF Filtering tab

# BIGapp 1.3.0

Introduced several feature updates and bug fixes for BIGapp including documentation improvements, enhanced data processing and visualization options, and updated UI text.

* Updated help file documentation with clearer examples and formatting.
* Enhanced data conversion in dosage2vcf with revised input validations and updated function calls, and improved PCA plotting with added shape options and dynamic aesthetics.
* Minor UI and version updates in home, GSAcc modules, and the application descriptor.

# BIGapp 1.2.0

* When a trait file is uploaded the user is prompted to select the missing data variable, and select the column containing the sample IDs
* Updated the VCF Filtering tab to allow sample subsetting in the Advanced Options
* Added option to subset samples for the PolyRAD dosage calling method
* Added support for genotype matrix to VCF conversion
* Fixed a UI bug when extracting target SNPs from MADC file
* Added new user warnings to PCA tab

# BIGapp 1.1.0

* The ability to get target + off-target SNPs from MADC files
* Support added for PolyRAD dosage calling (this is 100x faster than Updog at the expense of slightly more memory usage)
* Support for SVG figure download
* Export VCF as bgzipped files for interoperability with other tools
* UI improvements and bug fixes

This version will require R 4.4.0

# BIGapp 1.0.0

* Bug fixes
* Updated genomic selection methods to use relationship matrices
* Improved error handling
* Updated GWAS functions
* Updated PCA figure customization
* Improved Help tab information and links within each analysis tab
* Updated UI to be more consistent throughout

# BIGapp 0.6.1

* Minor bug fixes

# BIGapp 0.6.0

* Updated Dosage Calling tab to allow user to subset genotyping data and select parental genotypes depending on Updog model
* Genomic Selection can now use pedigree information for generating the relationship matrices

# BIGapp 0.5.0

* Original code was updated into the R Golem format
* System checks added
* Improved app stability
* New user warnings
* Updates to UI
* New plot added to PCA tab
* New Genomic Selection tab

# BIGapp 0.3.0

* Updated all analyses to use VCF
* Added VCF filtering functions
* Added ability to convert from DArT Dose Report and Counts files to VCF
* Updated Genomic Prediction analysis to use 5-fold cv
* Bug fixes for PCA figures

# BIGapp 0.2.0

* Added DAPC analysis tab
* Combined DAPC and PCA under a common "Population Structure" tab
* Updated the BIG logo

# BIGapp 0.1.0

This is the initial release for the BIG app. Several analysis tabs are still in development.
