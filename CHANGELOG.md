Pre processing
# Changelog
(Based on https://wiki.canterbury.ac.nz/download/attachments/58458136/CodeVersioning_v18p2.pdf?version=1&modificationDate=1519269238437&api=v2 )

## [20.3.2]
### Changed
 - Default alpha_rough (rough) for subduction and all 5.4.2 generated srfs to be 0. Previous behaviour can be achieved using a value of -1

## [20.3.1]
### Changed
- VMs are generated with respect to rrup rather than assuming rrup==rjb
- VM generation can specify a minimum size VM to generate.
- VMs which have a SRF Ztor larger than the rrup will not be generated
- VM generation defaultly writes to VMs rather than autovm
- If the VM folder already contains files - do not overwrite them
### Added
- VM generation - Added a flag for multithreading - this passes the number of threads to the NZVM binary 

## [19.5.6] - 2019-11-26 -- Added type 2 and single realisation to new source generation workflow
### Added
- Type 2 unperturbed faults are now available in the new workflow
- Additional logging for srf generation scripts
- Global and source specific parameters are usable with the realisation generation scripts
### Changed
- If an only has one event the event will be realised without a \_RELXX suffix, working in with the old perturbation workflow
- Type is now a required parameter for generate_realisation_from_gcmt
- Some binary stderr output is now piped to the log file instead of stderr/stdout
- Only printed output should be error messages and progress indicators

## [19.5.5] - 2019-11-13 -- Added gcmt to srf generation workflow with versioned perturbation system
### Added
- Script to generate realisation files for one or more gcmt based events
- Script to generate type 1 srfs from realisation files, either one at a time, or an entire simulation

## [19.5.4] - 2019-11-06 -- Added subduction srf generation
### Added
    - Logic to generate srfs for subduction faults using genslip v5.4.2
### Changed
    - Fixed nhm 2 srf generation using 'n' randomisation method

## [19.5.3] - 2019-10-07 -- Added logging and versioned uncertainty source generation
### Added
    - Added logging for gcmt2srf, nhm2srf and srfinfo2vm workflows
    - Added versioned uncertainty source params to srf gen params script

## [19.5.2] - 2019-08-14 -- Added MIT License
### Added
    - Added MIT License and badge 

## [19.5.3] - 2019-07-12 -- Extend rrup time
### Changed
    - The rrup value is now determined by the longest distance between each VM corner and its nearest srf corner instead of the length of the sides of the VM 

## [19.5.2] - 2019-05-13 -- Truncated log normal distribution
### Added
    - Truncated log normal distribution has been added. Mean must be given in "log(x)" format 

## [19.5.1] - 2019-05-08 -- Weibull distribution
### Added
    - The Weibull distribution is now available for parameter perturbation
        - A mean value is not required for this distribution
        - The variables k and scale_factor are optional, and have default values as per nhm2srf
    - Truncated normal distribution has also been added for parameter perturbation

## [19.4.1] - 2019-04-17 -- Initial Version
### Changed
Changes to srfinfo2vm:
- if the velocity model would be outside the bounds of the current model then it is translated so that it is within the bounds.

