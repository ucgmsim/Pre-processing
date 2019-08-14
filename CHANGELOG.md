Pre processing
# Changelog
(Based on https://wiki.canterbury.ac.nz/download/attachments/58458136/CodeVersioning_v18p2.pdf?version=1&modificationDate=1519269238437&api=v2 )

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

