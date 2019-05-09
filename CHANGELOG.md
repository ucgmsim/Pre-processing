Pre processing
# Changelog
(Based on https://wiki.canterbury.ac.nz/download/attachments/58458136/CodeVersioning_v18p2.pdf?version=1&modificationDate=1519269238437&api=v2 )

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

