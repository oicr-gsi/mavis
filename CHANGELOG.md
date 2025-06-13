
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.3.2] - 2025-06-13
### Fixed
- Fixed an issue where the file paths were not updated correctly due to a missing update during deployment 3.3.0.

## [3.3.1] - 2025-05-26
- Re-deployment to enable labels for optional outputs
- [GRD-948](https://jira.oicr.on.ca/browse/GRD-948)

## [3.3.0] - 2024-06-25
### Change
- Update mavis ensembl json file to Ensembl 110

## [3.2.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to medata only).

## [3.1.0] - 2023-07-18
### Added
- Add bin parameters to improve the raliability of config file creation.

## [3.0.4] - 2023-06-01
### Added
- Add parameter mavisQueue, modified for Ubuntu 20.

## [Unreleased] - 2022-12-13
### Changed
- Modified so that all modules for specific software versions are indicated in the WDL file.  

### Added
- Added argument to identify the genome used for analysis.   
- Added variables to define modules and resource files based on the provided genome version.    

## [3.0.1] - 2022-06-21
### Added
- Added tasks to do filtering of delly files.  Input struct has an optional boolean variable called doFilter.  if set to true for a delly file, then delly input will be filtered to keep only PASS calls..  Setting this flag for other SV types is currently not supported and will be ignored 

## [3.0.0] - 2022-05-05
### Changed
- Reversion to the 1.2 WDL, running as a single task. The rewrite in version 2.0 was designed to implement parallelization for large data that runs long. There were issues in the workflow that needed correction.  The current approach, as implemented for Marathon of Hope projects will used the simple single task workflow but limit the amount of data that can be processed. 
- Donor argument changed to SampleId.  The values given are a sample indicator, NOT a donor id
- Sample id needs to be provided to the mavis tool in a sanitized form to remove reserved charactes (eg. _ ).  The workflow however also uses this in final naming of the output.  A prefix argument is now included in the runMavis task which takes the provided SampleId without sanitization, and this is used for file naming
- Outputs were collected into arrays/globs.  This didn't appear to be necessary, and has been modified now so that the output is a) a single tab delimited file with mavis calls and b) a single zip of all the drawings

## [2.0.1] - 2021-05-28
### Changed
- Migrating to Vidarr

### Fixed
- Fixing issues with starFusion, bringing back runtime section for zip task
- [GC-8726](https://jira.oicr.on.ca/browse/GC-8726)
- [GDI-2051](https://jira.oicr.on.ca/browse/GDI-2051)

## [1.0.0] - 2021-06-01
### Changed
- Migrate to Vidarr

## [2.0] - 2021-02-01
### Changed
- Extensive rewrite.
- Break up the workflow into multiple tasks.
- Run all actions as WDL tasks, instead of submitting directly to the grid engine.
- Output File instead of Array[File?]
- Output .zip archive is a directory of files, not a flat collection, for safer unzipping.
- Detailed checks on output in `test/calculate.sh`

## [1.2] - 2020-10-14
### Changed
- Updated Arriba support.

## [1.1] - 2020-10-10
### Added
- Added Arriba support, exposed `MAVIS_TIME_LIMIT` env. variable for better control of the timeout for steps in MAVIS pipeline.

## [Unreleased] - 2020-05-14
### Added
- Added timeout parameter.

## [Unreleased] - 2020-02-12
### Fixed
- [GBS-1766](https://jira.oicr.on.ca/browse/GBS-1766) - Fix issue with cromwell not localizing related files as expected.

## [Unreleased] - 2020-02-03
### Changed
- Made zipped outputs optional (for cases when PASS SVs are not available the output will be nulls).

## [Unreleased] - 2019-12-20
### Changed
 - Incremented version to avoid possible collision with older version.

## [1.0] - 2019-12-19
### Added
 - Initial Release.
