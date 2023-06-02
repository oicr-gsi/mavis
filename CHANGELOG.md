## 3.0.4 - 2023-06-01
- add parameter mavisQueue, modified for Ubuntu 20
## 3.0.2 - 2022-12-13
- Modified so that all modules for specific software versions are indicated in the WDL file.  Added argument to identify the genome used for analysis.   Added variables to define modules and resource files based on the provided genome version.    
## 3.0.1 - 2022-06-21
- Added tasks to do filtering of delly files.  Input struct has an optional boolean variable called doFilter.  if set to true for a delly file, then delly input will be filtered to keep only PASS calls..  Setting this flag for other SV types is currently not supported and will be ignored 
## 3.0 - 2022-05-05
- Reversion to the 1.2 WDL, running as a single task. The rewrite in version 2.0 was designed to implement parallelization for large data that runs long. There were issues in the workflow that needed correction.  The current approach, as implemented for Marathon of Hope projects will used the simple single task workflow but limit the amount of data that can be processed. 
- Donor argument changed to SampleId.  The values given are a sample indicator, NOT a donor id
- Sample id needs to be provided to the mavis tool in a sanitized form to remove reserved charactes (eg. _ ).  The workflow however also uses this in final naming of the output.  A prefix argument is now included in the runMavis task which takes the provided SampleId without sanitization, and this is used for file naming
- Outputs were collected into arrays/globs.  This didn't appear to be necessary, and has been modified now so that the output is a) a single tab delimited file with mavis calls and b) a single zip of all the drawings
## 2.0.1 - 2021-05-28
- Migrating to Vidarr
- Fixing issues with starFusion, bringing back runtime section for zip task
- [GC-8726](https://jira.oicr.on.ca/browse/GC-8726)
- [GDI-2051](https://jira.oicr.on.ca/browse/GDI-2051)
## 2.0   - 2021-02-01
- Extensive rewrite
- Break up the workflow into multiple tasks
- Run all actions as WDL tasks, instead of submitting directly to the grid engine
- Output File instead of Array[File?]
- Output .zip archive is a directory of files, not a flat collection, for safer unzipping
- Detailed checks on output in `test/calculate.sh`
## 1.2   - 2020-10-14
- Updated Arriba support
## 1.1   - 2020-10-10
- Added Arriba support, exposed `MAVIS_TIME_LIMIT` env. variable for better control of the timeout for steps in MAVIS pipeline
## 1.0.4 - 2020-05-14
- Added timeout parameter
## 1.0.3 - 2020-02-12
- [GBS-1766](https://jira.oicr.on.ca/browse/GBS-1766) - Fix issue with cromwell not localizing related files as expected
## 1.0.2 - 2020-02-03
 - Made zipped outputs optional (for cases when PASS SVs are not available the output will be nulls)
## 1.0.1 - 2019-12-20
 - Incremented version to avoid possible collision with older version
## 1.0   - 2019-12-19
 - Initial Release
