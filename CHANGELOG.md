## 2.0.1 - 2021-05-28
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
