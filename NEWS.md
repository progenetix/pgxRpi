## Version: 0.99.0 (2023-08-25)

- Submitted to Bioconductor

## Version: 0.99.5 (2023-10-10)

- Updated data structure for frequency data, transitioning from a simple list to Bioconductor containers.
- Removed sections on survival analysis and frequency clustering analysis from the vignettes to align with the package's scope.

## Version: 0.99.7 (2023-10-20)

- Add `pgxFilter` function to expose all available filters.

## Version: 1.1.2 (2024-05-03)

- Add `segtoFreq` function to allow CNV frequency calculation from given segment data

## Version: 1.1.3 (2024-06-14)

- Add `pgxMetaplot` function to generate survival plots from metadata
- Add `num_cores` parameter for parallel query of variants

## Version: 1.1.5 (2024-07-30)

- Added `config/datatable_mappings.yaml` to define mapping rules between Beacon JSON responses and data tables.
- Modified metadata access to retrieve data directly from the Beacon API instead of using the `services/sampletable` API.
- Enabled querying of `analyses` information.
- Updated the `type` parameter in `pgxLoader` to align more closely with Beacon v2 model entities: biosamples, individuals, analyses, and g_variants.
- Added `entry_point` parameter to `pgxLoader`.
- Removed `filterLogic` parameter from `pgxLoader`.
- Optimized parallel query for variants.
- Cleaned up code and vignettes.

## Version: 1.1.6 (2024-08-05)

- Modified `extract_general_results` function to ensure it adapts correctly to arrays.
- Moved callset and cnvstats data from the "g_variant" type to "cnv_fraction" to better align with data types.
- Removed the `pgxCount` function and integrated its functionality into `pgxLoader` with the "sample_count" type, streamlining such query.


