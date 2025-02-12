## Version: 1.2.0 (2024-10-29)

bump 1.1.11 version to 1.2.0 prior to creation of RELEASE_3_20 branch

## Version: 1.2.1 (2025-01-14)

- Removed the `pgxFilter` function and incorporated its functionality into `pgxLoader` with the "filtering_terms" type, following the Beacon v2 response mapping.
- Updated vignette filenames for improved clarity.
- Made parameter checks more flexible for Beacon queries.

## Version: 1.2.2 (2025-02-12)

- Updated the "sample_count" extraction method from "services/collations" to beacon count response, expanding counts to include all available entities (analyses, biosamples, individuals) and changing the type from "sample_count" to "counts" in `pgxLoader`.
- Enabled parallel queries across multiple resource domains (in `pgxmetaLoader`).
- Optimized code for Beacon response mapping and updated YAML mapping rules.