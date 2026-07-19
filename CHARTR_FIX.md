# Patch: normalizePath instead of chartr

The previous `project_setup.R` used:

```r
chartr("\\\\", "/", normalizePath(...))
```

On this R installation, the `old` argument was interpreted as longer than the
single-character replacement, causing:

```text
'old' is longer than 'new'
```

The corrected code uses:

```r
normalizePath(JULIA_PROJECT, winslash = "/", mustWork = TRUE)
```

This is portable and does not require character translation.
