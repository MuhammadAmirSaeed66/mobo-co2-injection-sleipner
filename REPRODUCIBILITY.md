# Reproducibility notes

## Reported analysis

The final paper experiment uses ten paired seeds. All methods use the same
20-design initialization within each seed, followed by 80 sequential simulator
evaluations. Exact three-dimensional hypervolume and grid-only IGD+ are used for
reporting.

## Restart behavior

The main script writes seed–method results and checkpoints under
`results/reported/runs/`. Existing completed files are loaded on restart.

## Archived outputs

- `results/reported/` contains the extracted scientific outputs used in the
  manuscript.
- `results/archive/final_results_complete.zip` preserves the original output
  archive, including physical-cache files.

## External simulator

The Julia simulator must be the exact version used for the reported runs. Record
its commit hash or release tag in this file after adding it:

```text
Julia simulator repository: [ADD URL]
Julia simulator commit/tag: [ADD COMMIT OR TAG]
```

## Platform information

The original run's R session information is included at:

```text
results/reported/sessionInfo.txt
```
