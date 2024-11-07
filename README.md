Pipeline to call variants from high-coverage mitochondrial capture sequencing
1. UMI consensus calling
2. QC

To run the pipeline, clone the repo from github. In an interactive session, initialize the analysis directory with:
```
  bash [path to repo]/mitoseq/run_mitoseq --runmode=init --workdir=[directory]
```
Edit the config/samples.tsv and config/config.yaml files, including reference genome.
Perform a dryrun with:
```
  bash [path to repo]/mitoseq/run_mitoseq --runmode=dryrun --workdir=[directory]
```
If successful, launch the pipeline with:
```
  bash [path to repo]/mitoseq/run_mitoseq --runmode=run --workdir=[directory]
```

Note: Human and mouse samples must be processed in separate runs