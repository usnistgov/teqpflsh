## Build the docs

In this folder, to create the runner environment, build the docs and tear down the environment

```
conda env create
conda activate teqpflshdocs
make html
conda deactivate
conda env remove -n teqpflshdocs
```