# Strain Frequency Inference

Strain frequencies are inferred from SNV frequency trajectories across samples using a method developed in Roodgar et al., 2018 and extended in Wolff et al., 2023. 

To infer strain frequencies, execute the following code in the [`scripts/strain_inference/`](https://github.com/garudlab/Wasney-Briscoe-2024/tree/main/scripts/strain_inference) directory for the species of interest:

```
python strain_inference.py --species Bacteroides_vulgatus_57955
```

Ensure that the version of python 3 (not 2.7) is loaded. If `--species` is not passed to the script, it will by default process Bacteroides_vulgatus_57955.

the `strain_inference.py` script should produce the following files:


