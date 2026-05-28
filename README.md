# Detecting Causality with the Links--Gould Polynomial
**Authors:** Vladimir Chernov, Matthew Harper, Ben-Michael Kohli

This repository provides the official Python implementation of the paper \
[Detecting Causality with the Links--Gould Polynomial](https://arxiv.org/abs/2509.16868)

## Repository Description
We provide a script for computing the Links--Gould polynomial (LG) in the variables $$(s,q)$$ for Allen-Swenberg links.

## Scripts
For a given $$n$$, the `AS_compute.py` Python program file computes the Links--Gould polynomial of the $$n$$-th Allen-Swenberg link via the function `AS_link(n)`. The first four of these are pre-computed and stored in the `AS_link(n).pkl` files and are automatically loaded as `AS_link1`,...,`AS_link4` in the script. Additional `.pkl` files are included for the computation of higher $$n$$.

## Citation
If you use this data or script in your research, please cite our preprint.

```bibtex
@misc{LG-Allen-Swenberg,
      title = {Detecting causality with the Links--Gould polynomial}, 
      author = {Chernov, Vladimir and Harper, Matthew and Kohli, Ben-Michael},
      note = {\href{https://arxiv.org/abs/math/2605.15515}{arXiv:2605.15515}},
}
```
## Acknowledgments
BMK would like to thank Dartmouth College for its hospitality during his stay in Hanover in February 2026, when this project was initiated. MH was partially supported through the NSF-RTG grant \#DMS-2135960. The work of BMK is partially supported by the SNSF research program NCCR The Mathematics of Physics (SwissMAP), and the SNSF grant no. 200021-232258. This work was supported in part through computational resources and services provided by the Institute for Cyber-Enabled Research at Michigan State University.
