# Guidance to reproduce SUMMIT's results

## Simulation

### Figure 2(a) and 2(b)

```
Rscript /gpfs/research/chongwu/zichenzhang/SUMMIT-test/code/Simulation.R \
--h2_e A \
--h2_p 0.2 \
--p_causal B \
--n 31684 \
--sumstats TRUE \
--gene_ENSG ENSG00000258289 \
--UKB TRUE \
--folder_output SIM \
--t1e FALSE \
--seed C \
```

| A      | B | C |
| ----- | ----- | --- |
| 0.005 | 0.01 | 1  |
| 0.005 | 0.05 | 2  |
| 0.005 | 0.10 | 3  |
| 0.005 | 0.20 | 4  |
| 0.01 | 0.01 | 5  |
| 0.01 | 0.05 | 6  |
| 0.01 | 0.10 | 7  |
| 0.01 | 0.20 | 8  |
| 0.1 | 0.01 | 9   |
| 0.1 | 0.05 | 10  |
| 0.1 | 0.10 | 11  |
| 0.1 | 0.20 | 12  |

### Figure 2(c)

```
Rscript /gpfs/research/chongwu/zichenzhang/SUMMIT-test/code/Simulation.R \
--h2_e A \
--h2_p 0.2 \
--p_causal 0.01 \
--n B \
--sumstats TRUE \
--gene_ENSG ENSG00000258289 \
--UKB TRUE \
--folder_output SIM \
--t1e FALSE \
--seed C \
```

| A      | B | C |
| ----- | ----- | -- |
| 0.005 | 300 | 13 |
| 0.005 | 600 | 14 |
| 0.005 | 3000 | 15 |
| 0.005 | 10000 | 16 |
| 0.005 | 31684 | 17 |
| 0.01 | 300 | 18 |
| 0.01 | 600 | 19 |
| 0.01 | 3000 | 20 |
| 0.01 | 10000 | 21 |
| 0.01 | 31684 | 22 |
| 0.1 | 300 | 23 |
| 0.1 | 600 | 24 |
| 0.1 | 3000 | 25 |
| 0.1 | 10000 | 26 |
| 0.1 | 31684 | 27 |

### Supplementary Figure 1

```
Rscript /gpfs/research/chongwu/zichenzhang/SUMMIT-test/code/Simulation.R \
--h2_e A \
--h2_p 0.1 \
--p_causal B \
--n 31684 \
--sumstats TRUE \
--gene_ENSG ENSG00000258289 \
--UKB TRUE \
--folder_output SIM \
--t1e FALSE \
--seed C \
```

To reproduce this figure, use exact same options for Figure 2(a) and 2(b) but change the flag ```p_causal``` to 0.1, 0.5, and 0.8.

### Supplementary Figure 2

```
Rscript /gpfs/research/chongwu/zichenzhang/SUMMIT-test/code/Simulation.R \
--h2_e 0.05 \
--h2_p 0.2 \
--p_causal 0.2 \
--n A \
--sumstats B \
--gene_ENSG ENSG00000258289 \
--UKB TRUE \
--folder_output SIM \
--t1e FALSE \
--seed C \
```

| A      | B | C |
| ----- | ----- | -- |
| 300 | FALSE | 28 |
| 600 | FALSE | 29 |
| 3000 | FALSE | 30 |
| 10000 | FALSE | 31 |
| 31684 | FALSE | 32 |
| 31684 | TRUE | 33 |


## Notes

All the above codes are for replication purposes and the users may need to change the directory and install the necessary packages to run it smoothly. If you have any questions, feel free to contact us (Zichen Zhang, [zz17@fsu.edu](mailto:zz17@fsu.edu))



## Disclaimer

The codes are provided "as is" and the author disclaims all warranties with regard to these codes including all implied warranties of merchantability and fitness. In no event shall the author be liable for any special, direct, indirect, or consequential damages or any damages whatsoever resulting from loss of use, data or profits, whether in an action of contract, negligence or other tortious action, arising out of or in connection with the use or performance of these codes. 

### Author

Zichen Zhang, and Chong Wu
