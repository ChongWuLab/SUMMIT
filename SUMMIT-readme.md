# Guidance to reproduce SUMMIT's results

## Simulation

### Figure 2(a)

```
Rscript /gpfs/research/chongwu/zichenzhang/SUMMIT-test/code/Simulation.R \
--h2_e A \
--h2_p 0.2 \
--p_causal B\
--n 31684 \
--sumstats TRUE \
--gene_ENSG ENSG00000258289 \
--UKB TRUE \
--folder_output SIM \
--t1e FALSE \
--seed C \
```

| Syntax      | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |








## Notes

All the above codes are for replication purposes and the users may need to change the directory and install the necessary packages to run it smoothly. If you have any questions, feel free to contact us (Chong Wu, [cwu3@fsu.edu](mailto:cwu3@fsu.edu))



## Disclaimer

The codes are provided "as is" and the author disclaims all warranties with regard to these codes including all implied warranties of merchantability and fitness. In no event shall the author be liable for any special, direct, indirect, or consequential damages or any damages whatsoever resulting from loss of use, data or profits, whether in an action of contract, negligence or other tortious action, arising out of or in connection with the use or performance of these codes. 

### Author

Xinwei Ma, Jingshen Wang, and Chong Wu