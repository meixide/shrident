# NSW application: nsw_if_hapc_K5_s20261012

Blocks assembled: RCT, CPS-2, CPS-3, PSID-2, PSID-3. Persons: 3067; rows: 3624. Cross-fitting folds: 5.

## Fits

108 fits recorded; 100 with an interior CV minimum. hapc 2.3.0. Total fitting time 65.55 h.

Fits whose CV minimum is not interior:

| block | nuisance | fold | best_lambda | cv_index | cv_verdict |
|---|---|---|---|---|---|
| 0 | m_rct | 3 | 1.388e+04 | 30 | plateau-high (null model) |
| 1 | g_trt | 0 | 6.144e-06 | 0 | boundary-low |
| 1 | g_trt | 1 | 6.144e-06 | 0 | boundary-low |
| 1 | g_trt | 2 | 6.144e-06 | 0 | boundary-low |
| 1 | g_trt | 3 | 6.144e-06 | 0 | boundary-low |
| 1 | g_trt | 4 | 6.144e-06 | 0 | boundary-low |
| 1 | g_trt | 5 | 6.144e-06 | 0 | boundary-low |
| 4 | m_ec | 3 | 6.144e-06 | 0 | boundary-low |

Set reported in the article: `onestep_xfit`. Sets for the supplement: `onestep_full`, `plugin_eif`, `plugin_known`.

## Estimators: onestep_xfit

| label | block | psi | se | theta0 | ess | max_weight | w_eb_pair | w_gls |
|---|---|---|---|---|---|---|---|---|
| DIM-RCT | RCT | 1794 | 669.3 | 4555 |  |  | 0.08015 | 0.03346 |
| DR-RCT | RCT | 1890 | 680.6 | 4459 | 226.8 | 0.01089 | 0.0798 | 0.04331 |
| DIM-CPS-2 | CPS-2 | -3822 | 604.8 | 1.017e+04 |  |  | 0.08213 | 0.4941 |
| DRpool-CPS-2 | CPS-2 | 1391 | 603.6 | 4958 | 1567 | 0.003111 | 0.08217 | 0.2087 |
| DRext-CPS-2 | CPS-2 | 668.8 | 2355 | 5680 | 5.783 | 0.3621 | 0.03215 | 0.001508 |
| DIM-CPS-3 | CPS-3 | -635 | 675.6 | 6984 |  |  | 0.07995 | -0.01976 |
| DRpool-CPS-3 | CPS-3 | 1900 | 623.5 | 4449 | 504.9 | 0.004248 | 0.08157 | 0.2193 |
| DRext-CPS-3 | CPS-3 | 499.4 | 1263 | 5850 | 12.53 | 0.2581 | 0.06001 | -0.002865 |
| DIM-PSID-2 | PSID-2 | -3647 | 908.4 | 9996 |  |  | 0.07218 | -0.003587 |
| DRpool-PSID-2 | PSID-2 | 1770 | 648.1 | 4579 | 378.6 | 0.006589 | 0.08081 | 0.3131 |
| DRext-PSID-2 | PSID-2 | 2001 | 1305 | 4348 | 19.25 | 0.1642 | 0.0586 | 0.005321 |
| DIM-PSID-3 | PSID-3 | 1070 | 894.4 | 5279 |  |  | 0.07267 | 0.06478 |
| DRpool-PSID-3 | PSID-3 | 1644 | 664.2 | 4705 | 321.9 | 0.007412 | 0.08031 | -0.3649 |
| DRext-PSID-3 | PSID-3 | 2583 | 1339 | 3766 | 23.46 | 0.08353 | 0.0575 | 0.007599 |

## Estimators: onestep_full

| label | block | psi | se | theta0 | ess | max_weight | w_eb_pair | w_gls |
|---|---|---|---|---|---|---|---|---|
| DIM-RCT | RCT | 1794 | 669.3 | 4555 |  |  | 0.07737 | 0.0001636 |
| DR-RCT | RCT | 1746 | 675.3 | 4603 | 242.7 | 0.007395 | 0.0772 | 0.04953 |
| DIM-CPS-2 | CPS-2 | -3822 | 604.8 | 1.017e+04 |  |  | 0.07916 | 0.47 |
| DRpool-CPS-2 | CPS-2 | 1382 | 602.9 | 4968 | 1642 | 0.002946 | 0.07921 | 0.2344 |
| DRext-CPS-2 | CPS-2 | 451.4 | 1761 | 5898 | 3.87 | 0.4796 | 0.04476 | -0.001275 |
| DIM-CPS-3 | CPS-3 | -635 | 675.6 | 6984 |  |  | 0.07719 | 0.04231 |
| DRpool-CPS-3 | CPS-3 | 1847 | 625.3 | 4502 | 510.9 | 0.004019 | 0.0786 | 0.01244 |
| DRext-CPS-3 | CPS-3 | 244.8 | 1207 | 6104 | 7.298 | 0.3581 | 0.06059 | 0.009115 |
| DIM-PSID-2 | PSID-2 | -3647 | 908.4 | 9996 |  |  | 0.07011 | -0.009143 |
| DRpool-PSID-2 | PSID-2 | 1501 | 636.9 | 4849 | 403.3 | 0.005223 | 0.07828 | 0.3739 |
| DRext-PSID-2 | PSID-2 | 2595 | 1020 | 3754 | 25.77 | 0.1054 | 0.06656 | 0.01312 |
| DIM-PSID-3 | PSID-3 | 1070 | 894.4 | 5279 |  |  | 0.07056 | 0.04371 |
| DRpool-PSID-3 | PSID-3 | 1618 | 646.3 | 4731 | 335.2 | 0.006622 | 0.07802 | -0.2349 |
| DRext-PSID-3 | PSID-3 | 2122 | 1150 | 4227 | 24.83 | 0.08208 | 0.0624 | -0.003431 |

## Estimators: plugin_eif

| label | block | psi | se | theta0 | ess | max_weight | w_eb_pair | w_gls | psi_published | se_published | gap |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DIM-RCT | RCT | 1794 | 669.3 | 4555 |  |  | 0.07244 | 0.0001636 | 1794 | 671 | 0 |
| IPW-NB | RCT | 1751 | 675.3 | 4598 | 242.7 | 0.007395 | 0.07226 | 0.02477 | 1751 | 679.5 | -8.645e-06 |
| OR-NB | RCT | 1765 | 675.3 | 4584 | 242.7 | 0.007395 | 0.07226 | 0.02477 | 1765 | 670.2 | 1.81e-10 |
| DIM-CPS-2 | CPS-2 | -3822 | 604.8 | 1.017e+04 |  |  | 0.07439 | 0.47 | -3822 | 606.3 | 0 |
| OM-CPS-2 | CPS-2 | 1384 | 602.9 | 4965 | 1642 | 0.002946 | 0.07445 | 0.2344 | 1384 | 619.7 | -2.574e-10 |
| IPW-CPS-2 | CPS-2 | 1451 | 1761 | 4898 | 3.87 | 0.4796 | 0.03928 | -0.001275 | 1451 | 1618 | -0.0406 |
| DIM-CPS-3 | CPS-3 | -635 | 675.6 | 6984 |  |  | 0.07225 | 0.04231 | -635 | 677.2 | 0 |
| OM-CPS-3 | CPS-3 | 1839 | 625.3 | 4510 | 510.9 | 0.004019 | 0.07378 | 0.01244 | 1839 | 637 | 3.179e-10 |
| IPW-CPS-3 | CPS-3 | 426.6 | 1207 | 5923 | 7.298 | 0.3581 | 0.05484 | 0.009115 | 426.9 | 959.1 | -0.2538 |
| DIM-PSID-2 | PSID-2 | -3647 | 908.4 | 9996 |  |  | 0.06468 | -0.009143 | -3647 | 910.5 | 4.547e-13 |
| OM-PSID-2 | PSID-2 | 1474 | 636.9 | 4875 | 403.3 | 0.005223 | 0.07343 | 0.3739 | 1474 | 671.7 | 4.434e-11 |
| IPW-PSID-2 | PSID-2 | 1506 | 1020 | 4843 | 25.77 | 0.1054 | 0.06097 | 0.01312 | 1506 | 1060 | 0.08494 |
| DIM-PSID-3 | PSID-3 | 1070 | 894.4 | 5279 |  |  | 0.06515 | 0.04371 | 1070 | 897.4 | 0 |
| OM-PSID-3 | PSID-3 | 1550 | 646.3 | 4799 | 335.2 | 0.006622 | 0.07315 | -0.2349 | 1550 | 658.9 | 4.329e-10 |
| IPW-PSID-3 | PSID-3 | 1220 | 1150 | 5129 | 24.83 | 0.08208 | 0.05668 | -0.003431 | 1220 | 1209 | -0.008736 |

## Agreement with the published table

Largest gap 0.2538 (IPW-CPS-3); tolerance 0.5. All 15 estimates agree within the tolerance.

## Estimators: plugin_known

| label | block | psi | se | theta0 | ess | max_weight | w_eb_pair | w_gls | psi_published | se_published | gap |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DIM-RCT | RCT | 1794 | 669.3 | 4555 |  |  | 0.07174 | -0.4001 | 1794 | 671 | 0 |
| IPW-NB | RCT | 1751 | 677.8 | 4598 | 242.7 | 0.007395 | 0.07149 | -0.02005 | 1751 | 679.5 | -8.645e-06 |
| OR-NB | RCT | 1765 | 667.1 | 4584 |  |  | 0.07181 | 0.5551 | 1765 | 670.2 | 1.81e-10 |
| DIM-CPS-2 | CPS-2 | -3822 | 604.8 | 1.017e+04 |  |  | 0.07363 | 0.3955 | -3822 | 606.3 | 0 |
| OM-CPS-2 | CPS-2 | 1384 | 600.9 | 4965 |  |  | 0.07375 | 0.2811 | 1384 | 619.7 | -2.574e-10 |
| IPW-CPS-2 | CPS-2 | 1451 | 1617 | 4898 | 3.87 | 0.4796 | 0.04276 | -0.001071 | 1451 | 1618 | -0.0406 |
| DIM-CPS-3 | CPS-3 | -635 | 675.6 | 6984 |  |  | 0.07155 | 0.003056 | -635 | 677.2 | 0 |
| OM-CPS-3 | CPS-3 | 1839 | 621.2 | 4510 |  |  | 0.07316 | 0.07729 | 1839 | 637 | 3.179e-10 |
| IPW-CPS-3 | CPS-3 | 426.6 | 957.5 | 5923 | 7.298 | 0.3581 | 0.06257 | 0.01988 | 426.9 | 959.1 | -0.2538 |
| DIM-PSID-2 | PSID-2 | -3647 | 908.4 | 9996 |  |  | 0.06417 | -0.03601 | -3647 | 910.5 | 4.547e-13 |
| OM-PSID-2 | PSID-2 | 1474 | 646.3 | 4875 |  |  | 0.07243 | 0.2705 | 1474 | 671.7 | 4.434e-11 |
| IPW-PSID-2 | PSID-2 | 1506 | 1058 | 4843 | 25.77 | 0.1054 | 0.05929 | 0.01482 | 1506 | 1060 | 0.08494 |
| DIM-PSID-3 | PSID-3 | 1070 | 894.4 | 5279 |  |  | 0.06463 | 0.06939 | 1070 | 897.4 | 0 |
| OM-PSID-3 | PSID-3 | 1550 | 647.8 | 4799 |  |  | 0.07239 | -0.2248 | 1550 | 658.9 | 4.329e-10 |
| IPW-PSID-3 | PSID-3 | 1220 | 1204 | 5129 | 24.83 | 0.08208 | 0.05462 | -0.004451 | 1220 | 1209 | -0.008736 |

## Pooled estimates

Columns ending in `_zero` set the dispersion to zero: they are the estimates whose sandwich interval is valid when exact identification is maintained. `psi_convex` is the convex combination (weights minimising w'(Sigma + t I)w over the simplex); `support_convex` counts its non-zero weights. `dispersion_z` is the sample variance of the estimates standardised by its mean and standard deviation under exact identification, computed from Sigma.

`se_sandwich` uses the person-level influence matrix; `se_working` is the working-independence formula of arXiv v2 and is not a valid standard error.

| set | analysis | J | tau_pair | psi_eb_pair | se_sandwich_pair | se_working_pair | ci_lo_pair | ci_hi_pair | tau2_over_mean_v | dispersion_z | psi_eb_zero | se_sandwich_zero | ci_lo_zero | ci_hi_zero | psi_convex_zero | se_convex_zero | ci_lo_convex_zero | ci_hi_convex_zero | support_convex_zero | best_single | se_best_single | psi_gls | se_gls | sigma_rank | corr_mean |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| onestep_xfit | RCT-only | 2 | 0 | 1842 | 671.4 | 477.2 | 525.6 | 3158 | 0 | -0.3587 | 1842 | 671.4 | 525.6 | 3158 | 1803 | 669.2 | 491.6 | 3115 | 2 | DIM-RCT | 669.3 | 1803 | 669.2 | 2 | 0.9796 |
| onestep_xfit | Covariate-adjusted | 10 | 0 | 1704 | 623.7 | 247.6 | 481.2 | 2926 | 0 | -0.5818 | 1704 | 623.7 | 481.2 | 2926 | 1505 | 600.6 | 327.4 | 2682 | 6 | DRpool-CPS-2 | 603.6 | 1546 | 599.4 | 10 | 0.5489 |
| onestep_xfit | Full | 14 | 1723 | 606.2 | 632.4 | 523.2 | -633.3 | 1846 | 2.661 | 5.297 | 542.8 | 607.1 | -647.1 | 1733 | -1085 | 589.6 | -2240 | 70.93 | 9 | DRpool-CPS-2 | 603.6 | -960.1 | 587.6 | 14 | 0.5585 |
| onestep_full | RCT-only | 2 | 0 | 1771 | 670.6 | 475.4 | 456.3 | 3085 | 0 | -0.5302 | 1771 | 670.6 | 456.3 | 3085 | 1792 | 669.3 | 479.8 | 3103 | 2 | DIM-RCT | 669.3 | 1792 | 669.3 | 2 | 0.9899 |
| onestep_full | Covariate-adjusted | 10 | 0 | 1637 | 619 | 240.3 | 423.3 | 2850 | 0 | -0.06999 | 1637 | 619 | 423.3 | 2850 | 1486 | 599.3 | 311.6 | 2661 | 7 | DRpool-CPS-2 | 602.9 | 1483 | 599.2 | 10 | 0.5934 |
| onestep_full | Full | 14 | 1787 | 580.4 | 625.4 | 530.9 | -645.4 | 1806 | 3.774 | 9.317 | 543.4 | 605 | -642.4 | 1729 | -1091 | 588.6 | -2245 | 62.62 | 8 | DRpool-CPS-2 | 602.9 | -1100 | 588 | 14 | 0.5918 |
| plugin_eif | RCT-only | 3 | 0 | 1770 | 671.8 | 388.7 | 453.7 | 3087 | 0 | -0.5954 | 1770 | 671.8 | 453.7 | 3087 | 1792 | 669.3 | 480.5 | 3104 | 3 | DIM-RCT | 669.3 | 1792 | 669.3 | 2 | 0.9933 |
| plugin_eif | Covariate-adjusted | 11 | 0 | 1573 | 621 | 226.4 | 356.3 | 2791 | 0 | -0.7245 | 1573 | 621 | 356.3 | 2791 | 1440 | 599.3 | 265.1 | 2614 | 7 | OM-CPS-2 | 602.9 | 1433 | 599.2 | 10 | 0.616 |
| plugin_eif | Full | 15 | 1642 | 589.4 | 621.8 | 477.2 | -629.3 | 1808 | 3.285 | 8.515 | 586.7 | 606.4 | -601.8 | 1775 | -1105 | 588.6 | -2259 | 48.31 | 8 | OM-CPS-2 | 602.9 | -1104 | 588 | 14 | 0.6062 |
| plugin_known | RCT-only | 3 | 0 | 1770 | 669.8 | 387.6 | 457.6 | 3083 | 0 | -0.6533 | 1770 | 669.8 | 457.6 | 3083 | 1765 | 667.1 | 457.3 | 3072 | 1 | OR-NB | 667.1 | 1750 | 666.5 | 3 | 0.9931 |
| plugin_known | Covariate-adjusted | 11 | 0 | 1552 | 613.1 | 224.5 | 350.5 | 2754 | 0 | -0.7305 | 1552 | 613.1 | 350.5 | 2754 | 1448 | 593.5 | 285.1 | 2612 | 7 | OM-CPS-2 | 600.9 | 1445 | 593.4 | 11 | 0.6169 |
| plugin_known | Full | 15 | 1658 | 589.2 | 619.9 | 479 | -625.8 | 1804 | 3.591 | 10.54 | 583.5 | 604.3 | -600.9 | 1768 | -614.2 | 588.2 | -1767 | 538.6 | 9 | OM-CPS-2 | 600.9 | -476.4 | 587.6 | 15 | 0.6156 |

## Prediction intervals

Split conformal over 2000 random splits (seed 42); leave-one-out for comparison.

| set | analysis | J | split_level | n_train | n_cal | split_lo_first | split_hi_first | split_lo_median | split_hi_median | split_lo_sd | split_hi_sd | loo_level | loo_lo | loo_hi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| onestep_xfit | Covariate-adjusted | 10 | 0.6667 | 5 | 5 | 1171 | 2448 | 626.7 | 2214 | 355.9 | 229.1 | 0.8182 | 565.7 | 2484 |
| onestep_xfit | Full | 14 | 0.75 | 7 | 7 | -7881 | 2040 | -3820 | 2339 | 2471 | 216.7 | 0.8667 | -6632 | 2537 |
| onestep_full | Covariate-adjusted | 10 | 0.6667 | 5 | 5 | 1263 | 2037 | 564.1 | 2534 | 356.3 | 302.3 | 0.8182 | 472.4 | 2601 |
| onestep_full | Full | 14 | 0.75 | 7 | 7 | -6607 | 2490 | -3658 | 2449 | 2045 | 289.7 | 0.8667 | -6052 | 2691 |
| plugin_eif | Covariate-adjusted | 11 | 0.7143 | 5 | 6 | 1438 | 2052 | 707.3 | 1944 | 313.3 | 94.4 | 0.8333 | 641.2 | 2036 |
| plugin_eif | Full | 15 | 0.7778 | 7 | 8 | 394.7 | 1825 | -4080 | 1908 | 2078 | 100.1 | 0.875 | -6332 | 2003 |
| plugin_known | Covariate-adjusted | 11 | 0.7143 | 5 | 6 | 1420 | 2049 | 496 | 1932 | 416.9 | 94.3 | 0.8333 | 424.7 | 2033 |
| plugin_known | Full | 15 | 0.7778 | 7 | 8 | 418.5 | 1818 | -4048 | 1899 | 1996 | 95.41 | 0.875 | -6204 | 1990 |
