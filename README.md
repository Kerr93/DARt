# DARt
[![DOI](https://zenodo.org/badge/306890959.svg)](https://zenodo.org/badge/latestdoi/306890959)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Kerr93/DARt/blob/master/LICENSE)
![Safe](https://img.shields.io/badge/Stay-Safe-red?logo=data:image/svg%2bxml;base64,PHN2ZyBpZD0iTGF5ZXJfMSIgZW5hYmxlLWJhY2tncm91bmQ9Im5ldyAwIDAgNTEwIDUxMCIgaGVpZ2h0PSI1MTIiIHZpZXdCb3g9IjAgMCA1MTAgNTEwIiB3aWR0aD0iNTEyIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPjxnPjxnPjxwYXRoIGQ9Im0xNzQuNjEgMzAwYy0yMC41OCAwLTQwLjU2IDYuOTUtNTYuNjkgMTkuNzJsLTExMC4wOSA4NS43OTd2MTA0LjQ4M2g1My41MjlsNzYuNDcxLTY1aDEyNi44MnYtMTQ1eiIgZmlsbD0iI2ZmZGRjZSIvPjwvZz48cGF0aCBkPSJtNTAyLjE3IDI4NC43MmMwIDguOTUtMy42IDE3Ljg5LTEwLjc4IDI0LjQ2bC0xNDguNTYgMTM1LjgyaC03OC4xOHYtODVoNjguMThsMTE0LjM0LTEwMC4yMWMxMi44Mi0xMS4yMyAzMi4wNi0xMC45MiA0NC41LjczIDcgNi41NSAxMC41IDE1LjM4IDEwLjUgMjQuMnoiIGZpbGw9IiNmZmNjYmQiLz48cGF0aCBkPSJtMzMyLjgzIDM0OS42M3YxMC4zN2gtNjguMTh2LTYwaDE4LjU1YzI3LjQxIDAgNDkuNjMgMjIuMjIgNDkuNjMgNDkuNjN6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTM5OS44IDc3LjN2OC4wMWMwIDIwLjY1LTguMDQgNDAuMDctMjIuNjQgNTQuNjdsLTExMi41MSAxMTIuNTF2LTIyNi42NmwzLjE4LTMuMTljMTQuNi0xNC42IDM0LjAyLTIyLjY0IDU0LjY3LTIyLjY0IDQyLjYyIDAgNzcuMyAzNC42OCA3Ny4zIDc3LjN6IiBmaWxsPSIjZDAwMDUwIi8+PHBhdGggZD0ibTI2NC42NSAyNS44M3YyMjYuNjZsLTExMi41MS0xMTIuNTFjLTE0LjYtMTQuNi0yMi42NC0zNC4wMi0yMi42NC01NC42N3YtOC4wMWMwLTQyLjYyIDM0LjY4LTc3LjMgNzcuMy03Ny4zIDIwLjY1IDAgNDAuMDYgOC4wNCA1NC42NiAyMi42NHoiIGZpbGw9IiNmZjRhNGEiLz48cGF0aCBkPSJtMjEyLjgzIDM2MC4xMnYzMGg1MS44MnYtMzB6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTI2NC42NSAzNjAuMTJ2MzBoMzYuMTRsMzIuMDQtMzB6IiBmaWxsPSIjZmZiZGE5Ii8+PC9nPjwvc3ZnPg==)
[![GitHub stars](https://img.shields.io/github/stars/Kerr93/DARt.svg?style=social&label=Star&maxAge=2592000)](https://github.com/Kerr93/DARt/stargazers)

**DARt: Estimate Real-time Infection and the Time-varying Epidemiological Parameter R<sub>t**

Online platform: https://dsi-dart.shinyapps.io/dsi-covidv4/.

## Installation
```{r, eval = FALSE}
pip install DARtTool
```
    
## Quick Start
`DARtTool` is designed to be used with simple function calls, the core
functions of `DARtTool`are `DARt()` and `car_r()`. In the following section we give an overview of the simple use case for`DARt()`and`cal_r()`.

`DARt()` and `cal_r()` are the two single-call functions can be used on its own to infer the underlying infection case curve from reported cases and estimate Rt.

Firstly we need to define a DARt class, and initialize with three parameters: gt (generation time), ds(report of incubation delay) and an inputfile in csv format.

```
import DARtTool
dart = DARtTool.DARt(filename='./inputfile.csv')
```

The function cal_r() represents the core functionality of the package aiming to:

1. infer the underlying infection case curve from observations and estimate Rt.
2. provides visualisations of results.

```
dart.cal_r()
```

Estimating the underlying infection cases and Rt curve via smoothing is substantially computationally demanding than using filter but can provide reliable estimate.
    
## Input file example
    
Input file should contain date and corresponding
observations, the date should be sorted from the oldest to the latest. The first observation number should be nonzero.

See 'Example\_input.csv':


    date	    newCasesByPublishDate
    2021/1/1	53285
    2021/1/2	57725
    2021/1/3	54990
    2021/1/4	58784
    2021/1/5	60916
    2021/1/6	62322
    2021/1/7	52618
    2021/1/8	68053
    2021/1/9	59937
    2021/1/10	54940
    ...	        ...

## Citation
If you use this tool in your research, please cite our paper.
```bibtex
@misc{yang2020revealing,
      title={Revealing the Transmission Dynamics of COVID-19: A Bayesian Framework for $R_t$ Estimation}, 
      author={Xian Yang and Shuo Wang and Yuting Xing and Ling Li and Richard Yi Da Xu and Karl J. Friston and Yike Guo},
      year={2020},
      eprint={2101.01532},
      archivePrefix={arXiv},
      primaryClass={stat.AP}
}
```

## License
This source code is licensed under the [MIT](https://github.com/Kerr93/DARt/blob/master/LICENSE) license.

