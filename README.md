# Optimization-Techniques
Simulation Project for University course "Optimization Techniques" 

---


# 1. Unidimensional Optimization Algorithms in MATLAB

[cite_start]This repository contains a MATLAB project for the "Optimization Techniques" course[cite: 2]. [cite_start]The objective is to implement, analyze, and compare four different unidimensional search algorithms for finding the minimum of three given convex functions over a specified interval[cite: 4, 5].


## 📜 Description

[cite_start]The project explores derivative-free and derivative-based methods to iteratively narrow down an interval `[a, b]` to find the minimum of a function with a specified tolerance `l`[cite: 7, 38, 44]. [cite_start]Each algorithm is applied to three distinct functions, and its performance is analyzed based on the number of function evaluations and the number of iterations required for convergence[cite: 11, 12, 14].

## ⚙️ Algorithms Implemented

Four classical optimization algorithms have been implemented:

1.  [cite_start]**Bisection Method (Derivative-Free)**: A bracketing method that iteratively selects two internal points, $x_1$ and $x_2$, near the midpoint to decide which half of the interval to discard[cite: 39, 40, 41, 42]. This implementation is found in `thema_1.m`.
2.  [cite_start]**Golden Section Search**: An efficient bracketing method that narrows the interval by a constant factor related to the **golden ratio** ($\gamma \approx 0.618$)[cite: 134]. [cite_start]Its key advantage is reusing a function evaluation from the previous iteration, which reduces the number of function calls per iteration to one[cite: 137]. This implementation is found in `thema_2.m`.
3.  [cite_start]**Fibonacci Search**: Similar to the Golden Section search, but uses ratios of consecutive **Fibonacci numbers** to determine the new test points[cite: 386]. [cite_start]This method is optimal in the sense that it requires the minimum number of function evaluations for a fixed number of iterations[cite: 388, 390]. This implementation is found in `thema_3.m`.
4.  [cite_start]**Bisection Method (with Derivatives)**: This method uses the sign of the function's first derivative at the midpoint of the interval, $x_k = (a_k + b_k)/2$[cite: 560]. [cite_start]Based on the sign, it determines if the minimum lies in the left or right half, allowing it to discard the other half[cite: 560]. This implementation is found in `thema_4.m`.

## 📈 Functions Analyzed

[cite_start]All algorithms were used to find the minimum of the following three functions over the initial interval **[-1, 3]**[cite: 4, 6]:

* $f_1(x) = (x-2)^2 + x \ln(x+3)$
* $f_2(x) = e^{-2x} + (x-2)^2$
* $f_3(x) = e^x(x^3-1) + (x-1)\sin(x)$

## 🚀 How to Run

1.  Ensure you have MATLAB installed.
2.  Clone this repository to your local machine.
3.  Open MATLAB and navigate to the repository's root directory.
4.  Run any of the `thema_*.m` scripts to execute the corresponding algorithm and generate the analysis plots.

* `thema_1.m`: Bisection Method (Derivative-Free)
* `thema_2.m`: Golden Section Search
* `thema_3.m`: Fibonacci Search
* `thema_4.m`: Bisection Method (with Derivatives)

## 🗂️ Project Structure

* `thema_1.m`: MATLAB script implementing and testing the **Bisection Method (Derivative-Free)**.
* `thema_2.m`: MATLAB script implementing and testing the **Golden Section Search**.
* `thema_3.m`: MATLAB script implementing and testing the **Fibonacci Search**.
* `thema_4.m`: MATLAB script implementing and testing the **Bisection Method with Derivatives**.
* `Work1.pdf`: The detailed project report (in Greek) that contains in-depth analysis, all generated plots, and final conclusions.

## 📊 Results & Conclusions

[cite_start]The performance of the algorithms was evaluated based on two main criteria: the total number of function/derivative evaluations (`cnt`) and the number of iterations (`k`) required to reach a certain tolerance `l`[cite: 11, 12].

### Key Findings

* [cite_start]**Function Evaluations**: The **Golden Section** and **Fibonacci** methods are significantly more efficient than the standard Bisection method, requiring fewer function calculations to achieve the same precision[cite: 796]. [cite_start]The **Bisection method with derivatives** is the most efficient overall, needing roughly half the calculations of the Golden Section/Fibonacci methods[cite: 801].

* [cite_start]**Number of Iterations**: The **Bisection methods** (both with and without derivatives) converge in fewer iterations compared to the Golden Section and Fibonacci methods for a given tolerance `l`[cite: 797].

* **Algorithm Comparison**:
    * [cite_start]The **Golden Section** and **Fibonacci** methods are nearly equivalent in performance, though the Fibonacci method saves slightly more function calculations[cite: 799].
    * [cite_start]The **Bisection method with derivatives** is the best method for minimizing the number of function calculations, but it requires the derivative to be computed beforehand[cite: 801, 802].
    * [cite_start]For functions where the derivative is unknown or expensive to compute, the **Golden Section search** offers an excellent balance of simplicity and efficiency[cite: 796, 797].
