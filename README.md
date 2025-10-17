# Optimization-Techniques
Simulation Project for University course "Optimization Techniques" 

---


# 1. Unidimensional Optimization Algorithms in MATLAB
The `work1` folder contains a MATLAB project with the objective to implement, analyze, and compare four different unidimensional search algorithms for finding the minimum of three given convex functions over a specified interval.



## 📜 Description

This assignment explores derivative-free and derivative-based methods to iteratively narrow down an interval `[a, b]` to find the minimum of a function with a specified tolerance `l`. Each algorithm is applied to three distinct functions, and its performance is analyzed based on the number of function evaluations and the number of iterations required for convergence.

## ⚙️ Algorithms Implemented

Four classical optimization algorithms have been implemented:

1.  **Bisection Method (Derivative-Free)**: A bracketing method that iteratively selects two internal points, $x_1$ and $x_2$, near the midpoint to decide which half of the interval to discard. This implementation is found in `thema_1.m`.
2.  **Golden Section Search**: An efficient bracketing method that narrows the interval by a constant factor related to the **golden ratio** ($\gamma \approx 0.618$). Its key advantage is reusing a function evaluation from the previous iteration, which reduces the number of function calls per iteration to one. This implementation is found in `thema_2.m`.
3.  **Fibonacci Search**: Similar to the Golden Section search, but uses ratios of consecutive **Fibonacci numbers** to determine the new test points. This method is optimal in the sense that it requires the minimum number of function evaluations for a fixed number of iterations. This implementation is found in `thema_3.m`.
4.  **Bisection Method (with Derivatives)**: This method uses the sign of the function's first derivative at the midpoint of the interval, $x_k = (a_k + b_k)/2$. Based on the sign, it determines if the minimum lies in the left or right half, allowing it to discard the other half. This implementation is found in `thema_4.m`.

## 📈 Functions Analyzed

All algorithms were used to find the minimum of the following three functions over the initial interval **[-1, 3]**:

* $f_1(x) = (x-2)^2 + x \ln(x+3)$
* $f_2(x) = e^{-2x} + (x-2)^2$
* $f_3(x) = e^x(x^3-1) + (x-1)\sin(x)$

## 🚀 How to Run

1.  Ensure you have MATLAB installed.
2.  Clone this repository to your local machine.
3.  Open MATLAB and navigate to the repository's root directory.
4.  Run any of the `thema_*.m` scripts to execute the corresponding algorithm and generate the analysis plots.

## 🗂️ Project Structure

* `thema_1.m`: MATLAB script implementing and testing the **Bisection Method (Derivative-Free)**.
* `thema_2.m`: MATLAB script implementing and testing the **Golden Section Search**.
* `thema_3.m`: MATLAB script implementing and testing the **Fibonacci Search**.
* `thema_4.m`: MATLAB script implementing and testing the **Bisection Method with Derivatives**.
* `Work1.pdf`: The detailed project report (in Greek) that contains in-depth analysis, all generated plots, and final conclusions.

## 📊 Results & Conclusions

The performance of the algorithms was evaluated based on two main criteria: the total number of function/derivative evaluations (`cnt`) and the number of iterations (`k`) required to reach a certain tolerance `l`.

### Key Findings

* **Function Evaluations**: The **Golden Section** and **Fibonacci** methods are significantly more efficient than the standard Bisection method, requiring fewer function calculations to achieve the same precision. The **Bisection method with derivatives** is the most efficient overall, needing roughly half the calculations of the Golden Section/Fibonacci methods.

* **Number of Iterations**: The **Bisection methods** (both with and without derivatives) converge in fewer iterations compared to the Golden Section and Fibonacci methods for a given tolerance `l`.

* **Algorithm Comparison**:
    * The **Golden Section** and **Fibonacci** methods are nearly equivalent in performance, though the Fibonacci method saves slightly more function calculations.
    * The **Bisection method with derivatives** is the best method for minimizing the number of function calculations, but it requires the derivative to be computed beforehand.
    * For functions where the derivative is unknown or expensive to compute, the **Golden Section search** offers an excellent balance of simplicity and efficiency.

---

# 2. Multivariate Optimization Algorithms in MATLAB

The `work2` folder contains a MATLAB project focused on implementing and comparing three gradient-based algorithms for unconstrained multivariate optimization.


## 📜 Description

This assignment examines the performance of the Steepest Descent, Newton, and Levenberg-Marquardt methods. The core of the analysis involves testing these algorithms from different starting points and using various step-size (`gamma`) selection strategies to find the minimum of a two-variable function.

## ⚙️ Algorithms Implemented

1. **Steepest Descent Method**: An iterative method that moves in the direction of the negative gradient at each step. Three strategies for selecting the step size $\gamma_k$ are explored:
    * **Constant Step Size**: A fixed value of $\gamma$ is used for all iterations. 
    * **Optimal Step Size**: At each step, $\gamma_k$ is calculated to minimize the function along the descent direction, i.e., $\min_{\gamma} f(x_k - \gamma \nabla f(x_k))$. This is solved using the bisection method from the first assignment.
    * **Armijo Rule**: An adaptive line search method that ensures sufficient decrease in the function value without being excessively small.

2. **Newton's Method**: A second-order method that uses the Hessian matrix to find the descent direction. It typically converges much faster than first-order methods but requires the Hessian to be positive definite. 

3. **Levenberg-Marquardt (L-M) Method**: A modification of Newton's method designed to handle cases where the Hessian is not positive definite. It does this by adding a regularization term $\mu I$ to the Hessian, ensuring the resulting matrix is positive definite and the step is always a descent direction.

## 📈 Function Analyzed

The algorithms were used to find the minimum of the function $f(x, y) = x^5 e^{-x^2 - y^2}$ from three different starting points: **(0, 0)**, **(-1, 1)**, and **(1, -1)**.

## 🚀 How to Run

1.  Ensure you have MATLAB installed.
2.  Clone this repository to your local machine.
3.  Open MATLAB and navigate to the `work2` folder.
4.  Run the scripts to execute the algorithms and generate the analysis plots.

## 🗂️ Project Structure

* `D_THEMA_1_2.m`: Implements and tests the **Steepest Descent Method** with all three step-size strategies.
* `D_THEMA_3.m`: Implements and tests **Newton's Method**.
* `D_THEMA_4.m`: Implements and tests the **Levenberg-Marquardt Method**.
* `ergasia_2.pdf`: The detailed project report (in Greek) containing the analysis and conclusions.

## 📊 Results & Conclusions

The analysis demonstrated the high sensitivity of optimization algorithms to the starting point and the choice of step-size strategy.

### Key Findings

* **Importance of the Starting Point**: The choice of initial point proved critical.
    * From **(0, 0)**, the gradient is zero, so none of the algorithms could make progress. 
    * From **(1, -1)**, the constant step and Armijo methods became trapped in a local minimum near $x=0$, while the optimal step-size method successfully escaped it.
    * The starting point **(-1, 1)** was closest to the global minimum, and all methods converged successfully from here.

* **Algorithm Comparison**:
    * **Newton's Method** failed for this problem because the Hessian matrix was not positive definite at the given starting points, preventing the calculation of a valid descent direction.
    * **Steepest Descent** successfully converged with appropriate step-size rules but was generally slower than Levenberg-Marquardt due to its zig-zagging behavior.
    * **Levenberg-Marquardt** was the most robust and efficient method, successfully finding the minimum with fewer iterations by overcoming the limitations of Newton's method.

* **Step-Size Strategy Comparison**:
    * The **constant step size** was the least reliable and often failed to converge.
    * The **optimal step size** strategy (minimizing $\gamma$) was very fast and effective, especially for escaping local minima.
    * The **Armijo rule** proved to be the most reliable, guaranteeing convergence and offering a good balance between speed and stability. It is particularly useful for complex or high-dimensional problems.

---

# 3. Steepest Descent with and without Projection

The `work3` folder contains a MATLAB project analyzing the Steepest Descent method for both unconstrained and constrained optimization problems.


## 📜 Description

This assignment focuses on minimizing a quadratic function, $f(x, y) = \frac{1}{3}x^2 + 3y^2$. First, the standard **Steepest Descent** method is implemented, and its convergence is analytically studied with respect to the step size $\gamma$. Second, the **Projected Steepest Descent** method is used to minimize the same function subject to box constraints.

## ⚙️ Algorithms Implemented

1.  **Steepest Descent (Unconstrained)**: The classic gradient descent algorithm where the next point is found by moving in the direction of the negative gradient. The analysis heavily focuses on the impact of a constant step size $\gamma$ on the algorithm's stability and convergence speed.

2.  **Projected Steepest Descent**: An extension of gradient descent for constrained problems. After a standard gradient step, the resulting point is "projected" back onto the feasible set if it falls outside. This ensures that all iterates remain valid. The algorithm was applied to the feasible set defined by the box constraints: $-10 \le x \le 5$ and $-8 \le y \le 12$.

## 📈 Function and Constraints

* **Function**: $f(x, y) = \frac{1}{3}x^2 + 3y^2$
* **Constraints for Projection**: A convex set $X$ where $-10 \le x \le 5$ and $-8 \le y \le 12$.

## 🚀 How to Run

1.  Ensure you have MATLAB installed.
2.  Clone this repository to your local machine.
3.  Open MATLAB and navigate to the `work3` folder.
4.  Run the scripts to see the analysis for each method.

## 🗂️ Project Structure

* `t_thema_1.m`: Implements and analyzes the **Unconstrained Steepest Descent** method for various step sizes.
* `t_thema_2.m`: Implements and analyzes the **Projected Steepest Descent** method from different starting points.
* `Texnikes_3_NT_10718.pdf`: The detailed project report (in Greek) with mathematical analysis, plots, and conclusions.

## 📊 Results & Conclusions

The project highlights the critical importance of parameter selection for gradient-based methods and the utility of projection for handling constraints.

### Key Findings

* **Unconstrained Steepest Descent & Stability**:
    * The stability of the method was analytically determined to depend on the step size $\gamma$. For this specific function, the algorithm is guaranteed to converge only if **$0 \le \gamma < 1/3$**.
    * For $\gamma$ values within the stable range (e.g., 0.1, 0.3), the algorithm converged successfully. The closer $\gamma$ was to the upper bound of $1/3$, the more oscillatory the behavior became.
    * For $\gamma$ values outside the stable range (e.g., $\gamma > 1/3$), the algorithm diverged, with the iterates going to infinity.

* **Projected Steepest Descent**:
    * The projection method successfully **prevents divergence** by ensuring all points remain within the feasible box, even with parameter choices that would cause the unconstrained method to fail.
    * However, staying within the feasible set **does not guarantee convergence** to the minimum. For certain parameters, the algorithm was observed to oscillate indefinitely between points inside the box without reaching the minimum.
    * The convergence behavior was highly dependent on the choice of parameters ($\gamma$ and $s$) and the starting point. By carefully selecting these parameters, it was possible to achieve very fast convergence. For instance, tuning the product $s_k \gamma_k$ could make a variable converge to its optimal value in a single step once inside the feasible region.
