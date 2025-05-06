# **Angular Velocity Control of Hyperboloid-Shaped Spacecraft During Mars Descent** 🚀🛰️  

## **📌 Project Goals:**  
- **Objective:** Develop **control laws** (linear & nonlinear) to manage the **angular velocity** of a **single-sheet hyperboloid-shaped spacecraft** during Mars descent by adjusting its geometry.  
- **Key Focus:** Solve **nonlinear programming (NLP)** problems to optimize hyperboloid parameters (height, radii) for **min/max moments of inertia**, enabling **fuel-free stabilization** via shape morphing. 🔄📉  
- **Benefit:** Replace traditional thrusters → **reduce spacecraft mass**, **lower mission costs**, and **increase payload capacity** for Mars missions! 💰📦  

---

## **🛠️ Skills & Tools Used:**  
- **Control Theory:** Designed **linear** (1.3) and **nonlinear** (1.5) height-control laws for hyperboloid adjustment. 📐🔧  
- **Nonlinear Optimization:** Formulated NLP problems with constraints (volume, cross-section area) for hyperboloid as **solid body** (2.12) and **surface** (2.13). 🧮⚙️  
- **Software:**  
  - **MATLAB** 🖥️: Implemented **spatial grid traversal** to solve NLP, visualize results (e.g., inertia vs. height).  
- **Numerical Methods:** Analyzed inertia/velocity dynamics for real Mars landers (*Mars Polar Lander*, *Insight*, *Mars-3*). 📊🔢  

---

## **🌟 Key Results & Innovations:**  
1. **Solid vs. Surface Hyperboloid:**  
   - **Fixed Base Area (3 constraints):**  
     - Moment of inertia (**Iz**) remained constant → **no angular velocity change** (failed stabilization). 🚫📉  
   - **Variable Base Area (2 constraints):**  
     - **Iz** decreased with height, but angular velocity peaked at **<1.1 rad/s** (too low for practical use). 📉🔴  
2. **Practical Insight:** Shape morphing alone **insufficient** for effective stabilization; hybrid approaches (e.g., thrusters + shape) needed. 🤖🔧  

---

## **💡 Why This Matters:**  
- **Fuel-Free Concept Validated:** Proved feasibility of **geometry-based control** for spacecraft dynamics. 🌍🔬  
- **Lessons for Mars Missions:** Highlighted limitations of hyperboloid morphing → guides future designs for **Mars landers/rovers**. 🛰️🔴  
- **Toolchain Ready:** MATLAB pipelines adaptable for other **asymmetric spacecraft** optimization. 🚀✨  

---

## **🔬 Technical Highlights:**  
- **Constraints:**  
  - **Volume ≥ V_min** (ensure equipment space).  
  - **Cross-section ≤ S_max** (thermal shield limits).  
- **Optimization:** Minimized **Iz** while meeting constraints via MATLAB scripts.  
- **Visualization:** Plotted **Iz, ω(t), V, S, a, c vs. time/height** for all 3 landers. 📈📉  

---

# **Final Verdict:**  
This project merges **advanced control theory, optimization, and aerospace engineering** to explore fuel-free stabilization for Mars landers. While hyperboloid morphing alone fell short, the framework sets the stage for **hybrid systems** in future missions! 🎯🚀  

**#SpaceTech #ControlTheory #MATLAB #MarsMission #AerospaceEngineering #Innovation** 🛸🔴  
