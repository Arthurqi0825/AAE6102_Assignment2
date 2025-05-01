
# AAE6102 Assignment 1

### Satellite Communication and Navigation (2024/25 Semester 2)

#### The Hong Kong Polytechnic University  
**Department of Aeronautical and Aviation Engineering** 
# Task 1

     Differential GNSS Positioning Write a short essay (500–1000 words) comparing the pros and cons for the following GNSS techniques: Differential GNSS (DGNSS), Real-Time Kinematic (RTK), Precise Point Positioning (PPP), and PPP-RTK for **smartphone navigation.


Various localization techniques have been developed to enhance positioning accuracy for Global Navigation Satellite System (GNSS), including Differential GNSS (DGNSS), Real-Time Kinematic (RTK), Precise Point Positioning (PPP), and PPP-RTK. In recent years, smartphone navigation has become an integral part of daily life, assisting people on many perspectives from pedestrian guidance to autonomous vehicle operations. Achieving high-precision positioning on smartphones, however, remains challenging due to size of smart phone and hardware limitations. This part will compare the existing navigation methods on smartphone. 

## Differential GNSS(D-GNSS)
Differential GNSS (DGNSS) enhances positioning accuracy by utilizing correction data from ground-based reference stations, which continuously monitor GNSS signals and compute discrepancies between expected and actual measurements. These corrections, including Pseudorange Correction (PRC) and Rate of Pseudorange Correction (RRC), are transmitted to receivers and applied in real time. Thanks to its relatively simple mathematical model, DGNSS is suitable for real-time applications requiring moderate accuracy without heavy computational demands.

However, DGNSS performance degrades significantly in urban environments due to signal obstruction and multipath effects caused by tall buildings. Its accuracy also declines with increasing distance from the reference station because of spatially varying atmospheric errors. Limited ground station coverage further restricts its availability. Smartphones, equipped with low-cost GNSS chipsets and embedded antennas, are especially vulnerable to these challenges, resulting in reduced positioning accuracy and reliability when using DGNSS in complex environments.
| Pros | Cons |
| :--- | :--- |
| Simple and low computational complexity. | Accuracy degrades with distance from the reference station. |
| Suitable for real-time applications. | Performance significantly worsens in urban environments due to multipath. |
| No need for expensive equipment on the receiver side. | Sparse distribution of ground stations limits availability. |
| Provides immediate correction data. | Smartphone-grade GNSS hardware further reduces performance. |


## RTK

Real-Time Kinematic (RTK) positioning significantly enhances localization accuracy by utilizing carrier-phase measurements and a double-difference technique. By differencing measurements between receivers and satellites, RTK eliminates common GNSS errors like satellite clock and atmospheric delays, achieving centimeter-level precision. Its use of carrier signals, which operate at much higher frequencies than pseudoranges, allows for much finer measurement resolution, making RTK highly effective in precision applications.

However, RTK implementation on smartphones is challenging. It requires continuous, real-time communication with a nearby reference station, which is difficult to maintain in urban environments due to signal obstructions. Accuracy also degrades with increased distance from the base station, typically beyond 10–20 kilometers. Furthermore, smartphone-grade antennas capture significantly less signal strength—only around 8% compared to professional surveying antennas—making continuous carrier-phase tracking difficult and limiting RTK’s practical use on smartphones.

| Pros | Cons |
| :--- | :--- |
| Provides centimeter-level positioning accuracy. | Requires continuous connection to a nearby base station. |
| Eliminates many common GNSS errors through double-differencing. | Accuracy degrades with distance beyond 10–20 km from base station. |
| Rapid convergence time compared to PPP. | Urban environments cause signal blockage and instability. |
| High measurement precision using carrier phase. | Smartphone antennas severely limit signal quality and RTK reliability. |

## PPP

Precise Point Positioning (PPP) is a satellite-based technique achieving high-precision localization by using precise orbit and clock data without the need for local ground stations. It enables globally consistent, centimeter-level accuracy with a single receiver, making it ideal for remote or infrastructure-limited areas where traditional ground networks are unavailable.

However, PPP has notable limitations on smartphones. Achieving high-precision positioning requires long convergence times, often several minutes to hours. Moreover, smartphone-grade GNSS hardware—characterized by low-quality antennas and oscillators—suffers from low signal-to-noise ratios and high susceptibility to multipath and cycle slips, significantly degrading PPP performance and limiting real-time usability.

| Pros | Cons |
| :--- | :--- |
| Provides global coverage without reliance on ground stations. | Extremely long convergence time (minutes to hours). |
| Enables centimeter-level accuracy theoretically with a single receiver. | Smartphone hardware limitations (poor antennas, high noise). |
| Ideal for remote or infrastructure-limited regions. | Frequent cycle slips and signal interruptions degrade performance. |
| Globally consistent correction solutions. | Not suitable for real-time smartphone navigation. |

## PPP-RTK

Precise Point Positioning–Real-Time Kinematic (PPP-RTK) combines the wide accessibility of PPP with the fast convergence and high precision of RTK. It delivers centimeter-level positioning accuracy within seconds to minutes by leveraging precise satellite data along with atmospheric and hardware bias corrections, making it highly scalable for smartphone applications without requiring dense reference station networks.

Nevertheless, smartphone hardware presents serious limitations for PPP-RTK. Studies show that embedded antennas introduce severe multipath errors, leading to large pseudorange and phase residuals compared to external antennas. Furthermore, the computational demand of real-time PPP-RTK processing strains mobile devices, making consistent high-accuracy positioning challenging.

| Pros | Cons |
| :--- | :--- |
| Centimeter-level accuracy achievable within seconds to minutes. | Smartphone embedded antennas cause severe multipath errors. |
| Scalable over large areas with minimal infrastructure. | Computational demands strain smartphone processors. |
| Combines fast convergence with global accessibility. | Signal quality on smartphones remains a bottleneck. |
| Reduces dependency on nearby ground stations. | Residual errors significantly higher compared to external antennas. |


# Task 2

- Importing data from `NavSolution_Urban.mat` in Assignment 1
- Adding information from Sky mask information to assist on the localization

<!-- <img src= https://github.com/Arthurqi0825/AAE6102_Assignment2/blob/main/Q2/skymask.jpg alt="Skymask" width=500 height=350> -->

And plotting the result:
<img src= https://github.com/Arthurqi0825/AAE6102_Assignment2/blob/main/Q2/result.jpg alt="Result" width=500 height=350>

# Task 3
**1. Weighted Least Squares (WLS) Implementation**

The least square method for the observation can be written into 
$$
    \mathbf{\hat{z}} = \mathbf{H} \mathbf{x} + \epsilon 
$$

with LS Solution:
$$
    \hat{x} = (𝐇^𝑇 𝐇)^{−1}𝐇^𝑇 \mathbf{z} 
$$

As for the weight least squras (WLS), it can be re-written to

$$
    \hat{x} = (𝐇^𝑇 \mathbf{W} 𝐇)^{−1}𝐇^𝑇  \mathbf{W} \mathbf{z} 
$$

where:
- $H$ is the design matrix constructed from satellite positions,
- $W$ is the weighting matrix (initially an identity matrix assuming equal variance),
- $z$ represents the pseudorange measurements.

To ensure numerical stability, a small regularization term is added to the matrix during inversion. The **residuals** are computed as:

$$
r = z - A \hat{x}
$$

where $z$ is the vector of observed pseudoranges and $\hat{x}$ is the estimated position and clock bias.

### Weighted RAIM Algorithm

a. Input Data

The enhanced RAIM algorithm loads the dataset `navSolution_OpenSky.mat`, which contains corrected pseudoranges and satellite positions across multiple epochs.

b. Residual Computation

After solving for the navigation solution via WLS, the residual vector is computed for each epoch to quantify the measurement error:$$r = y - A \hat{x}$$

c. Test Statistic Computation

The chi-square test statistic is computed to evaluate the overall consistency of the residuals:

$$
\chi^2 = \frac{r^T W r}{\sigma_r^2}
$$

where $\sigma_r^2$ is the estimated variance of the residuals, normalized by the degrees of freedom based on the number of satellites.

d. Detection Threshold

Fault detection is performed using a chi-square test with the threshold derived from a chi-square distribution table:
- For a false alarm probability $ P_{fa} = 10^{-2}$, using statistic table
- For a missed detection probability $P_{md} = 10^{-7}$, using 5.33 $\sigma$:


e. Fault Detection

Each epoch’s computed $\chi^2$ value is compared against the critical threshold. A fault is declared if:

$$
\chi^2 > \chi^2_{critical}
$$

Epochs with detected faults are flagged accordingly, and the results are logged for further analysis.

f. Fault Exclusion

While the algorithm detects faults, an exclusion step (removing faulty measurements and recomputing) is not explicitly implemented in the provided code. However, it can be easily extended by re-solving WLS after excluding the detected faulty measurements.

### Protection Level (PL) Computation

Based on the calculation, the final Protection Level(PL) is *85.122m*

### Visualization and Stanford Chart Analysis

a. Visualizations

Several visualizations are generated to evaluate the RAIM performance:
- 3D Satellite Constellation across epochs.
- 2D (X-Y plane) Satellite Distribution.
- Chi-square test statistics over time with critical value indication.
- Protection Level progression over time.

These plots provide insights into the satellite geometry, test statistics, and protection level behavior.

<img src= https://github.com/Arthurqi0825/AAE6102_Assignment2/blob/main/Q3/Result.jpg alt="Result" width=500 height=350>

# Task 4
## Challenges and Difficulties of Using LEO Satellites for Navigation
```

- mode: ChatGPT - 4o
- URL: https://chatgpt.com/share/680f8c47-ac90-8005-b994-6c1cca11be96
```

| **Challenge** | **Simple Explanation** |
|:--------------|:------------------------|
| Lack of Navigation-Specific Signals | LEO satellites use communication signals not optimized for positioning. |
| Timing Instability | Many LEO satellites lack stable clocks, requiring external synchronization. |
| Large Doppler Shifts | High orbital speeds cause rapid Doppler shifts, complicating signal tracking. |
| Frequent Signal Reacquisition | Satellites move quickly, needing constant reacquisition by the receiver. |
| Complex Orbit Determination and Prediction | LEO orbits change rapidly and need frequent updates to maintain accuracy. |
| Coverage and Continuity Issues | Dense satellite networks are needed to avoid positioning gaps worldwide. |
| Signal Degradation in Urban and Indoor Environments | Multipath and blockage effects worsen positioning accuracy in cities. |
| High Receiver Complexity and Processing Demands | Extracting navigation data from LEO signals requires complex receivers. |


- **Lack of Navigation-Specific Signals**  
  Most commercial LEO satellites (e.g., Starlink, OneWeb) are not designed for navigation and do not broadcast dedicated navigation signals. Receivers must extract indirect information such as Doppler shifts or signal strength variations from communication signals [1], [2]. This non-cooperative approach complicates signal processing and reduces the precision compared to purpose-built GNSS signals.

- **Timing Instability**  
  Unlike GNSS satellites, many LEO satellites lack high-stability onboard atomic clocks. As a result, their timing accuracy is poor, requiring external synchronization methods such as GNSS-disciplined oscillators or network-based corrections [1], [3]. This adds complexity to the receiver design and introduces vulnerabilities to timing errors, which directly impact positioning accuracy.

- **Large Doppler Shifts**  
  Due to their high orbital velocities (up to 7.5 km/s), LEO satellites induce large Doppler frequency shifts—often tens of kilohertz [1], [4]. Tracking such rapidly changing signals demands highly dynamic frequency tracking loops, which are challenging to implement, especially on low-power mobile devices.

- **Frequent Signal Reacquisition**  
  LEO satellites move quickly across the sky, resulting in shorter visibility windows (typically a few minutes per satellite). Receivers must continuously reacquire and track new satellites, leading to increased computational demands and potential interruptions in navigation solutions [1].

- **Complex Orbit Determination and Prediction**  
  LEO satellites are subject to strong atmospheric drag and gravitational perturbations, which cause rapid changes in their orbits [1], [5]. Maintaining accurate orbit data (ephemeris) requires frequent updates, otherwise, ephemeris errors can accumulate quickly, severely degrading navigation performance.

- **Coverage and Continuity Issues**  
  Although LEO constellations consist of many satellites, ensuring that at least four satellites are continuously visible everywhere for real-time 3D positioning requires a dense, globally distributed network [1], [6]. Without hundreds of well-orchestrated satellites, gaps in coverage or degraded geometry can occur, particularly at high latitudes or in obstructed environments.

- **Signal Degradation in Urban and Indoor Environments**  
  LEO signals, like GNSS, are susceptible to multipath effects and signal blockages caused by tall buildings and walls [2], [4]. Since Doppler-based techniques are sensitive to small phase and frequency errors, reflections and obstructions introduce significant biases and noise, deteriorating positioning accuracy.

- **High Receiver Complexity and Processing Demands**  
  To extract navigation observables from non-cooperative signals, receivers must integrate wideband front-ends, software-defined radio (SDR) architectures, and real-time signal processing algorithms [1]. These requirements increase hardware complexity, power consumption, and cost, posing challenges for practical deployment in mass-market mobile devices.
# Task 5

```
Model: ChatGPT - 4o
URL: https://chatgpt.com/share/6812f888-5130-8005-8508-0274344b572b
```
Global Navigation Satellite Systems (GNSS) are most widely recognized for their roles in positioning, navigation, and timing. However, over the past few decades, GNSS technology has also emerged as a transformative tool in the field of remote sensing. GNSS remote sensing leverages signals transmitted from navigation satellites and their interactions with the Earth's atmosphere, land, and oceans to retrieve valuable environmental information. Among the various GNSS-based remote sensing techniques, GNSS Reflectometry (GNSS-R) stands out as an innovative and rapidly growing method, offering unique capabilities for global environmental monitoring.

GNSS Reflectometry (GNSS-R) utilizes the signals emitted by GNSS satellites that are reflected off the Earth's surface. Conceptually, GNSS-R can be considered a passive bistatic radar system, where the GNSS satellite acts as the transmitter and a separate receiver captures the reflected signals. The characteristics of these reflections—such as signal strength, delay, and Doppler shift—encode information about the reflecting surface, including its roughness, moisture content, and dielectric properties. Unlike conventional active radar systems, GNSS-R does not require a dedicated transmitter, making it significantly more cost-effective and energy-efficient. This passive nature enables the deployment of GNSS-R receivers on small satellites, aircraft, and unmanned aerial vehicles, supporting broad and frequent Earth observations.

One of the fundamental tools used in GNSS-R is the Delay Doppler Map (DDM), which records the power distribution of the reflected signals across different delays and Doppler frequencies. The shape and power of the DDM are directly influenced by the surface reflectivity and roughness. For instance, calm water surfaces produce sharp, mirror-like reflections resulting in concentrated DDM peaks, while rougher surfaces like choppy oceans or vegetated land scatter the signals, broadening the DDM features. Through the analysis of DDM parameters such as Peak Power, Delay Spread, and Leading Edge Slope, scientists can infer critical environmental properties.

The applications of GNSS-R are diverse and significant. Over oceans, GNSS-R is used to estimate wind speeds by observing how wind-induced roughness affects surface reflectivity. This application is exemplified by NASA’s Cyclone Global Navigation Satellite System (CYGNSS) mission, launched in December 2016, which utilizes a constellation of eight microsatellites equipped with GNSS-R receivers to monitor tropical cyclone development. In addition to ocean wind retrieval, GNSS-R is employed to measure sea surface height, monitor floods by detecting increases in surface reflectivity, assess lake and river levels, and observe sea ice extent and type. Furthermore, over land, GNSS-R can track soil moisture levels and vegetation biomass, providing valuable inputs for hydrological studies and agricultural monitoring.

The advantages of GNSS-R over traditional remote sensing techniques are substantial. GNSS signals, operating in the L-band microwave frequencies, can penetrate clouds, rain, and vegetation, allowing for continuous, all-weather observations. The availability of multiple GNSS constellations (such as GPS, GLONASS, Galileo, and BeiDou) ensures a dense and dynamic network of signals, increasing the temporal and spatial resolution of measurements. Additionally, GNSS-R systems are lightweight, consume less power, and offer a much lower cost per unit area of observation compared to active radar missions.

Nevertheless, GNSS-R also faces several challenges. The reflected signals are generally weaker compared to those from active radars, necessitating sophisticated receiver designs and advanced signal processing algorithms to extract meaningful information. Surface properties such as roughness, vegetation cover, and dielectric constant can introduce complexities into the interpretation of the DDM, requiring accurate geophysical models and calibration strategies. Moreover, because GNSS-R is an opportunistic remote sensing technique relying on the availability of GNSS satellites, the geometry and coverage of observations can be less flexible compared to dedicated radar missions.

Despite these challenges, the progress in GNSS-R research and technology continues to enhance its reliability and application scope. With ongoing missions like CYGNSS and future planned constellations, GNSS-R is poised to become an indispensable tool in global environmental monitoring. Its ability to provide high-frequency, all-weather observations at low cost makes it particularly valuable for climate change studies, disaster monitoring, and resource management.

In conclusion, GNSS remote sensing represents a powerful extension of navigation satellite capabilities beyond traditional positioning services. Among the different techniques, GNSS Reflectometry offers a unique and effective method for observing the Earth’s surface and atmosphere. As technological advancements continue to overcome existing limitations, GNSS-R will undoubtedly play a central role in addressing the pressing environmental and societal challenges of the future.



### References for writing tasks
```
[1] T. Walter, J. Blanch, P. Enge, P. Thevenon, and M. Meurer, "PNT through LEO satellites: A survey on current status, challenges, and opportunities," *IEEE Journal of Selected Topics in Signal Processing*, vol. 16, no. 2, pp. 247–263, Mar. 2022.  
[2] A. Albrecht, A. Konovaltsev, M. Walter, T. Hien, and M. Meurer, "LEO satellite signals for Doppler positioning in urban scenarios," in *Proceedings of the ION GNSS+*, Sep. 2021, pp. 2266–2278.  
[3] M. Amrullah, et al., "A study on LEO satellite clock stability and its impact on navigation performance," *Sensors*, vol. 21, no. 10, May 2021.  
[4] R. DiEsposti, M. Navarro, and D. Margaria, "Signal tracking and positioning using low Earth orbit satellite signals of opportunity," *NAVIGATION*, vol. 67, no. 2, pp. 377–392, 2020.  
[5] Xona Space Systems, "PULSAR Navigation System Overview," Technical White Paper, 2022.
```