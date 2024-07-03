

# PRISMA-CH4

**Code Available from 7th July**

This repository contains the code and resources for detecting and quantifying methane (CH4) emissions using PRISMA hyperspectral satellite imagery.

The methodology, originally developed for the CLEAR-UP project funded by ASI, leverages the PRISMA Hyperspectral sensor to monitor methane emissions using an enhanced matched filter technique known as the Cluster Tuned Matched Filter (CTMF). This is applied to L1C radiance images. The unit absorption spectrum, obtained using multiple MODTRAN6® simulations, serves as the target in the matched filter, enabling methane detection and quantification. Key steps include utilizing PRISMA L1C radiance images in the 2300 nm short-wave infrared methane absorption window, adapting the matched filter to account for background clutter, and applying clustering for refined background spectrum calculations. Methane concentration is estimated by linearizing the Beer-Lambert absorption law and using the unit absorption spectrum, which represents the unit methane absorption as a function of wavelength. Automation through a lookup table (LUT) based on precomputed radiances for various conditions enables real-time applicability. Results computed over 200 PRISMA acquisitions revealed methane plumes at various sites of interest in countries including Algeria, Argentina, Brazil, China, India, Mexico, Turkmenistan, and the USA. For instance, the Natural Gas compression plant in Kamyshlydzha, Turkmenistan, consistently exhibited methane plumes across all 13 images analyzed. Several large landfills were examined using available PRISMA images or through specific acquisition requests. Evidence of methane concentration enhancement was observed in the Buenos Aires landfill and the Pirana landfill in Ahmedabad.

**Citation:**

- Ferrari, A., Laneve, G., Pampanoni, V., Carvajal, A., Rossi, F. (2024, July). Monitoring Methane Emissions from Landfills Using PRISMA Imagery. In IGARSS 2024-2024 IEEE International Geoscience and Remote Sensing Symposium (soon to be published). IEEE. (GitHub code repository available here)

**Links:**

- [Docker CTMF_v4_6](https://drive.google.com/file/d/197rUulwsnqs67gWljOwPchdRw3-Yu87Y/view?usp=sharing)
- [LUT CH4](https://drive.google.com/file/d/196adGp_XCcTXAk3SRjiOnBJxUhDANNvn/view?usp=sharing)
- [DEM](https://drive.google.com/file/d/10e1VtibryVxcHT4-Gb0ryhyk17JCF04f/view?usp=sharing)
- [PRISMA Sample Images](https://drive.google.com/drive/folders/1IqJE_szLeWtHDRRdjURAl2JIBeQL_9zy?usp=sharing)

## Methodology

### In the SWIR spectrum, the absorption windows of methane are centered about 1700nm and 2300nm.

The radiance measured by a hyperspectral sensor is modeled as a superposition of signal \( b \), scaled by its strength \( \alpha \), an average background radiance \( L_{\text{mean}} \) (averaged over the entire image), and a zero-mean noise or clutter term \( \epsilon \) \[12\]\[19\].  

\[ r = \alpha b + L_{\text{mean}} + \epsilon \]                

\( \epsilon \) represents both sensor noise and scene clutter (non-desirable components). For an uncorrelated background, the optimal filter is the basic matched filter (i.e., the target itself, scaled to make the variance of the filtered image equal to one). Considering the real case of background clutter with correlation between spectral channels, the optimal filter is “matched” both to the signature of the target and the background:

**Clutter Matched Filter (CMF):**     

\[ q = \frac{C^{-1} b}{\sqrt{b^T C^{-1} b}} \]                

Here, \( q \) is normalized to ensure that in the absence of signal, the variance of the matched filter image, \( q^T r_i \) (for the i-th pixel) is one. Values larger than one indicate strong evidence of the signature's presence, quantified as a "number of sigmas", assuming that the clutter itself is Gaussian. When the clutter's covariance matrix \( C \) is accurately known, the equation defines the optimal matched filter, maximizing the signal-to-clutter ratio. However, this covariance is seldom known beforehand and is typically deduced from the data, often by averaging the outer product of the mean-subtracted radiance across all pixels.

\[ C \approx \frac{1}{N} \sum_{i=1}^N (L_i - L_{\text{mean}})(L_i - L_{\text{mean}})^T \]          

The CTMF \[12\], performs image clustering before applying the matched filter \[14\]\[15\]\[16\]\[17\]\[18\]. The average background spectrum and covariance matrix are calculated for each class. The CTMF score for the pixel (x,y) in the j-th class is given by:

\[ \alpha_j (x,y) = q_j (L_{(x,y)_j} - L_{\text{mean}_j}) \]  

\[ q_j = \frac{b_j^T C_j^{-1}}{\sqrt{b_j^T C_j^{-1} b_j}} \]            

Where the optimal filter \( q_j \) for the j-th class is expressed as \[12\], \( L_{\text{mean},j} \) is the average radiance of class j and \( C_j \) is the related covariance matrix \[18\] and \( b_j \) is the target signal. An important assumption is made approximating the off-plume radiance with the average radiance of class j, \( L_{\text{mean},j} \), to determine the expression of \( b \).

\[ L_{\text{no\_pl}} = L_{\text{mean},j} \]     

### Methane Concentration Estimation

The impact of methane enhancement is described by the Beer-Lambert absorption law; linearization of the Beer-Lambert law using a first-order Taylor series approximation allows for the formulation of a linear inversion problem to estimate methane concentration:

\[ \rho_{LM} = \frac{(L_{(x,y)_j} - L_{\text{mean}_j}) d_j^T C_j^{-1}}{d_j^T C_j^{-1} d_j} \]          

\[ d_j = -A \left(1 + \frac{1}{\cos \theta}\right) L_{\text{mean},j} \]  

Where \( \theta \) is the solar zenith angle, \( A \) is the absorption coefficient at the instrument's wavelengths, representing the specific absorption of methane. The linear estimation method allows calculating the concentration of the gas based on the difference between the observed radiance \( L_{(x_i,y_i)} \) and the mean radiance \( L_{\text{mean},j} \). This approach is consistent with the methodology outlined in previous studies \[11\] and allows calculating the concentration of the gas as the CTMF score \( \alpha_j \) divided by the optimal filter applied to the target signal \( (q_j \cdot b_j) \).

### Scene-Specific Target Spectrum Automatic Generation

The unit absorption spectrum is derived from multiple simulations with varying methane gas concentrations using radiative transfer models such as MODTRAN6®. This spectrum represents the methane absorption as a function of wavelength, normalized per unit path length and concentration. The unit for this measurement is typically (ppm·m)\(^{-1}\), indicating the inverse of the product of methane concentration (in parts per million) and the path length (in meters) over which the absorption occurs.

Subsequently, the simulated radiance spectra are convolved with the PRISMA spectral response function, and a regression analysis is carried out to obtain the unit absorption spectrum for each PRISMA spectral channel. The slope of the regression line, which best represents the relationship between the concentration-path product and the natural logarithm of the radiance, provides the unit absorption value for each specific wavelength band. This unit absorption spectrum is scaled by element-wise multiplication with the average radiance at each wavelength to generate the target spectrum, which is then used in the matched filter. To automate the generation of the unit target spectrum and make it independent of real-time MODTRAN6® runs, a lookup table (LUT) has been developed based on the methodology described by Foote et al. (2020) \[4\]. This LUT contains precomputed at-sensor radiances for a wide range of atmospheric and geometric parameters, including variations in sensor altitude, water vapor content, ground elevation, sun zenith angle, and methane concentration. For each PRISMA acquisition, the Level 2C PRISMA image is also used to read the water vapor content, and a Digital Elevation Model (DEM) is used to read the ground elevation. By interpolating within the LUT, it is possible to quickly generate a unit target spectrum that accurately matches the conditions of a given satellite pass.

![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/631d12fe-f2eb-424c-9e27-da3ace35469a)

![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/160778eb-f03e-477b-83be-0ce0200c0bbb)

![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/d471e0b0-3343-4c26-a5ef-7cb7f16bfa4d)
